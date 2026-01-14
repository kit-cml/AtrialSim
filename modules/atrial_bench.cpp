#include "atrial_bench.hpp"

#include <types/cellmodels/grandi_2011_atrial_with_meta.hpp>

#include <types/solvers/forwardeulersolver.hpp>

#include <functions/inputoutput.hpp>
#include <functions/helper_cvode.hpp>
#include <functions/helper_math.hpp>
#include <types/cml_consts.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <cvode/cvode.h>



int atrial_bench(const Parameter *p_param)
{
  // simulation parameters.
  // Assigned to the constant variables for
  // the sake of simplicity.
  const double cycle_length = p_param->cycle_length;
  const short number_pacing = p_param->number_pacing;
  short number_pacing_write = p_param->number_pacing_write;
  const char *user_name = p_param->user_name;
  const double time_step_min = p_param->time_step_min;
  const double time_step_max = p_param->time_step_max;
  const double writing_step = p_param->writing_step;
  const char *solver_type = p_param->solver_type;
  const double stimulus_duration = p_param->stimulus_duration;
  const double stimulus_amplitude_scale = p_param->stimulus_amplitude_scale;

  const double gks_scale = p_param->gks_scale;
  const double gkr_scale = p_param->gkr_scale;
  const double gk1_scale = p_param->gk1_scale;
  const double gto_scale = p_param->gto_scale;
  const double gna_scale = p_param->gna_scale;
  const double gnal_scale = p_param->gnal_scale;
  const double pnab_scale = p_param->pnab_scale;
  const double pcab_scale = p_param->pcab_scale;


  if(time_step_min > writing_step){
    mpi_printf(cml::commons::MASTER_NODE,"%s\n%s\n",
    "WARNING!!! The writing_step values is smaller than the timestep!",
    "Simulation still run, but the time series will use time_step_min as writing step.");
  }
  else if(number_pacing_write >= number_pacing){
    number_pacing_write = number_pacing;
    mpi_printf(cml::commons::MASTER_NODE,"%s\n%s\n",
    "WARNING!!! The number_pacing_write is larger than the number_pacing",
    "All period will be printed");
  }

  // this is the cellmodel initialization part
  Cellmodel *p_cell;
  short cell_type;
  const char *cell_model = p_param->cell_model;
  if( strstr(cell_model,"endo") != NULL ) cell_type = 0;
  else if( strstr(cell_model,"epi") != NULL ) cell_type = 1;
  else if( strstr(cell_model,"myo") != NULL ) cell_type = 2;
  else cell_type = 9;
  mpi_printf(cml::commons::MASTER_NODE,"Using %s cell model\n", cell_model);

  p_cell = new grandi_2011_atrial_with_meta();
  p_cell->initConsts();

  // apply stimulus protocol from user input
  p_cell->CONSTANTS[BCL] = cycle_length;
  p_cell->CONSTANTS[duration] = stimulus_duration;
  p_cell->CONSTANTS[amp] *= stimulus_amplitude_scale;

  // apply user input conductance scale
  p_cell->CONSTANTS[gKs_jn] *= gks_scale;
  p_cell->CONSTANTS[gKs_sl] *= gks_scale;
  p_cell->CONSTANTS[gKr] *= gkr_scale;
  p_cell->CONSTANTS[gK1] *= gk1_scale;
  p_cell->CONSTANTS[gto] *= gto_scale;
  p_cell->CONSTANTS[gNa] *= gna_scale;
  p_cell->CONSTANTS[gNaL] *= gnal_scale;
  p_cell->CONSTANTS[gNaB] *= pnab_scale;
  p_cell->CONSTANTS[gCaB] *= pcab_scale;

  // variables for I/O
  char buffer[255];
  FILE* fp_time_series;
  FILE* fp_last_states;
  FILE* fp_qnet;
  snprintf(buffer, sizeof(buffer), "%s/time_series_%s.csv", cml::commons::RESULT_FOLDER, cell_model);
  fp_time_series = fopen( buffer, "w" );
  if(fp_time_series == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr,"Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
  snprintf(buffer, sizeof(buffer), "%s/gates_%s.csv", cml::commons::RESULT_FOLDER, cell_model);
  snprintf(buffer, sizeof(buffer), "%s/last_states_%hdpaces_%s.dat", cml::commons::RESULT_FOLDER, number_pacing, cell_model);
  fp_last_states = fopen( buffer, "w" );
  if(fp_last_states == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr,"Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
  snprintf(buffer, sizeof(buffer), "%s/qnet_%s.csv", cml::commons::RESULT_FOLDER, cell_model);
  fp_qnet = fopen(buffer, "w");
  if(fp_last_states == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr,"Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
  fprintf(fp_qnet,"%s,%s\n","Pace_count","Qnet");

  fprintf(fp_time_series,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
      "Time(ms)","Vm(mV)","dVm/dt(mV/ms)","cai(nM)",
      "INa(A_per_F)","INaL_jn(nA_per_F)","ICaL(A_per_F)",
      "Ito(A_per_F)","IKr(A_per_F)","IKs(A_per_F)",
      "IK1(A_per_F)","Inet(A)","Inet_AUC(C)");

  double time_step = time_step_min;
  double tcurr = 0.;
  double tmax = number_pacing*cycle_length;
  double time_point = 25.0;
  short pace_count = 0;
  short last_print_pace = number_pacing - number_pacing_write;
  double start_time = cycle_length * last_print_pace;
  double next_output_time = start_time;
  double tprint = 0.;

  // CVode solver.
  CVodeSolverData *p_cvode;
  int cvode_retval;
  if(strncasecmp(solver_type, "Euler", sizeof(solver_type)) == 0){
    mpi_printf(0, "Using Euler Solver.\n");
  }
  else if( strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0){
    mpi_printf(0, "Using CVode Solver.\n");
    p_cvode = new CVodeSolverData();  
    init_cvode(p_cvode, p_cell, tcurr);
    set_dt_cvode(p_cvode, tcurr, time_point, cycle_length,
    time_step_min, time_step_max, &time_step);
  }
  else{
    mpi_fprintf(0, stderr, "Solver type %s is undefined! Please choose the available solver types from the manual!\n", solver_type);
    return 1;
  }

  double inet = 0.;
  double inet_auc = 0.;

  while(tcurr < tmax){ // begin of computation loop
    // compute and solving part
    // and execute the function
    // when reaching end of the cycle.
    //
    // Different solver_type has different
    // method of calling computeRates().
    if(strncasecmp(solver_type, "Euler", sizeof(solver_type)) == 0){
      p_cell->computeRates(tcurr,
          p_cell->CONSTANTS,
          p_cell->RATES,
          p_cell->STATES,
          p_cell->ALGEBRAIC);
      solveEuler(time_step_min, p_cell);
      // increase the time based on the time_step.
      tcurr += time_step_min;
    }
    else if(strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0){
      cvode_retval = solve_cvode(p_cvode, p_cell, tcurr+time_step, &tcurr);
      if( cvode_retval != 0 ){
        return 1;
      }
      set_dt_cvode(p_cvode, tcurr, time_point, cycle_length,
      time_step_min, time_step_max, &time_step);
    }

    inet = p_cell->ALGEBRAIC[INaL] + p_cell->ALGEBRAIC[ICaL] + p_cell->ALGEBRAIC[Ito] + p_cell->ALGEBRAIC[IKr] + p_cell->ALGEBRAIC[IKs] +
           p_cell->ALGEBRAIC[IK1];
    inet_auc += inet * time_step;

    if(floor(tcurr/cycle_length) != pace_count){
      end_of_cycle_funct(&pace_count, last_print_pace, &inet_auc, fp_qnet);
      //mpi_printf(0,"Entering pace %hd at %lf msec.\n", pace_count, tcurr);
    }


    // output part.
    if( tcurr >= next_output_time - cml::math::EPSILON ){
      // relative time since writing began
      tprint = next_output_time - start_time;
      snprintf(buffer, sizeof(buffer),
          "%.4lf,%.4lf,%.4lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
          p_cell->STATES[V], p_cell->RATES[V], p_cell->STATES[cai]*cml::math::MILLI_TO_NANO,
          p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL_jn],
          p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito],
          p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs],
          p_cell->ALGEBRAIC[IK1], inet, inet_auc);
      fprintf(fp_time_series, "%.4lf,%s", tprint, buffer);

      // schedule next output
      next_output_time += writing_step;
    }
   
  } // end of computation loop
 
  // writing last states to be used for
  // initial condition to the drug simulation.
  short idx;
  for(idx = 0; idx < p_cell->states_size; idx++){
    fprintf(fp_last_states, "%.16lf\n", p_cell->STATES[idx]);
  }
 
  // cleaning, cleaning....
  if( strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0 ){
    clean_cvode(p_cvode);
    delete p_cvode;
  }
  fclose(fp_last_states);
  fclose(fp_time_series);
  fclose(fp_qnet);
  delete p_cell;
  
  return 0;
}

void end_of_cycle_funct(short *pace_count, short last_print_pace, double *inet_auc, FILE *fp_qnet)
{
  *pace_count += 1;
  if( *pace_count > last_print_pace ){
    double qnet = *inet_auc / 1000.;
    fprintf(fp_qnet, "%hd,%lf\n", *pace_count, qnet);
  }
  *inet_auc = 0.;
}
