#ifndef ATRIAL_BENCH_HPP
#define ATRIAL_BENCH_CPP

#include <types/parameter.hpp>

#include <cstdio>

int atrial_bench(const Parameter *p_param);
void end_of_cycle_funct(short *pace_count, short last_print_pace, double *inet_auc, FILE *fp_qnet);

#endif
