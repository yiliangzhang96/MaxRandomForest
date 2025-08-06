#ifndef COUNT_H
#define COUNT_H
#include "sampling.h"
void Count(shrink_reg_data_short & sreg, int* Tmp_count, bool* Tmp_unique, int dim, bool discrete, parameters & p);
void Calculate_Weight(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count, int samplesize, int RegionID);
#endif
