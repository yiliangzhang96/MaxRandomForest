#include "stl.h"
#include "count.h"
#include "sampling.h"
using namespace std;

void Count(shrink_reg_data_short & sreg, int* Tmp_count, bool* Tmp_unique, int dim, bool discrete, parameters & para){
    double* Tmp_upbound;
    Tmp_upbound=new double[dim];
    usint_mask Tmp_reg_code;

//*** Get upperbound for each dimension of the current region (start) + by FJKu    
    for(int d=0; d< dim; d++){
        Tmp_count[d]=0;
        Tmp_reg_code.x = (sreg.reg_code[d].x << 1);
        Tmp_reg_code.mask = (sreg.reg_code[d].mask << 1) + 1;
        pair<double, double> ranges = convert_ranges(Tmp_reg_code);
        Tmp_upbound[d]=ranges.second;
        if(discrete) Tmp_unique[d] = true;
        else         Tmp_unique[d] = false;
    }
//*** Get upperbound for each dimension of the current region (end) + by FJKu    
//*** Counting and check unique (start) + by FJKu    
    for(int i=0; i<sreg.num; i++){
        for(int j=0; j< dim; j++){
            //if(sreg.regdata[i][j] <= Tmp_upbound[j]){
            if(para.data1D[para.pt_start[sreg.region_start+i]+j] <= Tmp_upbound[j]){
                ++Tmp_count[j];
            }
            //if(Tmp_unique[j] && fabs(sreg.regdata[0][j] - sreg.regdata[i][j]) > 1e-10)
            if(Tmp_unique[j] && fabs(para.data1D[para.pt_start[sreg.region_start]+j] - para.data1D[para.pt_start[sreg.region_start+i]+j]) > 1e-10)
                Tmp_unique[j] = (Tmp_unique[j] && false);
        }
    } 
//*** Counting and check unique (end) + by FJKu  
    delete [] Tmp_upbound;
}

void Calculate_Weight(vector<double> & c_original, bool* Tmp_unique, int dim, int* Tmp_count, int samplesize, int RegionID){
    for(int i=0; i < dim; ++i){
        c_original[dim*RegionID+i] = 0;
        double beta=0.5;
        if ((double) samplesize / 200.0 < 0.5) beta = max(0.1, (double) samplesize / 200.0); // heuristic pseudo ct.
        if(Tmp_unique[i])
            c_original[dim*RegionID+i] = -MAXNUM;
        else{
            c_original[dim*RegionID+i] += lgamma_c(beta + Tmp_count[i]);
            c_original[dim*RegionID+i] += lgamma_c(beta + (samplesize-Tmp_count[i]));
        }
        if(c_original[dim*RegionID+i] > 1 - MAXNUM){
            c_original[dim*RegionID+i] -= lgamma_c(beta + samplesize);
            c_original[dim*RegionID+i] += samplesize * log(2.0);
        }
    }
}
