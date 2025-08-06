

#ifndef SAMPLING_H
#define	SAMPLING_H

#include "math_utils.h"
#include "data_store.h"
#include <math.h>


static const double pi=3.141596;
static const double b = 0;
static const int MAXNUM=100000000;




struct usint_mask{      //region location, thchiu
	unsigned int x;
	unsigned int mask;                    // 00000111111, with number of '1' equals to the length of x
};


struct shrink_reg{
	vector<usint_mask> reg_code;
	int num;
};

struct shrink_reg_data{             //don't use anymore, thchiu
	vector<usint_mask> reg_code;
	vector<vector<double> > regdata;
	int num;
};

struct shrink_reg_data_short{
	vector<usint_mask> reg_code;    //region location & length, use convert_ranges() to get the boundary, thchiu
	data_store<vector<double> > regdata;
	int num;
        int region_start;       //+ by thchiu

        shrink_reg_data_short(){
            num = 0;
        }

        shrink_reg_data_short(const shrink_reg_data_short& a){
            reg_code = a.reg_code;
            //regdata = a.regdata;  //- by ed520
            num = a.num;
            region_start = a.region_start; //+ by thchiu + by ed520
        }

        void operator=(const shrink_reg_data_short& a){
            reg_code = a.reg_code;
            //regdata = a.regdata;  //- by ed520
            num = a.num;
            region_start = a.region_start; //+ by thchiu + by ed520
        }
};


struct OnePartition{        //don't use anymore, thchiu
	vector<shrink_reg> sregs;
	double w;                              // log(weight)
};

struct OnePartition_data {
    vector<shrink_reg_data_short> sregs;
//  vector<shrink_reg_data> sregs;
    double w; // log(weight)
    int PartitionID;    //+ by thchiu

    OnePartition_data(){
        w = 0.0;
    }

    // no need since everything is default -- john mu
    //OnePartition_data(const OnePartition_data& a){
    //    sregs = a.sregs;
    //    w = a.w;
    //}
};

struct HDlBPscore{
    double HD;
    double score;
};



struct mixnormalparas{
    double mu1;
    double mu2;
    double mu3;
    double mu4;
    double sd1;
    double sd2;
    double ratio;
};

// a way to compare different regions(avoid duplication)
// s1 and s2 should have the same length

struct CompairSReg {

    bool operator()(const vector<usint_mask> s1, const vector<usint_mask> s2) const {
        for (int d = 0; d < (int) s1.size(); d++) {
     //       cerr<<"s1[" << d << "].mask="<< s1[d].mask
      //              <<" s2.mask["<< d <<"]="<<s2[d].mask
     //               <<" s1["<< d <<"].x="<<s1[d].x
     //               <<" s2["<< d <<"].x="<<s2[d].x<<"\n";
            if (s1[d].mask > s2[d].mask){
                //cerr<<"false\n";
                return false;}
            if (s1[d].mask < s2[d].mask) {
                //cerr<<"true\n";
                return true;}
            if (s1[d].x > s2[d].x){
                //cerr<<"false\n";
                return false;}
            if (s1[d].x < s2[d].x) {
                //cerr<<"true\n";
                return true;}
        }
      // cerr<<"false\n";
        return false;
    }
};


struct CompairSReg2 {

    int operator()(const vector<usint_mask> s1, const vector<usint_mask> s2) const {
        for (int d = 0; d < (int) s1.size(); d++) {
            if (s1[d].mask > s2[d].mask){
                return 1;}
            if (s1[d].mask < s2[d].mask) {
                return -1;}
            if (s1[d].x > s2[d].x){
                return 1;}
            if (s1[d].x < s2[d].x) {
                return -1;}
        }
        return 0;
    }
};



inline bool myfunction3 (const shrink_reg &i, const shrink_reg &j) {
    CompairSReg mouse_compair;
    return mouse_compair(i.reg_code, j.reg_code);
}


// compare two binary partitions. the sub regions should be sorted(using CompairSReg or ~2?) before compare
struct CompairSRegs
{
	 bool  operator()(const OnePartition P1, const OnePartition P2) const
	{
             CompairSReg2 mouse_compair;

             int num1 = (int)P1.sregs.size();
             int num2 = (int)P2.sregs.size();
             if(num1< num2) return true;
             for(int d=0; d< (int)P1.sregs.size(); d++){
                int comp = mouse_compair(P1.sregs[d].reg_code, P2.sregs[d].reg_code);
		if( comp==1) return true;
                if( comp==-1) return false;
             }
             return false;
	}
};


struct CompairSRegs_data
{
	 bool  operator()(const OnePartition_data P1, const OnePartition_data P2) const
	{
             CompairSReg2 mouse_compair;

             int num1 = (int)P1.sregs.size();
             int num2 = (int)P2.sregs.size();
             if(num1< num2) return true;
             for(int d=0; d< (int)P1.sregs.size(); d++){
                int comp = mouse_compair(P1.sregs[d].reg_code, P2.sregs[d].reg_code);
		if( comp==1) return true;
                if( comp==-1) return false;
             }
             return false;
	}
};


struct parameters{
    int dim;
    int n;          //# of partition, thchiu
    vector<vector<double> > data;
    vector<vector<double> > testdata;
    int levels;     //# of cut, thchiu
    int resampling; //Flag for turning on resample(), parameter, thchiu
    int steps;      //do resample() every "steps", thchiu
    int samplesize; //# of data point, thchiu
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    double maxpercentage;   //??, thchiu
    double* data1D; //+ by thchiu
    int* pt_start; //+ by thchiu   + ed520
};

//inline bool myfunction2 (const OnePartition_data &i, const OnePartition_data &j) {
//	return (i.sumlog>j.sumlog);
//}


inline bool myfunction4 (const pair<double, double> &i, const pair<double, double> &j) {
	return (i.second>j.second);
}



inline double lprod_usint_mask(const vector<usint_mask>& reg_code) {                 // return log of the area
	double result = 0;
	for (int i = 0; i < (int)reg_code.size(); i++) {
		result -= log((double)(reg_code[i].mask+1));
	}
	return result;
}




// only for comparison for the same level onepartition_data.
//inline double lBPscore(OnePartition_data P){  //- by thchiu
inline double lBPscore(OnePartition_data& P){    //+ by thchiu
    double score=0;
    for(int i=0; i< (int)P.sregs.size(); i++){
        score += lgamma_c(0.5+P.sregs[i].num);
        score -= P.sregs[i].num * lprod_usint_mask(P.sregs[i].reg_code);

    }
  //  cout<<"lbpscore="<<score<<"\n";
    return score;


}


// mimik log(factorial(n))
inline double stirling(int n){
	if(n<=1) return 0;
	if(n<10) {
		double r=0;
		for(int i=2; i<=n; i++){
			r+= log((double)i);
		}
		return r;

	}
	double r=0.5*log((double)(2.0*pi*n))+n*log((double)n)-n+1/12/n;
	return (r);
}


inline double logsum(vector<double> la){
	int len = (int)la.size();
	if(len==0) return 0;

	double m = la[0];
	for(int i=0; i<len; i++){
		if(la[i]>m)  m=la[i];
	}
	double smallerthanm = 0;
	for(int i=0; i<len; i++){
            if(la[i]==-MAXNUM) continue;
            smallerthanm += exp(la[i]-m) ;
	}
	return (m + log(smallerthanm));
}


// from usint_mask to a pair<double, double> range with its small and large boundary

inline pair<double, double> convert_ranges(usint_mask root) {
    double up = 0.5;
    //		double low=0;
    unsigned int u = root.x;
    unsigned int v = root.mask;
    while (v > 0) {
        bool ind = (u & 1u);
        if (ind) {
            up += 0.5;
            //				low += 0.5;
        }
        //			low/=2;
        up /= 2;
        u = u >> 1;
        v = v >> 1;
    }
    up = up * 2;
    double low = up - 1.0 / (double) (root.mask + 1);
    pair<double, double> range;
    range.first = low;
    range.second = up;
    return range;
}
/*
inline void PtsInReg(shrink_reg_data &sr, vector<usint_mask> reg_code, vector<vector<double> >& data, int samplesize, int ind, bool & unique, bool keep_data, bool discrete) {
    int NumInReg = 0;
    pair<double, double> ranges = convert_ranges(reg_code[ind]);
    if(discrete) unique = true;
    else unique = false;               // if there is no need to check
    if (samplesize >= 1) {
        double uniq = data[0][ind];                   // check unique
        for (int i = 0; i < samplesize; i++) {
            if (unique && fabs(uniq - data[i][ind]) > 1e-10) unique = false;  // check unique
            if (data[i][ind] <= ranges.second && data[i][ind] > ranges.first) {
                NumInReg += 1;
                if(keep_data)  sr.regdata.push_back(data[i]);
            }
        }
    }
    sr.num = NumInReg;
    sr.reg_code = reg_code;
}
*/

inline void PtsInReg_short(shrink_reg_data_short &sr, vector<usint_mask> reg_code, data_store<vector<double> >& data, int samplesize, int ind, bool & unique, bool keep_data, bool discrete) {
    int NumInReg = 0;
    pair<double, double> ranges = convert_ranges(reg_code[ind]);
    if(discrete) unique = true;
    else unique = false;               // if there is no need to check
    if (samplesize >= 1) {
        double uniq = data[0][ind];                   // check unique
        for (int i = 0; i < samplesize; i++) {
            if (unique && fabs(uniq - data[i][ind]) > 1e-10) unique = false;  // check unique
            if (data[i][ind] <= ranges.second && data[i][ind] > ranges.first) {
                NumInReg += 1;
                //if(keep_data)  sr.regdata.push_back(data.get_ptr(i)); //- by thchiu
                if(keep_data)  sr.regdata.push_back(data[i]);   //+ by thchiu
            }
        }
    }
    sr.num = NumInReg;
    sr.reg_code = reg_code;
}


inline vector<double> get_all_density (OnePartition_data& P){
	double tot = 0;
	vector<double> den;
	for(int i=0; i<(int)P.sregs.size(); i++){
		tot += P.sregs[i].num;
	}
	for(int i=0; i<(int)P.sregs.size(); i++){
		den.push_back( P.sregs[i].num / exp(lprod_usint_mask(P.sregs[i].reg_code)) /tot );
	}

	return den;
}

inline int find(OnePartition_data& P, vector<double>& sample){
	int dim = P.sregs[0].reg_code.size();
	for(int i=0; i<(int)P.sregs.size(); i++){
		bool is_in=true;

		for(int d=0; is_in && d<dim; d++){
			pair<double,double> range = convert_ranges(P.sregs[i].reg_code[d]);
			if(sample[d]>range.second || sample[d]<range.first) {
				is_in=false;
			}
		}
		if(is_in)  return i;
	}
	cout<<"wrong in finding this: "<<endl;
        for(int i=0; i<(int)sample.size(); i++)  cout<<sample[i]<<" ";
        cout<<endl;
        return 0;
}

// lattic<=0: return the density on the point sample
// lattic>0: return the mean of the densities of the points which is 1unit(length=lattic) far from sample coordinately
inline double one_density(OnePartition_data& P, vector<double>& sample, double lattic){
	vector<double> dens = get_all_density (P);
        if(lattic<=0.000001){
                int ind = find(P,sample);
                return dens[ind];
        } else{
            int dim = (int)sample.size();
            double meanden=0;                     //;
            int times=0;
            vector<double> candidate= sample;
            for(int i=0; i<dim; i++){                    // build the grid for averaging the densities
                candidate[i]=sample[i]-lattic;
                if(candidate[i]>0){
                    times++;
                    meanden += dens[find(P, candidate)];
                }
                candidate[i]=sample[i]+lattic;
                if(candidate[i]<1){
                    times++;
                    meanden += dens[find(P, candidate)];
                }
                candidate[i]=sample[i];
                if(meanden != meanden) cout<<"i="<<i<<" is.nan in one_density!";

            }
            if(times>0){
                meanden = ( meanden + (double)times*dens[find(P, sample)]) ;
                meanden = meanden /(double)times/2.0;
                return meanden;
            }
            return dens[find(P, sample)];

        }
}
//+by cy
struct region{
	double* bound;
	int time;
	
};	



//void generate_map (shrink_reg sr, vector<double> mid,vector<double> half, int dim, map<vector<usint_mask>, shrink_reg, CompairSReg>& newmap, vector<vector<double> >& data, int samplesizee, int masklen);

void get_info_from_onetree(OnePartition_data &P, vector<double> & c_original,vector<double> &c,int cutreg, double maxpercentage, double &smooth, double *time_array, bool discrete, parameters &para);


//void rand_choose_one_nextleveltree(OnePartition_data &newP, OnePartition_data &P, vector<double> &c, int &cutreg, double smooth, bool discrete);  //- by thchiu
void rand_choose_one_nextleveltree(OnePartition_data &P, vector<double> &c, int &cutreg, double smooth, bool discrete, parameters &para);  //+ by thchiu


bool sample_trees(vector<OnePartition_data>& Ps, vector<vector<double> >& c_originals, vector<int>&cutregs,  double maxpercentage, double* time_array, bool discrete, parameters& para);

void getmax_corr(vector<OnePartition_data>& newPs,  vector<int> & modeinds, vector<double>& modewts, double & ess);


double getmaxlBP(vector<OnePartition_data>& newPs, int & modelBP);


// Andy 
//double getmaxlBP(vector<Partition>& newPs , int & modelBP  );
// end of Andy 




double path_log_count(OnePartition & P) ;
double path_log_count(OnePartition_data & P,vector< map< OnePartition , double, CompairSRegs> >& pathctmap);

//void resample(vector<OnePartition_data> &Ps, int halfweight, vector<vector<double> >& c_originals, vector<int> &cutregs) ;    // - by thchiu
void resample(vector<OnePartition_data>* &PsPtr, parameters& p, vector<vector<double> >* &c_originalsPtr, vector<int>* &cutregsPtr) ;    // + by thchiu

double densities_from_partition(vector<double> &testdata,  OnePartition_data &P, int n);
double densities_from_mixpartition(vector<double> &testdata,  vector<OnePartition_data> &Ps, vector<double> wts, int n);
double densities_from_mixpartition(vector<double> &testdata,  vector<OnePartition_data> &Ps, vector<double> wts, int n, double lattic);

//double densities_from_partition(vector<double> &testdata, vector<vector<double> > P); //- by ed520
double densities_from_partition(vector<double> &testdata, vector<vector<double> >& P);  //+ by ed520

// count the number of different partitions:keep a list of different partitions
//(OnePartition_data P_i, double Sum_i,//double Prod_i, vector<int> indexes). when a new one come check its sum and product of the hesh value of all the subregions and compare these two numbers with those of the partitions in the list.
//if overlap occurs, check rigirously(compairRegs) if the partition is already in the list, if yes, add the index of it to the corresponding partition.// use
// use this list for HD_betweenPs.

#endif	/* SAMPLING_H */

