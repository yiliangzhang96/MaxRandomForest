/*
start from \Omega, N trees
1. update trees to add one more cut (based on inexact OPT)
2. calculate w_t = w_t-1*gamma(0.5+n_k^1)gamma(0.5+n_k^2)/gamma(0.5+n_k)(\mu(A_k^1)\mu(A_k^2)/\mu(A_k))^\beta|A_k|^n_k/||||/\alpha_k\lambda_j
3. if necessary, do resampling
*/




#include "sampling.h"
#include "output.h"
#include "count.h"  //+ by thchiu
#include "timer_ed520.h" //+ by ed520
extern TIMER timer;  //+ by ed520

// consider the discrete case: for any region, it will never cut the dimension in which the value of the sample points all agree. If there is no region no dimension to be cut, 'c' will be empty.
// P has t regions: P.sregs.size()=t, then c_original is of length dim*(t-1), cutreg=0,1,...,t-2
void get_info_from_onetree(OnePartition_data &P, vector<double> & c_original,vector<double> &c,int cutreg, double maxpercentage, double &smooth, double *time_array, bool discrete, parameters & para){ //+ by thchiu
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerStart_cpu(timer.sTime_get_info_C);
    else
        timer.TimerStart_cpu(timer.sTime_get_info_NC);
#endif
    int dim = P.sregs[0].reg_code.size();

    int *Tmp_count;
	  Tmp_count = new int[dim];
    bool *Tmp_unique;
	  Tmp_unique = new bool[dim];

//*** Initial partition counting and weight-calculating (start) + by FJKu
    if ((int) c_original.size() == 0) {
        for (int d=0; d<dim; d++)   c_original.push_back(0);
  
        Count(P.sregs[0], Tmp_count, Tmp_unique, dim, discrete, para);    //+ by thchiu

        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[0].num, 0);  //+ by thchiu
    }
//*** Initial partition counting and weight-calculating (end) + by FJKu
//*** Other partition counting and weight-calculating (start) + by FJKu
    else {

        if ((int) c_original.size() != ((int) P.sregs.size() - 1) * dim) {
            cout<<" c_original.size() ="<< c_original.size() <<endl;
            cout<<c_original[0]<<endl;
            cout<<"P.sregs.size()="<<P.sregs.size()<<endl;
            cout << "get_info wrong c_original!\n";
            exit(5);
        }

        for (int d=0; d<dim; d++)   c_original.push_back(0);
        //shrink_reg_data_short & sreg = P.sregs[cutreg];   //- by thchiu

        //Below 4 lines are added by thchiu
        Count(P.sregs[cutreg], Tmp_count, Tmp_unique, dim, discrete, para);
        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[cutreg].num, cutreg);
        Count(P.sregs[(int) P.sregs.size() - 1], Tmp_count, Tmp_unique, dim, discrete, para);
        Calculate_Weight(c_original, Tmp_unique, dim, Tmp_count, P.sregs[(int) P.sregs.size() - 1].num, (int) P.sregs.size() - 1);

        /*{{{*/
        /*vector<usint_mask> new_reg_code;      //- by thchiu
        vector<int> updateindex;
        updateindex.push_back(cutreg);
        updateindex.push_back((int) P.sregs.size() - 1);

        for (int kk = 0; kk < 2; kk++) {
            int i = updateindex[kk];


            //shrink_reg_data_short sreg(P.sregs[i]); // if the ith region is updated   //- by thchiu
            shrink_reg_data_short & sreg = P.sregs[i]; // + by thchiu

            for (int d = 0; d < dim; d++) { // if the ith region is splitted according to the dth dimension
                new_reg_code = sreg.reg_code;
                beta=0.5;
                if ((double) sreg.num / 200.0 < 0.5) beta = max(0.1, (double) sreg.num / 200.0); // heuristic pseudo ct.
                c_original[dim * i + d] = 0;
                for (int j = 0; j < 2; j++) {
                    new_reg_code[d].x = (sreg.reg_code[d].x << 1) + j; // generate a child's regcode
                    new_reg_code[d].mask = (sreg.reg_code[d].mask << 1) + 1;
                    bool unique = false;
                    shrink_reg_data_short  sr;

                    PtsInReg_short(sr, new_reg_code, sreg.regdata, sreg.num, d, unique, false, discrete);

                    if (unique) {// || new_reg_code[d].mask >= 16) {  second part : for letter data(cont)
                        c_original[dim * i + d] = -MAXNUM;
                        break;
                    }
                    c_original[dim * i + d] += lgamma(beta + sr.num);
                }
                if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] -= lgamma(beta + sreg.num);
                if (c_original[dim * i + d] > 1 - MAXNUM) c_original[dim * i + d] += sreg.num * log(2.0);
            }

        }*//*}}}*/

    }
//*** Other partition counting and weight-calculating (start) + by FJKu
//*** manipulation of outputs of weight-calculating (start) + by FJKu
    c = c_original;
    double cmax= c[0];
    int maxind= 0;

     for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] > cmax ){
            cmax = c[i];
            maxind=i;
        }
    }

    if (cmax <1 -MAXNUM) {
        c.clear();
        cout << "no regions no dimensions to be split!!" << endl;
      return;
    }
    double cmax2 = -2 * MAXNUM; // cmax2 is the second maximum of c
    for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] > cmax2 && c[i] < cmax) cmax2 = c[i];
    }
    if (cmax2 <1 -MAXNUM) return;

   for(int i = 0; i < (int) c.size(); i++) {
       if(c[i]!=c[i]) {
                     cout<<"i="<<i<<"c[i]=NAN at a"<<endl;
       }
    }

    if (maxpercentage < 0) smooth = 20;
    else smooth = (cmax - cmax2) / log(maxpercentage / (1 - maxpercentage)); // smooth=1: no smooth
    if(smooth > 1000000) cout<<"smooth==inf! cmax="<<cmax<<" cmax2="<<cmax2<<endl;
    if (smooth > 1) {
        for (int i = 0; i < (int) c.size(); i++) {
            if( c[i] > 1-MAXNUM )           c[i] = c[i] / smooth;
        }
    }
     for(int i = 0; i < (int) c.size(); i++) {
       if(c[i]!=c[i]) {
                     cout<<"i="<<i<<"c[i]=NAN at b"<<endl;
                     cout<<"smooth="<<smooth<<endl;
       }
    }
//*** manipulation of outputs of weight-calculating (end) + by FJKu
#ifdef PROFILE
     if(para.dim == 1)
         timer.TimerFinish_cpu(timer.tTime_get_info_C , timer.sTime_get_info_C);
     else
         timer.TimerFinish_cpu(timer.tTime_get_info_NC , timer.sTime_get_info_NC);
#endif
}

void rand_choose_one_nextleveltree(OnePartition_data &P, vector<double> &c, int &cutreg, double smooth, bool discrete, parameters& para) {// + by thchiu
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerStart_cpu(timer.sTime_rand_choose_C);
    else
        timer.TimerStart_cpu(timer.sTime_rand_choose_NC);
#endif
//*** Generate samples for random cut (start) + by FJKu
    int dim = P.sregs[0].reg_code.size();
    vector<double> ec(dim * P.sregs.size(), 0);
    double cmax = max(c);
    for (int i = 0; i < (int) c.size(); i++) {
        if (c[i] < 1 - MAXNUM) {
            ec[i] = 0;
        } else ec[i] = exp(c[i] - cmax + 5);	// ec comes from calculate weight + by FJKu
        if (c[i]!=c[i]) cout << "is.nan(ec[" << i << "]), c[i]=" << c[i] << endl;
    }
//*** Generate samples for random cut (end) + by FJKu
//*** Do random cut (start) + by FJKu
    int cut = rand_int(ec);
    int cutvar = cut % dim;
    cutreg = (cut - cutvar) / dim;
//*** Do random cut (end) + by FJKu

    shrink_reg_data_short & sreg = P.sregs[cutreg]; //reference, + by thchiu

//*** Update old region and add new region (start) + by FJKu
    //Begin + by thchiu
    double Tmp_upbound;
    sreg.reg_code[cutvar].x = (sreg.reg_code[cutvar].x << 1);
    sreg.reg_code[cutvar].mask = (sreg.reg_code[cutvar].mask << 1) + 1;
    pair<double, double> ranges = convert_ranges(sreg.reg_code[cutvar]);
    Tmp_upbound=ranges.second;

    shrink_reg_data_short sr;


    int data_head = sreg.region_start;
    int data_tail = data_head+sreg.num;
    int Tmp_Count=0;
    for(;data_head != data_tail;){
        if(para.data1D[para.pt_start[data_head]+cutvar] > Tmp_upbound){
            --data_tail;
            int tmp=para.pt_start[data_head];
            para.pt_start[data_head]=para.pt_start[data_tail];
            para.pt_start[data_tail]=tmp;
            //Tmp_Count++;
        }
        else{
            Tmp_Count++;
            ++data_head;
        }
    }

    if(data_head < sreg.region_start + sreg.num)
        if(para.data1D[para.pt_start[data_head]+cutvar] <= Tmp_upbound) // + by ed520
            Tmp_Count++;
    
    sr.region_start = sreg.region_start + Tmp_Count;

    
    sr.num = sreg.num - Tmp_Count;
    sr.reg_code = sreg.reg_code;
    sr.reg_code[cutvar].x ++;

    
    sreg.num = Tmp_Count;
    P.sregs.push_back(sr);
    //End + by thchiu
//*** Update old region and add new region (end) + by FJKu

    if (smooth > 1) P.w = P.w + (smooth - 1) * c[cut] + logsum(c); // + logsum(c); + by thchiu
    else P.w = P.w + logsum(c); //+ by thchiu
    
#ifdef PROFILE
    if(para.dim == 1)
        timer.TimerFinish_cpu(timer.tTime_rand_choose_C , timer.sTime_rand_choose_C);
    else
        timer.TimerFinish_cpu(timer.tTime_rand_choose_NC , timer.sTime_rand_choose_NC);
#endif
}


bool sample_trees(vector<OnePartition_data>& Ps, vector<vector<double> >& c_originals, vector<int>&cutregs,  double maxpercentage, double* time_array, bool discrete, parameters& para) {   //+ by thchiu
    vector<double> c;
    double smooth = 0.0;
    //OnePartition_data newP;   //- by thchiu

    for (int i = 0; i < (int) Ps.size(); i++) {  // loops for M times, M is the number of partition, default: 200
//*** counting and weight-calculating (start) + by FJKu
        get_info_from_onetree(Ps[i], c_originals[i], c, cutregs[i], maxpercentage, smooth, time_array, discrete, para);   //+ by thchiu

        if (max(c) < 1 - MAXNUM) {
            return true;
        }
//*** counting and weight-calculating (end) + by FJKu
//*** random cut (start) + by FJKu
        rand_choose_one_nextleveltree(Ps[i], c,cutregs[i], smooth, discrete, para); //+ by thchiu
//*** random cut (end) + by FJKu
    }
//*** calculate the weights of each partition, for resampling compuation (start) + by FJKu
    double largestwt = Ps[0].w;
    for(int i=0; i<(int)Ps.size(); i++)  {  //problem5: start from 1 is OK, thchiu
        if(Ps[i].w>largestwt)  largestwt= Ps[i].w;
    }
    for (int i = 0; i < (int) Ps.size(); i++) { // normalization, avoid too small or large weight
        Ps[i].w = Ps[i].w - largestwt + 5;

    }
//*** calculate the weights of each partition, for resampling compuation (end) + by FJKu
    return false;
}




// get the best partition so far: the one with the largest weight*number of appearance in the samplings.

// return the maximum of the score among the partitions in the level; modelBP stores the index of the maximum
double getmaxlBP(vector<OnePartition_data>& newPs, int & modelBP) {

    double maxlBP = lBPscore(newPs[0]);
    modelBP=0;
    for (int i = 1; i < (int) newPs.size(); i++) {
            if(lBPscore(newPs[i])>maxlBP){
                 maxlBP = lBPscore(newPs[i]);
                 modelBP=i;
            }
    }
    return (maxlBP);
}

// get the best partition so far: the one with the largest weight*number of appearance in the samplings.
// used for testing only, no need to run in the real data until the last run
void getmax_corr(vector<OnePartition_data>& newPs,  vector<int> & modeinds, vector<double>& modewts, double & ess) {

    vector<double> wts;
    //int dim = newPs[0].sregs[0].reg_code.size();  //commented by thchiu
    for (int i = 0; i < (int) newPs.size(); i++) { // normalization, avoid too small or large weight
        //       newPs[i].w = newPs[i].w - largestwt + 5;
        if((exp(newPs[i].w))!=(exp(newPs[i].w)))     wts.push_back(0.000001);
        else  wts.push_back(exp(newPs[i].w));

    }
    CompairSRegs_data mcompare0;
    sort(newPs.begin(), newPs.end(), mcompare0 );
    CompairSRegs_data mcompair;
    vector<double> wtprodtimes;
    int diffpartitionuntilnow=0;             // also the index of wtprodtimes;
    wtprodtimes.push_back( wts[0] ); //mouse
    modeinds.clear();
    modewts.clear();
    modeinds.push_back(0);

    for (int i = 1; i < (int) newPs.size(); i++) {
        if (!mcompair(newPs[i], newPs[i-1]) && !mcompair(newPs[i-1], newPs[i])) {
            wtprodtimes[diffpartitionuntilnow] += wts[i];
        } else {
            wtprodtimes.push_back(wts[i]);
            modeinds.push_back(i);
            diffpartitionuntilnow++;
        }

    }
    vector<pair<double, double> > prodweightsinfo;
    pair<double, double> onepair;
    for(int i=0; i<=diffpartitionuntilnow; i++){
        onepair.first = modeinds[i];
        onepair.second = wtprodtimes[i];
        prodweightsinfo.push_back(onepair);

    }
    sort(prodweightsinfo.begin(), prodweightsinfo.end(), myfunction4);
    modeinds.clear();
    modewts.clear();
    double largestweight = prodweightsinfo[0].second;
    bool more = true;
    for (int i = 0; more && i <= diffpartitionuntilnow; i++) {
        if ((prodweightsinfo[i].second / largestweight) > 0.0001) {
            modeinds.push_back((int)prodweightsinfo[i].first);
            modewts.push_back(prodweightsinfo[i].second);
        } else more = false;
    }

    cout<<"includeweight.size()="<<(int)modeinds.size()<<'\n';

    return;

}



//return the log of number of paths
/*
   double path_log_count(OnePartition & P) {
   cerr<<"count!\n";
   if (P.sregs.size() <= 3) return 0;
   int dim = P.sregs[0].reg_code.size();


   vector<int> avail_dim(dim, 1);
   for (int i = 0; i < dim; i++) {
   for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
   if (P.sregs[bri].reg_code[i].mask == 0) {
   avail_dim[i] = 0;
   break;
   }
   }
   }

   vector<double> log_count; // log(N) for each dimension, if this dimension is not avail, =0

   for (int i = 0; i < dim; i++) {
   if (avail_dim[i] == 1) {
   OnePartition P1, P2;
   int total1 = -1;
   int total2 = -1;
   for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
   shrink_reg_data_short temp(P.sregs[bri]);
   if ((~(temp.reg_code[i].mask >> 1)) & temp.reg_code[i].x) { // start from 0, represents the first half
   temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
   temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
   P1.sregs.push_back(temp);
   total1 += 1;
   } else { // start from 1, the other half
   temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
   temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
   P2.sregs.push_back(temp);
   total2 += 1;
   }
   }
   log_count.push_back(path_log_count(P1) + path_log_count(P2) + stirling(total1 + total2) - stirling(total1) - stirling(total2));

   }
   }
   return (logsum(log_count));
   }
 */



double path_log_count(OnePartition_data & P, vector< map< OnePartition, double, CompairSRegs> >& pathctmap) {
    int Pregsize = (int) P.sregs.size();
    if (Pregsize <= 3) return 0;
    if (Pregsize <= (int) pathctmap.size()) {
        OnePartition Pwodata;
        for (int k = 0; k < Pregsize; k++) {
            shrink_reg sr;
            sr.reg_code = P.sregs[k].reg_code;
            sr.num = 0;
            Pwodata.sregs.push_back(sr);
        }
        Pwodata.w = 0;

        sort(Pwodata.sregs.begin(), Pwodata.sregs.end(), myfunction3);
        map<OnePartition, double, CompairSRegs>::iterator it;
        it = pathctmap[Pregsize - 1].find(Pwodata);

        //        for( it2 = pathctmap[Pregsize - 1].begin(); it2!= pathctmap[Pregsize - 1].end(); it2++){
        //            CompairSRegs mouse_compare;
        //            cout<<mouse_compare(Pwodata, it2->first)<<" "<< mouse_compare(it2->first, Pwodata)<<"\n";

        //        }

        if (it == pathctmap[Pregsize - 1].end()) {
            cout << "not found in the map!\n";
            cout << "Pregsize-1=" << Pregsize - 1 << "\n";
            cout << "pathctmap[Pregsize-1].size=" << pathctmap[Pregsize - 1].begin()->first.sregs.size() << "\n";
            cout << "to find:\n";
            print_partition(Pwodata);
            cout<<"last one:\n";
            print_partition(pathctmap[Pregsize - 1].begin()->first);
        } else return it->second;
    }
    int dim = P.sregs[0].reg_code.size();

    vector<int> avail_dim(dim, 1);
    for (int i = 0; i < dim; i++) {
        for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
            if (P.sregs[bri].reg_code[i].mask == 0) {
                avail_dim[i] = 0;
                break;
            }
        }
    }
    

    vector<double> log_count; // log(N) for each dimension, if this dimension is not avail, =0

    for (int i = 0; i < dim; i++) {
        if (avail_dim[i] == 1) {
            OnePartition_data P1, P2;
            int total1 = -1;
            int total2 = -1;
            for (int bri = 0; bri < (int) P.sregs.size(); bri++) {
                shrink_reg_data_short temp = P.sregs[bri];
                if ((~(temp.reg_code[i].mask >> 1)) & temp.reg_code[i].x) { // start from 0, represents the first half
                    temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
                    temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
                    P1.sregs.push_back(temp);
                    total1 += 1;
                }
                else { // start from 1, the other half
                    temp.reg_code[i].mask = temp.reg_code[i].mask >> 1;
                    temp.reg_code[i].x = temp.reg_code[i].x & temp.reg_code[i].mask;
                    P2.sregs.push_back(temp);
                    total2 += 1;
                }
            }
            log_count.push_back(path_log_count(P1, pathctmap) + path_log_count(P2, pathctmap) + stirling(total1 + total2) - stirling(total1) - stirling(total2));
            //			cout<<"dim="<<dim<<" log_count="<<log_count[log_count.size()-1]<<endl;
        }
    }
    return (logsum(log_count));
}

/*
//up=mid+half
// low=mid-half

void generate_map(shrink_reg_data sr, vector<double> mid, vector<double> half, int dim, map<vector<usint_mask>, shrink_reg_data, CompairSReg>& newmap, vector<vector<double> >& data, int samplesize, int masklen) {

vector<int> subsubsize(2, 0);
vector<vector<int> > subsize(dim, subsubsize);

vector<vector<double> > onedata;
vector<vector<vector<double> > > d1(2, onedata);
vector< vector<vector<vector<double> > > > subdata(dim, d1);

for (int i = 0; i < samplesize; i++) {
for (int d = 0; d < dim; d++) {
if (data[i][d] < mid[d]) {
subsize[d][0] += 1;
subdata[d][0].push_back(data[i]);
} else {
subsize[d][1] += 1;
subdata[d][1].push_back(data[i]);
}
}
}

for (int d = 0; d < dim; d++) {
shrink_reg_data sr0 = sr;
shrink_reg_data sr1 = sr;
sr0.reg_code[d].x = sr.reg_code[d].x << 1;
sr1.reg_code[d].x = (sr.reg_code[d].x << 1) + 1;

sr0.reg_code[d].mask = (sr.reg_code[d].mask << 1) + 1;
sr1.reg_code[d].mask = sr0.reg_code[d].mask;

vector<double> mid0 = mid;
vector<double> mid1 = mid;
vector<double> half1 = half;

mid0[d] -= half[d]/2;
mid1[d] += half[d]/2;
half1[d]/=2;



if(subsize[d][0]>1 && newmap.find(sr0.reg_code)==newmap.end()){
sr0.num=subsize[d][0];
newmap[sr0.reg_code] = sr0;
if(half1[d]>pow(0.5,masklen) )   generate_map ( sr0, mid0, half1, dim, newmap, subdata[d][0], subsize[d][0], masklen);
}
if(subsize[d][1]>1 && newmap.find(sr1.reg_code)==newmap.end()){
sr1.num=subsize[d][1];
newmap[sr1.reg_code] = sr1;
if(half1[d]>pow(0.5,masklen) )  generate_map ( sr1, mid1, half1, dim, newmap, subdata[d][1], subsize[d][1], masklen);
}

}
return;
}


 */

//void resample(vector<OnePartition_data> &Ps, int halfweight, vector<vector<double> >& c_originals, vector<int> &cutregs) {    //commented by thchiu
void resample(vector<OnePartition_data>* &PsPtr, parameters& p, vector<vector<double> >* &c_originalsPtr, vector<int>* &cutregsPtr) {     //added by thchiu
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerStart_cpu(timer.sTime_resample_C);
    else
        timer.TimerStart_cpu(timer.sTime_resample_NC);
#endif
    //////////////  Ps is changed to *PsPtr in this function, thchiu ///////////////
    int halfweight = p.resampling;
    int size = (int) (*PsPtr).size();
    vector<double> ec(size, 0);
    for (int i = 0; i < size; i++) {
        if (halfweight == 2) {
            (*PsPtr)[i].w /= 2.0;
            ec[i] = exp((*PsPtr)[i].w);
        } else {
            ec[i] = exp((*PsPtr)[i].w);
            (*PsPtr)[i].w = 1;
        }
    }
    //	cout<<endl;
    //vector<OnePartition_data> newPs;  //commented by thchiu
    vector<OnePartition_data>* newPsPtr;  //added by thchiu
    newPsPtr = new vector<OnePartition_data>; //added by thchiu

    //vector<vector<double> > newc_originals; //Problem12: optimized, thchiu    //- by thchiu
    //vector<int> newcutregs;   //- by thchiu
    vector<vector<double> >* newc_originalsPtr; //+ by thchiu
    newc_originalsPtr = new vector<vector<double> >;// + by thchiu
    vector<int>* newcutregsPtr; //+ by thchiu
    newcutregsPtr = new vector<int>; //+ by thchiu

    int* Newdata1D_Index;
    Newdata1D_Index = new int[p.samplesize * p.n];

    for (int i = 0; i < size; i++) {
        int index = rand_int(ec);
        //newPs.push_back((*PsPtr)[index]); //problem3: data copy?, thchiu, commented by thchiu
        //newc_originals.push_back(c_originalsP[index]); //problem3: data_copy?, - by thchiu
        //newcutregs.push_back(cutregs[index]); //- by thchiu
        (*newPsPtr).push_back((*PsPtr)[index]); //+ by thchiu
        (*newc_originalsPtr).push_back((*c_originalsPtr)[index]); //+ by thchiu
        (*newcutregsPtr).push_back((*cutregsPtr)[index]); //+ by thchiu
        (*newPsPtr)[i].PartitionID=i; //+ by thchiu

        for(int j=0; j<p.samplesize; j++){ //+ by thchiu
            Newdata1D_Index[p.samplesize * i + j] = p.pt_start[p.samplesize*index+j];
        }

        for(unsigned j=0; j<(*newPsPtr)[i].sregs.size(); j++){
            (*newPsPtr)[i].sregs[j].region_start = ((*PsPtr)[index].sregs[j].region_start%p.samplesize) + i*p.samplesize;
        }
    }
    //Ps=newPs; // - by thchiu
    //c_originals=newc_originals;   //- by thchiu
    //cutregs=newcutregs;   //- by thchiu
    delete PsPtr;   //+ by ed520
    PsPtr = newPsPtr; // + by thchiu
    //delete c_originalsPtr;     //+ by ed520
    delete c_originalsPtr;  //+ by ed520
    c_originalsPtr = newc_originalsPtr; //+ by thchiu
    //delete cutregsPtr;      //+ by ed520
    delete cutregsPtr;  //+ by ed520
    cutregsPtr = newcutregsPtr; //+ by thchiu
    delete [] p.pt_start;   //+ by thchiu
    p.pt_start = Newdata1D_Index; //+ by thchiu

#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerFinish_cpu(timer.tTime_resample_C , timer.sTime_resample_C);
    else
        timer.TimerFinish_cpu(timer.tTime_resample_NC , timer.sTime_resample_NC);
#endif
    return;
}

// n is the number of points in P
double densities_from_partition(vector<double> &testdata,  OnePartition_data &P, int n){

    int regnum = (int)P.sregs.size();
    int dim = (int) P.sregs[0].reg_code.size();


    pair<double, double> onerange;

    for (int i = 0; i < regnum; i++) {
        bool ind = true;
        for (int d = 0; d < dim; d++) {

            onerange = convert_ranges(P.sregs[i].reg_code[d]);
            if (testdata[d] > onerange.second || testdata[d] <= onerange.first) {
                ind = false;
                break;
            }
        }
        if (ind == true) {
            return ((double) P.sregs[i].regdata.size() / n / exp(lprod_usint_mask(P.sregs[i].reg_code)));

        }
    }
    cout << "wrong densities_from_partition!" << endl;
    for (int d = 0; d < dim; d++) {
        cout<<testdata[d]<<" ";
    }
    exit(10);
    return -1;


}

//the density of a point from a mixture partition density.
// Ps.size()==wts.size(); wts need not to be normalized.

double densities_from_mixpartition(vector<double> &testdata, vector<OnePartition_data> &Ps, vector<double> wts, int n) {

    double den = 0;
    double sumwts = sum(wts);
    for (int i = 0; i < (int) Ps.size(); i++) den = den + densities_from_partition(testdata, Ps[i], n) * wts[i] / sumwts;
    return den;
}

double densities_from_mixpartition(vector<double> &testdata, vector<OnePartition_data> &Ps, vector<double> wts, int n, double lattic) {

    double den = 0;
    double sumwts = sum(wts);
    for (int i = 0; i < (int) Ps.size(); i++) den = den + one_density(Ps[i], testdata, lattic) * wts[i] / sumwts;
    return den;
}


//double densities_from_partition(vector<double> &testdata, vector<vector<double> > P) {  //- by ed520
double densities_from_partition(vector<double> &testdata, vector<vector<double> >& P) { //+ by ed520
    double map_den = -1.0;
    int dim = (int)testdata.size();
    bool found = false;
    for (int j = 0; !found && j < (int) P.size(); j++) {
        found = true;
        for (int d = 0; found && d < dim; d++) {

            double low = P[j][2 * d];
            double high = P[j][(2 * d) + 1];

            if (testdata[d] <= low || testdata[d] > high) {
                found = false;
            }
        }

        if (found) {
            map_den = P[j][2 * dim];
        }
    }

    if (!found) {
        cerr << "Sample not found, idx:\n";
    }

    return map_den;
}

