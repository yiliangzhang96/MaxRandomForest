
#include "SISfunctions.h"
#include "timer_ed520.h"    //+ by ed520
extern TIMER timer;

void Create(parameters & p, vector<OnePartition_data>& OrigPs){
    usint_mask um;
    um.mask = 0;
    um.x = 0;
    vector<usint_mask> root(p.dim, um);    //(0,0) should be the root of the sample space.
    shrink_reg_data_short sr;
    sr.num = p.data.size();
    //data_store<vector<double> > originaldata(p.data); //- by ed520
    //sr.regdata = originaldata;    //- by ed520
    sr.reg_code = root;
    OnePartition_data Pt;
    Pt.w = 0;
    Pt.sregs.push_back(sr);
    for (int i = 0; i < p.n; i++) {
        Pt.sregs[0].region_start = i * p.samplesize;    //+ by thchiu + by ed520
        Pt.PartitionID=i;   //+ by thchiu
        OrigPs.push_back(Pt);
    }
 }

bool is_milestone(int i, int plevels){

    if(i==(int)((double)plevels/10.0)) return true;
    if(i==(int)((double)plevels/5.0)) return true;
    if(i==(int)((double)plevels/3.0)) return true;
    if(i==(int)((double)plevels/2.0)) return true;
    if(i==(int)((double)plevels*2.0/3.0)) return true;
    if(i==(int)((double)plevels*4.0/5.0)) return true;
    return false;

}

void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, double beta)   //added by thchiu
{
    vector<double> mmax = vector<double>();
    vector<double> mmin = vector<double>();

    SIS(p, discrete, P_bestlevel, "", mmax, mmin, beta);
}

void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, string ofilename, vector<double>& mmax,vector<double> &mmin, double beta) { //added by thchiu
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerStart_cpu(timer.sTime_C_SIS);
    else
        timer.TimerStart_cpu(timer.sTime_NC_SIS);
#endif
//*** Initialize partition parameters (start) + by FJKu
    cerr<<"beta="<<beta<<endl;
    double* time_array = NULL;
    vector<vector<double> >* c_originals = new vector<vector<double> >;   //+ by ed520
    (*c_originals).resize(p.n);  //+ by ed520
    vector<vector<double> >* c_originalsPtr;    //+ by thchiu
    c_originalsPtr = c_originals;  //+ by thchiu
    vector<int>* cutregs = new vector<int>;    //the region to be cut, thchiu + by ed520
    (*cutregs).resize(p.n);  //+ by ed520
    vector<int>* cutregsPtr; //+ by thchiu
    cutregsPtr = cutregs;  //+ by thchiu
    vector<OnePartition_data> *newPsPtr = new vector<OnePartition_data>; // + by ed520
    Create(p, *newPsPtr);       //+ by ed520
    vector<OnePartition_data>* PsPtr;   //added by thchiu
    PsPtr = newPsPtr;     //added by thchiu + by ed520

    int modelBP=0;
    double currlBPscore = -1.0; // negative one means no cut...


    // added by johnmu: need the case when you don't cut
    P_bestlevel = (*PsPtr)[0];

    double maxlBPscore=currlBPscore;
    int maxlevel=0;
//*** Initialize partition parameters (end) + by FJKu
//*** SIS (start) + by FJKu
    for (int i = 0; i < p.levels; i++) {  // loops for at most level cuts, default: at most 1000 cuts + by FJKu
        cout << "\n\n#i=" << i << endl;
        if ((i > 0) && (p.resampling) && (i % p.steps == 0)) { //(ess<10) &&  with getmax included.
            resample(PsPtr, p, c_originalsPtr, cutregsPtr); //added by thchiu
        }
//*** sample tree (start) + by FJKu: including counting, weight-calculating, and random cut
        bool exhaust = sample_trees(*PsPtr, *c_originalsPtr, *cutregsPtr, p.maxpercentage, time_array, discrete, p);  //+ by thchiu
//*** sample tree (end) + by FJKu: including counting, weight-calculating, and random cut

        if (exhaust) break; // in case there are no more regions to split(letter)
        
//*** Manipulate BPScore (start) + by FJKu
        currlBPscore = -lgamma_c(p.samplesize + (i + 2) / 2.0)- (i + 2) * lgamma_c(0.5) + lgamma_c(0.5 * (i + 2)) - beta * (i + 2);
        
        
        currlBPscore = currlBPscore + getmaxlBP(*PsPtr, modelBP);    //+ by thchiu
        cout << currlBPscore << "\n";

        cout<<"currlBPscore="<<currlBPscore<<endl;    //+ by Yiliang 1/9/18
        cout<<"maxlBPscore="<<maxlBPscore<<endl;    //+ by Yiliang 1/9/18
        
        if (currlBPscore > maxlBPscore - 30) {    // edited by Yiliang 1/9/18
            maxlBPscore = currlBPscore;
            if(currlBPscore > maxlBPscore) maxlBPscore = currlBPscore;
            
            P_bestlevel = (*PsPtr)[modelBP];//+ by thchiu
            cout<<"modelBP="<<modelBP<<endl;
            
            maxlevel=i;
            cout<<"newmax level="<<i+2<<endl;
        }
        if( is_milestone(i, p.levels) && ofilename!="" ) {
             string ofilename2 =ofilename+"_level_"+toStr<int>(i)+".txt";
             ofstream outfile (ofilename2.c_str());
             print_partition(outfile,P_bestlevel, mmax, mmin);
        }
//*** Manipulate BPScore (end) + by FJKu
//*** Check Stopping Criteria (start) + by FJKu
        if(i>maxlevel+10) {
            cout<<"no improvement in 10 levels!\n";
            
            break;
        }
//*** Check Stopping Criteria (end) + by FJKu
		
    }
//*** SIS (end) + by FJKu
    
    delete PsPtr;       //+ by ed520
    delete c_originalsPtr;//+ by ed520
    delete cutregsPtr; //+ by ed520
#ifdef PROFILE
    if(p.dim == 1)
        timer.TimerFinish_cpu( timer.tTime_C_SIS , timer.sTime_C_SIS );
    else
        timer.TimerFinish_cpu(timer.tTime_NC_SIS ,timer.sTime_NC_SIS);
#endif
    return;
}



double Output_distance_2 ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect){

    vector<double>   modewts;
    vector<OnePartition_data> modePs;
    vector<int> modeinds;
    double ess;
    if (pathctcorrect) {
         cout << "P_bestlevel after pathct correction:" << endl;
         for (int par = 0; par < (int) P_bestlevel.size(); par++) {
             P_bestlevel[par].w -= path_log_count(P_bestlevel[par], pathctmap);
         }
    }

    getmax_corr(P_bestlevel, modeinds, modewts, ess);
    modePs.clear();
    
    for (int lala = 0; lala < (int) modeinds.size(); lala++) {
        modePs.push_back(P_bestlevel[modeinds[lala]]);
    }
 //   vector<double> oneHDs = getvarioushellingerdist(p, smoothneighbordist, density, modePs, modewts);

    vector<double> oneKLs = getvariousKLdist(p, smoothneighbordist, density, modePs, modewts);

 //   cout << " unifHD=" << oneHDs[0] << " HD=" << oneHDs[1] << "\n";
 //   cout << " unifHDsmooth=" << oneHDs[2] << " HDsmooth=" << oneHDs[3] << "\n";

    cout << " unifKL=" << oneKLs[0] << " KL=" << oneKLs[1] << "\n";
    cout << " unifKLsmooth=" << oneKLs[2] << " KLsmooth=" << oneKLs[3] << "\n";
    cout << " unifKLsmooth*2=" << oneKLs[4] << " KLsmooth*2=" << oneKLs[5] << "\n";
    return oneKLs[1];

}

void Output_distance ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect){
        Output_distance_2 (p,P_bestlevel, density, smoothneighbordist, pathctmap,  pathctcorrect);
        return;
}


void Output_distance_mixnormal ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp,  double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect,vector<OnePartition_data> SepPs ){

    vector<double>  modewts;
    vector<int> modeinds;

    vector<OnePartition_data> modePs;
    double ess;

    if (pathctcorrect) {
        cout << "P_bestlevel after pathct correction:" << endl;

        for (int par = 0; par < (int) P_bestlevel.size(); par++) {
            P_bestlevel[par].w -= path_log_count(P_bestlevel[par], pathctmap);
        }
    }
    getmax_corr(P_bestlevel, modeinds, modewts, ess);
    modePs.clear();
    for (int lala = 0; lala < (int) modeinds.size(); lala++) {
        modePs.push_back(P_bestlevel[modeinds[lala]]);
    }

//    vector<double> oneHDs = getvarioushellingerdist_mixnormal(p, smoothneighbordist, density, mixp, modePs, modewts, SepPs);
    vector<double> oneKLs = getvariousKLdist_mixnormal(p, smoothneighbordist, density, mixp, modePs, modewts, SepPs);
 //   cout << " unifHD=" << oneHDs[0] << " HD=" << oneHDs[1] << "\n";
 //   cout << " unifHDsmooth=" << oneHDs[2] << " HDsmooth=" << oneHDs[3] << "\n";

    cout << " unifKL=" << oneKLs[0] << " KL=" << oneKLs[1] << "\n";
    cout << " unifKLsmooth=" << oneKLs[2] << " KLsmooth=" << oneKLs[3] << "\n";

    cout << " smoothneighbordist = " << smoothneighbordist << endl;
//    for (int k = 0; k < (int) optmodeinds.size(); k++) {
//        cout << "wt=" << optmodewts[k] << '\n';
//        print_partition(optmodeinds[k]);
//    }

}

// - by thchiu
/*
// for the plot: BPScore is linear of KL divengence
void SISforlinearity( parameters & p, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp) {

    double* time_array = NULL;
    vector<vector<double> > c_originals;
    c_originals.resize(p.n);
//    vector<double> scores;
    vector<int> cutregs;
    cutregs.resize(p.n);
    vector<double> onelevelscores;
    vector<OnePartition_data> newPs;
    Create(p, newPs);
    vector<OnePartition_data>  EmptySepPs;
    double currlBPscore = - 1;

    for (int i = 0; i < p.levels; i++) {
        cout << "\n\n###i=" << i << endl;
        if ((i > 0) && (p.resampling) && (i % p.steps == 0)) { //(ess<10) &&  with getmax included.
           resample(newPs, p.resampling, c_originals, cutregs);
        }
        //bool exhaust = sample_trees(newPs, c_originals, cutregs, p.maxpercentage, time_array,false);  //commented by thchiu
        sample_trees(newPs, c_originals, cutregs, p.maxpercentage, time_array,false);

        for(int j=0; j<(int)newPs.size(); j++){
        currlBPscore = -lgamma(p.samplesize + (i + 2) / 2.0)- (i + 2) * lgamma(0.5) + lgamma(0.5 * (i + 2)) - 1 * (i + 2);
        currlBPscore = currlBPscore + lBPscore(newPs[j]);

        cout << currlBPscore << "\t";
        cout<< KL_fromf_mixnormal(density, mixp, p.testdata, newPs[j], p.samplesize,0,EmptySepPs)<<endl;
        }

    }
    return;
}
*/
