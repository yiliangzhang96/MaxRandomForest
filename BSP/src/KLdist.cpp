#include "KLdist.h"
#include <cmath>
#include <iostream>


vector<double> getvariousKLdist( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x),vector<OnePartition_data>& modePs, vector<double>& modewts) {
    
//    vector<double> HDs;
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, 0.0));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, 0.0));
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, smoothneighbordist));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, smoothneighbordist));
//    HDs.push_back(KL_fromf(density, p.testdata, modePs[0], p.samplesize, 2 * smoothneighbordist));
//    HDs.push_back(KL_mix(density, p.testdata, modePs, modewts, p.samplesize, 2 * smoothneighbordist));
    return all_KL(density, p.testdata, modePs, modewts, p.samplesize,  smoothneighbordist);
}

vector<double> getvariousKLdist_mixnormal( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp, vector<OnePartition_data>& modePs, vector<double>& modewts,vector<OnePartition_data>& SepPs) {
    vector<double> HDs(6,0);
    HDs[0]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0], p.samplesize,0.0,SepPs);
    HDs[1]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize,0.0,SepPs);
    HDs[2]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0], p.samplesize, smoothneighbordist,SepPs);
    HDs[3]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize, smoothneighbordist,SepPs);
    HDs[4]=KL_fromf_mixnormal(density,mixp, p.testdata, modePs[0], p.samplesize, 2 * smoothneighbordist,SepPs);
    HDs[5]=KL_fromf_mixnormal_mix(density,mixp, p.testdata, modePs, modewts, p.samplesize, 2 * smoothneighbordist,SepPs);
    return HDs;
}


double KL_fromf_mixnormal_mix(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double>& wts, int samplesize, double lattic, vector<OnePartition_data>& SepPs ) {
    int dim = (int)testdata[0].size();
    double result = 0;
    int k=0;
    vector<vector<double> > Ftestdata;
    if ((int) SepPs.size() > 0)  Ftransform(testdata, Ftestdata, SepPs);
    double d1, d2;
    for (int i = 0; i < (int) testdata.size(); i++) {
        d1 = density(testdata[i], mixp);
        if ((int) SepPs.size() > 0) {
            d2 = densities_from_mixpartition(Ftestdata[i], Ps, wts, samplesize, lattic);
            for (int d = 0; d < dim; d++) {
                vector<double> onepiece(1, testdata[i][d]);
                d2 *= one_density(SepPs[d], onepiece, -1);
            }
        } else d2 = densities_from_mixpartition(testdata[i], Ps, wts, samplesize, lattic);
        if(d2>0 && d1>0 && !(log(d2/d1) != log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        }
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}


double KL_fromf_mixnormal(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, OnePartition_data &Ps, int samplesize, double lattic, vector<OnePartition_data>& SepPs) {
    double result = 0;
    double d1, d2;
    int dim = (int)testdata[0].size();
    vector<vector<double> > Ftestdata;
     if ((int) SepPs.size() > 0)  Ftransform(testdata, Ftestdata, SepPs);
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        d1 = density(testdata[i],mixp);
        if((int)SepPs.size()>0){
            d2 = one_density(Ps, Ftestdata[i], lattic);
            for(int d = 0; d<dim; d++){
                vector<double> onepiece(1,testdata[i][d]);
                d2 *= one_density(SepPs[d], onepiece, -1);
            }
        }
        else d2= one_density(Ps, testdata[i], lattic);
        if(d2>0 && d1>0 && !(log(d2/d1) != log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        }
    }
//    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}



double KL_fromf(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, OnePartition_data &P, int samplesize, double lattic) {


    double result = 0;
    double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = density(testdata[i]);
        double d2= one_density(P, testdata[i], lattic);
        if(d2>0 && d1>0 && !(log(d2/d1) != log(d2/d1)))  {
            result+= log(d2/d1);
            k++;
        }
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;
}


double KL_mix(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic) {

    double result = 0;
     double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = density(testdata[i]);
        double d2 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic);
        if(d2>0 && d1>0&& !(log(d2/d1) != log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }

    cout<<"k="<<k<<endl;
    result /=(double)k;
    result = -1.0*result;
    return result;
}

// the best partition should be Ps[0]
// potentially can be faster if combine lattice=0 and >0
vector<double> all_KL(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic) {

    vector<double> dists(6,0);
    double k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {

        double d1 = 1;//density(testdata[i]);
        double d2 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic);
        double d3 = one_density(Ps[0], testdata[i], lattic);
        double d4 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, 0.0);
        double d5 = one_density(Ps[0], testdata[i], 0.0);

        double d6 = densities_from_mixpartition(testdata[i], Ps,wts, samplesize, lattic*2);
        double d7 = one_density(Ps[0], testdata[i], lattic*2);
        if(d2>0 && d1>0 && d3>0 && d4>0 && d5>0 && d6>0&& d7>0&&!(log(d2/d1) != log(d2/d1))) {
             k++;
             dists[0]+= log(d5/d1);
             dists[1]+= log(d4/d1);
             dists[2]+= log(d3/d1);
             dists[3]+= log(d2/d1);
             dists[4]+= log(d7/d1);
             dists[5]+= log(d6/d1);
        }
    }

    cout<<"k="<<k<<endl;
    for(int i=0; i<6; i++) {
        dists[i] /=(double)k;
        dists[i]= -dists[i];
    }
    return dists;
}


double KL_fromf_mixnormal_compare(double (*density)(const vector<double> &x, mixnormalparas &mixp),mixnormalparas &mixp, vector<vector<double> >& testdata, vector<double>  matlabresult) {

    double result = 0;
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        double d1 = density(testdata[i], mixp);
        double d2 = matlabresult[i];
        if(d2>0 && d1>0&& !(log(d2/d1) != log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }
    cout<<"k="<<k<<endl;
    result /=(double)k;
    result =-1.0*result;
    return result;
}

double KL_fromf_compare(double (*density)(const vector<double> &x), vector<vector<double> >& testdata, vector<double>  matlabresult) {

    double result = 0;
    int k=0;
    for (int i = 0; i < (int)testdata.size(); i++) {
        double d1 = density(testdata[i]);
        double d2 = matlabresult[i];
        if(d2>0 && d1>0&& !(log(d2/d1) != log(d2/d1))) {
             k++;
             result+= log(d2/d1);
        }
    }
    result /=(double)k;
    result =-1.0*result;
    return result;
}

