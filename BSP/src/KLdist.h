/*
 * File:   KLdist.h
 * Author: Mouse
 *
 * Created on May 21, 2012, 8:21 PM
 */

#ifndef KLDIST_H
#define	KLDIST_H

#include "sampling.h"
#include "output.h"
//#include "hellingerdist.h"


vector<double> getvariousKLdist( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x),vector<OnePartition_data>& modePs, vector<double>& modewts);

vector<double> getvariousKLdist_mixnormal( parameters & p,double smoothneighbordist, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp, vector<OnePartition_data>& modePs, vector<double>& modewts,vector<OnePartition_data>& SepPs);

vector<double> all_KL(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic);

double KL_fromf_mixnormal_mix(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double>& wts, int samplesize, double lattic, vector<OnePartition_data>& SepPs );

double KL_fromf_mixnormal(double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas &mixp, vector<vector<double> >& testdata, OnePartition_data &Ps, int samplesize, double lattic, vector<OnePartition_data>& SepPs) ;

double KL_fromf(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, OnePartition_data &P, int samplesize, double lattic);

double KL_mix(double (*density)(const vector<double> &x),vector<vector<double> >& testdata, vector<OnePartition_data> &Ps, vector<double> & wts, int samplesize, double lattic);

double KL_fromf_mixnormal_compare(double (*density)(const vector<double> &x, mixnormalparas &mixp),mixnormalparas &mixp, vector<vector<double> >& testdata, vector<double>  matlabresult);

double KL_fromf_compare(double (*density)(const vector<double> &x), vector<vector<double> >& testdata, vector<double>  matlabresult);

#endif	/* KLDIST_H */

