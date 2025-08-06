/*
 * File:   SISfunctions.h
 * Author: Mouse
 *
 * Created on April 5, 2012, 1:28 PM
 */

#ifndef SISFUNCTIONS_H
#define	SISFUNCTIONS_H

//#include "hellingerdist.h"
#include "KLdist.h"
#include "data_store.h"

void Create(parameters & p, vector<OnePartition_data>& OrigPs);


//void SIS( parameters & p, bool discrete, vector<OnePartition_data>& P_bestlevel, int & bestmodelBP, double beta=1.0);  //commented by thchiu
void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, double beta = 1.0);   //added by thchiu

void SIS( parameters & p, bool discrete, OnePartition_data & P_bestlevel, string ofilename, vector<double> &mmax,vector<double> &mmin,  double beta=1.0);   //added by thchiu

void SISforlinearity( parameters & p, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp);


double Output_distance_2 ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect);


void Output_distance ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x), double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect);

void Output_distance_mixnormal ( parameters & p,vector<OnePartition_data> P_bestlevel, double (*density)(const vector<double> &x, mixnormalparas &mixp), mixnormalparas& mixp,  double smoothneighbordist,vector< map< OnePartition, double, CompairSRegs> >& pathctmap, bool pathctcorrect,vector<OnePartition_data> SepPs );

#endif	/* SISFUNCTIONS_H */

