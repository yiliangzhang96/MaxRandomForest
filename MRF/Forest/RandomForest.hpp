//
//  RandomForest.hpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright  2017 1. All rights reserved.
//

#ifndef RandomForest_hpp
#define RandomForest_hpp

using namespace std;

#include <stdio.h>
#include <vector> 

#define dim 10//28 for higgs   //dimension
#define sub_dim 5 //14 for higgs  //dimension subsample
#define Size 5000     //sample size
#define E_size 5000  //Monte Carlo sample size in H-distance calculation
//#define gg 0         // proportion of density transportation
#define threshold_n0 3 /* modify the density if the depth is > certain threshold, in order to reduce overfittinig */
#define importance_depth 10       // first couple of splits only on important features
#define theta 0.001  //coefficient in discrepancy
#define Epsilon 0.001   //non-sense
#define Cov_threshold 0.05//0.05    //coefficient in dimension filtering
#define N_tree 10    //number of trees
#define NN 5     //number of trees in marginal pdf estimation

#define k1 0.9   // resample ratio
#define k2 1     // resample constant

#define K1 0.05  // resample ratio in marginal estimation
#define K2 1     // resample constant

//#define Csize 80000
//const int Bsize = 100000;
//const int Csize = 80000;

//int row,column;
//int r, f;




typedef struct TNode
{
    double density;
    double Fcord;
    double Division_position;
    int Division_cooordinate;
    struct TNode* left;
    struct TNode* right;
    double volume;
    double bound[2*dim];
}TNode;

int RandomForest(TNode *RT, vector< vector<double> > dat, double n, vector<double> margins, double volume, int Depth, double gamma, vector<int> important_features);
double Star_Discrepancy(vector< vector<double> > dat, vector<double> a, vector<double> b, double Row, int P, vector<int> S);
double One_Star_Discrepancy(vector< vector<double> > dat, int P, double a, double b, double Row);
double Imbalance(vector< vector<double> > dat,double a, double b, int d, int k, double volume, double Row, int ndivide);
vector<int> Imb(vector< vector<double> > dat, vector<double> selected_min, vector<double> selected_max, double n_sample, int ndivide, vector<int> feature_candidates);
double One_Dimension_Discrepancy(vector< vector<double> > dat, vector<double> a, vector<double> b, double Row, vector<int> S);
double randomBeta( double alpha, double beta);
double Density_Estimate(vector<double> dat, TNode* Forest[N_tree]);
double Marginal_Estiamte(double dat, TNode* Forest[N_tree], int P);
double Travel_Node(double dat, TNode* D, double sum, int P);
double Eliminate_Forest(TNode* Forest[N_tree]);
void Eliminate_Node(TNode* D);
double Hellinger_Error(TNode* Forest[N_tree], vector< vector<TNode*> > Forestm, vector<double> mmin, vector<double> mmax, double V);
double Marginal_error(double x[10000], double y[10000]);
double Norm_pdf(double dat[dim]);

double Imbalance2(vector<double> dat,double a, double b, int k, double volume, double Row, int ndivide);
vector<double> CopulaForest(TNode *RT, vector<double> dat, double n, vector<double> D, double length, int Deep, double gamma);
vector<double> Copula(TNode* Forest[NN], vector<double> dat, double mmax, double mmin);
//double
double get_cord(double dat, TNode* Forest[NN]);
double get_density(double dat, TNode* Forest[NN]);

vector<int> covariance(vector< vector<double> > dat);


double phi(double x);

#endif /* RandomForest_hpp */
