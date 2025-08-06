//
//  RandomForest.cpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright 2017 1. All rights reserved.
//

#include "RandomForest.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <math.h>
#include <cmath>
#include <vector> 
#include "stdlib.h"

using namespace std;

//#define STACKSIZE 419430400
//#define ndivide 3
//#define sub_dim 5
std::ofstream outfile2;
std::ofstream outfile3;

//int Row;
//int deepth = 0;

//** original rule for spliting **//
// a: left bound   
// b: right bound   
// d: dimension  
// k:number of (quantile) volume: volume of the cube
double Imbalance(vector< vector<double> > dat, double a, double b, int d, int k, double volume, double Row, int ndivide){
    int i;
    double p=0,q = 0;
    for(i = 0; i < Row; i++){
        if((dat[i][d]) < a + k*(b-a)/ndivide){p++;}
        else{q++;}
    }

    //return fabs((p*k/ndivide-q*(1-k/ndivide))/((p*k/ndivide + q*(1-k/ndivide))*volume));
    return fabs((1.0*p*ndivide/k-q*(1.0*ndivide/(ndivide - k)))/((1.0*p*ndivide/k + q*(1.0*ndivide/(ndivide - k)))*(pow(-log2(b-a),2)+1)));
}


//** new rule for spliting **//
double Imbalance2(vector<double> dat,double a, double b, int k, double volume, double Row, int ndivide){
    int i;
    double p=0,q = 0;
    for(i = 0; i < Row; i++){
        if((dat[i]) < a + k*(b-a)/ndivide){p++;}
        else{q++;}
    }
    //return fabs((p*k/ndivide-q*(1-k/ndivide))/((p*k/ndivide + q*(1-k/ndivide))*volume));
    return fabs((1.0*p*ndivide/k-q*(1.0*ndivide/(ndivide - k)))/((1.0*p*ndivide/k + q*(1.0*ndivide/(ndivide - k)))*(pow(-log2(b-a),2)+1)));
}


//** ??? **//
vector<int> Imb(vector< vector<double> > dat, vector<double> selected_min, vector<double> selected_max, double n_sample, int ndivide, vector<int> feature_candidates){
    int i, k, s, num = 0;
    double prod, max = 0;
    vector<int> res(2);
    res[0] = -1;
    res[1] = -1;
    vector<double> test(ndivide);
    vector<double> c(ndivide);

    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(k = 0; k < sub_dim; k++){
        for(i = 1; i < ndivide; i++){
            test[i] = 1.00*i/ndivide;
            num = 0;
            for(s = 0; s < n_sample; s++){
                //k = 0;
                if(((dat[s][feature_candidates[k]]-selected_min[feature_candidates[k]])) < (selected_max[feature_candidates[k]]-selected_min[feature_candidates[k]])*test[i]){
                    num++;
                }
            }
            prod = test[i];
            c[i] = fabs((num/n_sample)-prod);
            //std::cout << "ci: " << c[i] << std::endl;
            if(c[i]>max){max = c[i]; res[0] = feature_candidates[k]; res[1] = i;}
        }
    }
    //std::cout << " dim: " << res[0] << " posit: "<< res[1] << std::endl;
    return res;
}



//** ??? **//
double One_Star_Discrepancy(vector< vector<double> > dat, int P, double a, double b, double Row){
    vector<double> test(100);
    vector<double> c(100);
    int i, j, s, num;
    double prod, max;
    num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<100; i++){
        test[i] = (i+1)*0.01;
        num = 0;
        for(s = 0; s<Row; s++){
            //k = 0;
            if(((dat[s][P]-a)) < (b-a)*test[i]){
                num++;
            }
        }
        
        prod = test[i];
        c[i] = fabs((num/Row)-prod);
        //std::cout << "ci: " << c[i] << std::endl;
        if(c[i]>max){max = c[i];}
    }
    
    return max;
}


//** calculate the 1-dimensional star discrepancy ? **//
double One_Dimension_Discrepancy(vector< vector<double> > dat, vector<double> a, vector<double> b, double Row, vector<int> S){
    int i;
    double n, max = 0;
    for(i = 0; i<sub_dim; i++){
        n = One_Star_Discrepancy(dat, S[i], a[S[i]], b[S[i]], Row);
        //std::cout << "N: " << n << std::endl;
        if(n > max){
            max = n;
        }
    }
    
    return max;
}


//** calculate the star discrepancy in a general case **//
// a,b:boundary for the cube;  
// Row: number of rows needed in calculation
double Star_Discrepancy(vector< vector<double> > dat, vector<double> a, vector<double> b, double Row, int P, vector<int> S){
    vector< vector<double> > test(50000, vector<double>(dim));
    vector<double> c(50000);
    int i, j, k, s, num;
    double prod, max;
    //std::cout << & num << std::endl;
    //num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[20000];  j loop to initialize test;  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<50000; i++){
        //srand(time(NULL));
        for(j = 0; j<sub_dim; j++){
            test[i][S[j]] = (rand()%999)+1;
            test[i][S[j]] = test[i][S[j]]/1000;
        }
        
        num = 0;
        for(s = 0; s<Row; s++){
            k = 0;
            for(j = 0; j<sub_dim; j++){
                if(((dat[s][S[j]]-a[S[j]])/(b[S[j]]-a[S[j]])) < test[i][S[j]]){
                    k++;
                }
            }
            if(k==sub_dim){
                num++;
            }
        }
        prod = 1;
        for(j = 0; j< sub_dim; j++){
            prod = prod * test[i][S[j]];
        }
        c[i] = fabs((num/Row)-prod);
        if(c[i]>max){max = c[i];}
    }
    //std::cout << "MAX: " << max <<std::endl;
    return max;
}




//** construct Max-Random Forest **//
//RT:tree  
//dat:data that going to generate branch  
//n:number of data inside   
//D:cube information
int RandomForest(TNode *RT, vector< vector<double> > dat, double n, vector<double> margins, double volume, int Depth, double gamma, vector<int> important_features)
{
    int i,j,k;
    double best_score = 0, imbalance_score;
    int best_feature = -1, n_left, n_right;
    double best_position = -1; 
    int ndivide = 3;
    vector<int> feature_candidates(sub_dim);

    Depth++;

    //vector< vector<double> > sample(Size, vector<double>(dim));
    
    /* if no samples left in the current one */
    // if(n==0){
    //     RT->density = 0;
    //     if(Depth > threshold_n0*dim+1){
    //         RT->density = gamma/(Size*volume);
    //         for(i = 0;i < 2*dim; i++){RT->bound[i] = margins[i];}
    //         RT->volume = volume;
    //     }
    //     //std::cout << "n==0! Depth = " << Depth << "; " << "v = "<< volume << "; " << "Density = " << RT->density << std::endl;
    //     return 0;
    // }
    
    // if(Depth == 0){
    //     outfile3.open ("00.txt");
    // }

    //parameter manipulation
    vector<double> bound_min(dim), bound_max(dim); // A, B
    for(i = 0; i < 2*dim; i++){
        if(i%2 == 0){bound_min[i/2] = margins[i];}
        else{bound_max[(i-1)/2] = margins[i];}
    }
    RT->density = -1;

    /* copy the important features */
    for(i = 0; i < important_features[0]; i++){feature_candidates[i] = important_features[i+1];}
    
    // If Depth > importance_depth then will be back to original dimension
    int Dim = important_features[0];
    if(Depth > importance_depth){
        Dim = sub_dim;
        i = 0;
        while(i < sub_dim/*J[0]*/){
            k = rand()%dim;
            j = 0;
            while(j < i && k != feature_candidates[j]){
                j++;
            }
            if(j == i){
                feature_candidates[i] = k;
                i++;
            }
        }
    }
    //std::cout << "  First of all, v = " << volume << std::endl;
    
    ///*** calculate the imbalance score ***///
    int selected_feature;
    double selected_min, selected_max;
    for(k = 0; k < Dim/*J[0]*/; k++){
        for(i = 1; i < ndivide; i++){
            selected_feature = feature_candidates[k];
            selected_min = margins[2*feature_candidates[k]];
            selected_max = margins[2*feature_candidates[k]+1];
            imbalance_score = Imbalance(dat, selected_min, selected_max, selected_feature, i, volume, (int)(k1*n+k2), ndivide);
            std::cout << imbalance_score << std::endl;
            if(best_score < imbalance_score && (selected_max - selected_min > 0*pow(0.5,4))){ //Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n)
                best_score = imbalance_score;//Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n);
                best_feature = selected_feature; //split dimension
                best_position = i; //split position
            }
        }
    }

    /* DSP?
    vector<double> S = Imb(dat, selected_min, selected_max, n, ndivide, feature_candidates);
     Nm = S[0];
     Im = S[1];
    */



    /* the important features are selected completely at random */
    // best_feature = feature_candidates[rand()%sub_dim];
    // best_position = rand()%(ndivide-1) + 1;
    

    /* corner case: all directions are completely even, stop splitting */
    if(best_feature == -1 && Depth > 1){
        RT->density = gamma*n/(Size*volume);
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = margins[i];
        }
        RT->volume = volume;
        std::cout << "Even! Depth = " << Depth << "; " << "v = " << volume << "; " << "Density = " << RT->density << std::endl;
        return 0;
    }
    
    //////****************************************************************/////
    //////******                      Split                     **********/////
    //////****************************************************************/////
    
    /* in practice, seems that a fixed depth is also fine?  */

    double threshold = theta * pow(k1*Size,0.5)/(k1*n+k2);
    if(((Depth < 10) /*&& ((threshold < Epsilon) || (One_Dimension_Discrepancy(sample, A, B, (int)(k1*n+k2),s) > threshold ) || (Star_Discrepancy(sample, A, B, (int)(k1*n+k2), Nm,s) > threshold  ))*/)){
        std::cout << "split!";
        //double left_data[Size][dim];
        vector< vector<double> > left_data(Size, vector<double>(dim));
        //double right_data[Size][dim];
        vector< vector<double> > right_data(Size, vector<double>(dim));
        
        ///*** divide the data according to the split position ***///
        double Va = (best_position * bound_max[best_feature] + (ndivide - best_position) * bound_min[best_feature])/ndivide;
        n_left = n_right = 0;
        for(i = 0; i < n; i++){
            if(dat[i][best_feature] < Va){
                for(j = 0; j < dim; j++){
                    left_data[n_left][j] = dat[i][j];
                }
                n_left++;
            }
            else{
                for(j = 0; j < dim; j++){
                    right_data[n_right][j] = dat[i][j];
                }
                n_right++;
            }
        }
        
        /* corner case: either one of the child node has no data, stop */
        if((n_right == 0 || n_left == 0)){
            RT->density = n/(Size*volume);
            //RT->bound[0] = D[2*Nm];
            for(i = 0;i < 2*dim; i++){
                RT->bound[i] = margins[i];
            }
            RT->volume = volume;
            std::cout << "Leaf with n ->> 0! Depth = " << Depth << "; " << "Density = " << RT->density << std::endl;
            return 0;
        }
        
        /* coefficient used when no sample is left */
        double lgamma = gamma,rgamma = gamma;
        
        ///*** grow two children ***///
        TNode* left;
        TNode* right;
        
        left = (TNode*)malloc(sizeof(TNode)*1);
        right = (TNode*)malloc(sizeof(TNode)*1);
        
        left->left = NULL;
        left->right = NULL;
        right->left = NULL;
        right->right = NULL;
        
        RT->left = left;
        RT->right = right;
        //n_left to calculate the sample in left child;   left_data to be passed into left
        
        //0: coordinate;    1: position
        RT->Division_cooordinate = best_feature;
        RT->Division_position = Va;
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = margins[i];
        }
        std::cout << " At " << best_feature << std::endl;
        
        ///*** left child ***///
        vector<double> left_margins(2*dim);
        for(i=0;i<2*dim;i++){left_margins[i] = margins[i];}
        left_margins[2*best_feature+1] = Va;
        RandomForest(left, left_data, n_left, left_margins, volume*best_position/ndivide, Depth,lgamma, important_features);
        
        ///*** right child ***///
        vector<double> right_margins(2*dim);
        for(i=0;i<2*dim;i++){right_margins[i] = margins[i];}
        right_margins[2*best_feature] = Va;
        RandomForest(right, right_data, n_right, right_margins, volume*(1-(best_position/ndivide)), Depth,rgamma, important_features);
    }
    
    ///*** not meet the split condition, stop and calculate the density ***///
    else{
        RT->density = gamma*n/(Size*volume);
        //RT->bound[0] = D[2*Nm];
        for(i = 0;i < 2*dim; i++){
            RT->bound[i] = margins[i];
        }
        RT->volume = volume;
        std::cout << "Leaf! Depth = " << Depth << "; " << "Density = " << RT->density << std::endl;
    }
    
    // if(Depth == 1){
    //     outfile3.close ();
    // }
    // if(RT->density == 0){
    //     int aaa;
    //     aaa = 1+1;
    // }

    return 0;
}




//** generate Beta-distributed variables **//
double randomBeta( double alpha, double beta)
{
    /*Johnk's beta generator*/
    double u, v;
    double x, y;
    
    //std::uniform_real_distribution<double> rb(0,5);
    do
    {
        //std::default_random_engine random1;
        u = rand() %3000;
        std::cout << u/1000 << std::endl;
        //std::default_random_engine random2;
        v = rand() %3000;
        std::cout << v/1000 << std::endl;
        x=pow(u/1000,1/alpha);
        y=pow(v/1000,1/beta);
    } while (x+y>1);
    return x/(x+y);
}


//**  get estimated density from the forest **//
double Density_Estimate(vector<double> dat, TNode* Forest[N_tree])
{
    int i;
    double density=0;
    for (i = 0;i < N_tree; i++){
        TNode* T = Forest[i];
        while(T->left != NULL){
            if(dat[T->Division_cooordinate] < T->Division_position){
                T = T->left;
            }
            else{
                T = T->right;
            }
        }
        density = density + T->density;
    }
    return density/(N_tree);
}

//** ??? **//
double Marginal_Estiamte(double dat, TNode* Forest[N_tree], int P){
    int i;
    double M = 0, marg[N_tree];
    for(i = 0; i < N_tree; i++){
        TNode* T = Forest[i];
        //sum = 0
        marg[i] = Travel_Node(dat, T, 0, P);
        //if(i == 9){
        M = M + marg[i];
        //}
        if(dat == 0.001 || dat == 0.998){
            
        }
    }
    
    return M/(N_tree);
}

//** delete the forest **//
double Eliminate_Forest(TNode* Forest[N_tree]){
    int i;
    outfile2.open ("calculate.txt");
    for(i = 0; i < N_tree; i++){
        TNode* T = Forest[i];
        Eliminate_Node(T);
    }
    outfile2.close();
    return 0;
}

//** eliminate the subtree with the given node **//
void Eliminate_Node(TNode* D){
    if(D->left->density != -1){
        free(D->left);
        outfile2 << 1 << std::endl;
    }
    else{
        Eliminate_Node(D->left);
    }
    
    if(D->right->density != -1){
        free(D->right);
    }
    else{
        Eliminate_Node(D->right);
    }
    
    free(D);
}



//** traverse the node and get the sum of the nodes **//
double Travel_Node(double dat, TNode* D, double sum, int P){
    if(D->density == -1){
        if(D->Division_cooordinate == P){
            if(dat < (D->Division_position)){
                sum = Travel_Node(dat, D->left, sum, P);
            }
            
            else{
                sum = Travel_Node(dat, D->right, sum, P);
            }
        }
        else{
            sum = Travel_Node(dat, D->left, sum, P);
            sum = Travel_Node(dat, D->right, sum, P);
        }
    }
    
    else{
        if(D->density == 0){
            sum = sum;
        }
        else{
            sum = sum + (D->density)*(D->volume)/(D->bound[2*P+1] - D->bound[2*P]);
        }
    }
    return sum;
}

//** calculate the star discrepancy **//
double Star_Discrepancy2(double dat[Size][dim],double a[dim], double b[dim], double Row, int P){
    vector<double> test(100);
    vector< vector<double> > c(dim, vector<double>(100));
    int i, j, k, s, num;
    int prod, max;
    num = 0;
    max = 0;
    
    //generate vector test
    //i loop for test[100];  s loop for dat[];  sj loop for element in dat[][]
    for(i = 0; i<100; i++){
        test[i] = (i+1)/100;
        for(k = 0; k<dim; k++){
            for(s = 0; s<Row; s++){
                //k = 0;
                if(((dat[s][P]-a[P])/b[P]) < test[i]){
                    num++;
                }
            }
            
            prod = test[i];
            c[k][i] = fabs((num/Row)-prod);
            if(c[k][i]>max){max = c[k][i];}
        }
    }
    
    return max;
}

//** ??? **//
double Marginal_error(double x[10000], double y[10000]){
    int i,j;
    double f, E = 0, F = 0;
    for(i = 0; i < 10000;i++){
        f = (cos(x[i]*3.14*3.5)+1)/(sin(3.14*3.5)/(3.14*3.5) + 1);
        E = E + pow(f*y[i],0.5);
        F = F + pow(f-y[i],2);
    }
    
    std::cout << "Marginal_HD: " << 1 - E/10000 << std::endl;
    std::cout << "Marginal_MISE: " << F/10000 << std::endl;
    
    return 0;
}

//** calculate the hellinger distance **//
double Hellinger_Error(TNode* Forest[N_tree], vector< vector<TNode*> > Forestm, vector<double> mmin, vector<double> mmax, double V){
    int i,j,k;
    double F_out, E_out;
    vector<double> dat(dim),datm(dim-1);
    double E = 0, F = 0;
    double f,g,fm;
    double dat1,dat2;
    double gde,fde,Em = 0, Hm = 0;
    
    //KL-data
    std::ifstream outfile4;
    //outfile4.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/100dim/10j+90m/kl.txt");
    outfile4.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/test_data.txt");
    
    // //???
    // std::ifstream outfile6;
    // outfile6.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/100dim/2j+98m/kl1.txt");
    
    // //specific dimension
    // std::ofstream outfile2;
    // outfile2.open ("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/100dim/2j+98m/kl2.txt");

    //true pdf for KL data
    std::ifstream outfile7;
    //outfile7.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/100dim/10j+90m/pdf.txt");//pdf2 for 98m
    outfile7.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/test_data_pdf.txt");
    
    for(i = 0; i < E_size/*200000*/; i++){
        f = 1;
        g = 1;
        fm = 1;
        for(j = 0; j < dim; j++){
            //input the KL data
            outfile4 >> dat[j];
            
        }

        //input the true pdf of KL data
        outfile7 >> fm;
        
        TNode* ForestM[NN];
        /*for(j = 0; j < dim; j++){
            ForestM[j] = Forestm[i][j];
        }*/
        //for(j = 0; j < dim-1; j++){
        //      outfile6 >> datm[j];
        //}
        
        double marginal_pdf;
        //anti-copula
        for(j = 0; j<dim;j++){
            for(k = 0; k < NN; k++){
                ForestM[k] = Forestm[j][k];
            }
            //to prevend ttt be 0
            marginal_pdf = (get_density(dat[j], ForestM) + 0.0000)/1.0000;
            if(marginal_pdf == 0){
                marginal_pdf = 0.001;
            }
            
            //if KL data is outside [min,max]
            if(dat[j] > mmax[j] || dat[j] < mmin[j]){
                g = fm;
                break;
            }
            
            //calculate the coordinate after copula
            dat[j] = get_cord(dat[j], ForestM);
            
            //if(j < 2){marginal_pdf = 1;}
            //if(j == 10){
            //    outfile2 << marginal_pdf << std::endl;
            //}
            // multiplies the ith marginal density
            g = g*marginal_pdf;
            //if(j == 0||j == 1){g = g*ttt;}
        }
        
        fde = datm[0]/*(dat1*dat2)*/;
        //take the copula density
        gde = Density_Estimate(dat, Forest);
        //density equals to copula pdf multiplies n marginal pdfs
        g = g*gde;
        std::cout << " f: " << 1*fm << " Est:" << g << std::endl;
        //outfile2 << f << "; " << g << std::endl;
        fm = fm*1;
        
        // is it possible for fm to be 0 ?
        if(fm == 0){
            Em = Em + 1;
            Hm = Hm + 1/fm;
            E = E + 1;  //pow(g*1, 0.5);//0.5*pow((sqrt(g) - sqrt(f)),2)/f;  //pow(g*g,0.5);
            F = F;
        }
        else{
            Hm = Hm + 1/fm;
            Em = Em + pow(g/fm,0.5);
            E = E + pow(g/fm,0.5);  //Hellinger distance //pow(g*1, 0.5);//0.5*pow((sqrt(g) - sqrt(f)),2)/f;
            F = F + log(fm/g); //KL-divergence
        }
    }
    
    //outfile2.close();
    outfile4.close();
    //outfile7.close();
    //outfile6.close();

    std::cout << "ADG: " << F/(1*E_size) << std::endl;
    
     F_out = F/E_size;
     E_out = E/E_size;
    
    return 1 - 1*E/(1*E_size);
}



//** calculate the pdf after normalization? **//
double Norm_pdf(double dat[dim]){
//for(j = 0; j < dim; j++){
//dat[j] = (rand()%10000)*0.0001;}
 /*   int i,j;
    std::ifstream outfile4;
    outfile4.open("/Users/zhang/Desktop/831test.txt");

    double dat[200000][dim];
    for(i = 0;i<200000;i++){
        for(j = 0;j<dim;j++){
            outfile4 >> dat[i][j];
        }
    }
 */
    
return 1;
}


//** copula forest **//
vector<double> CopulaForest(TNode *RT, vector<double> dat, double n, vector<double> D, double length, int Depth, double gamma)
{
     int i,j,k,h;
        double max = 0, imbalance_score;
        int n_left, n_right;
        double Im = -1,test,rou;
        vector<double> TEST(dim);
        vector<double> sample(Size);
        int after;
        int ndivide = 3;
    
    if(Depth == 0){
        rou = 0;
        RT->Fcord = 0.5;
        //rou = 0;
    }
        if(Depth == 0){
            outfile3.open ("00.txt");
        }
        for(k = 0; k < (int)(K1*n + K2); k++){
            h = rand()%(int)(n);
            sample[k] = dat[h];
            outfile3 << dat[h] << std::endl;
        }
    
        RT->density = -1;
        Depth++;
        //parameter manipulation
        double A, B;
        A = D[0];
        B = D[1];
    
        //RT->density = -1;
        //subsample
        //std::default_random_engine random(time(NULL));
        //std::uniform_int_distribution<int> unif3(0,dim-1);

        i = 0;
        std::cout << "    Cop, v = " << length << std::endl;       
        //find the imbalance information
            
            for(i = 1; i < ndivide; i++){
                //std::cout << n << "haha";
                imbalance_score = Imbalance2(sample, D[0], D[1], i, length, (int)(K1*n+K2), ndivide);
                std::cout << imbalance_score << std::endl;
                if(max < imbalance_score && (D[1] - D[0] > 0*pow(0.5,4))){ //Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n)
                    max = imbalance_score;//Imbalance(dat, D[2*s[k]-1], D[2*s[k]], s[k], i, volume, n);
                    //Nm = s[k]; //split dimension
                    Im = i; //split position
                }
            }

        /*
         int* S = Imb(dat, A, B, n, ndivide,s);
         Nm = S[0];
         Im = S[1];
         */
        //Nm = s[rand()%sub_dim];
        //Im = rand()%(ndivide-1) + 1;
        
        if(max == 0 && Depth > 1){
            RT->density = gamma*n/(Size*length);
            for(i = 0;i < 2; i++){
                RT->bound[i] = D[i];
            }
            RT->volume = length;
            //std::cout << "Even! Depth = " << Depth << "; " << "v = " << length << "; " << "Density = " << RT->density << std::endl;
            for(i = 0; i < n; i++){
                dat[i] = RT->Fcord + RT->density * (dat[i] - (RT->bound[1] + RT->bound[0])*0.5);//+ 0*(((double) rand() / (RAND_MAX)) * 2 - 1)*0.5*RT->volume;
                    }
            return dat;
        }
        
        double threshold = theta * pow(K1*Size,0.5)/(K1*n+K2);
        
        //One_Star_Discrepancy(dat, Nm, A[Nm], B[Nm], n);
        //Star_Discrepancy(dat, A, B, n, Nm);
        //std::cout << "  N: " << n <<"  One: " << One_Dimension_Discrepancy(dat, A, B, n) << "  Star: " << Star_Discrepancy(dat, A, B, n, Nm) << "  Threshold: " << threshold << std::endl;
        
        //////*****      Split Condition     *****///////(One_Star_Discrepancy(dat, Nm, A[Nm], B[Nm], n) > threshold )
        //(n > 5 && Depth < 80) ||
        if(((Depth < 8) /*&& ((threshold < Epsilon) || (One_Dimension_Discrepancy(sample, A, B, (int)(k1*n+k2),s) > threshold ) || (Star_Discrepancy(sample, A, B, (int)(k1*n+k2), Nm,s) > threshold  ))*/)){
            //std::cout << "split!";
            
            vector<double> left_data(Size);
            vector<double> right_data(Size);
            
            // Va split absolute position
            double Va = 1.0*(Im * D[1] + (ndivide - Im) * D[0])/ndivide;
            n_left = n_right = 0;
            for(i = 0; i < n; i++){
                if(dat[i] < Va){
                    left_data[n_left] = dat[i];
                    n_left++;
                }
                else{
                    right_data[n_right] = dat[i];
                    n_right++;
                }
            }
            
            double lgamma = gamma,rgamma = gamma;
            // if(Depth > threshold_n0*dim && (n_right == 0)){
            //     rgamma = gg*gamma*n;
            //     lgamma = (1-gg)*gamma;
                
            // }
            
            // if(Depth > threshold_n0*dim && (n_left == 0)){
            //     lgamma = gg*gamma*n;
            //     rgamma = (1-gg)*gamma;
                
            // }
            
            vector<double> dat_left = left_data;
            vector<double> dat_right = right_data;
            
            //std::cout << "  One: " << One_Dimension_Discrepancy(dat, A, B, n) << "  Star: " << Star_Discrepancy(dat, A, B, n, Nm) << "  Threshold: " << threshold << std::endl;
            //build left child; right child
            TNode* left;
            TNode* right;
            left = (TNode*)malloc(sizeof(TNode)*1);
            right = (TNode*)malloc(sizeof(TNode)*1);
            left->left = NULL;
            left->right = NULL;
            right->left = NULL;
            right->right = NULL;
            RT->left = left;
            RT->right = right;
            
            //n_left to calculate the sample in left child;   left_data to be passed into left
            //0: coordinate;    1: position

            RT->Division_position = Va;
            for(i = 0;i < 2; i++){
                RT->bound[i] = D[i];
            }

            //left branch
            //std::cout << "left: " << volume << " Im:" << Im << " prod: " << Im/ndivide << " In: " << volume*Im/ndivide << std::endl;
            vector<double> LD(2);
            for(i=0;i<2;i++){LD[i] = D[i];}
            LD[1] = Va;

            if(n_left > 0){
            left->Fcord = RT->Fcord - 0.5*n/Size + 0.5*n_left/Size;
            //?
            dat_left = CopulaForest(left, left_data, n_left, LD, length*Im/ndivide, Depth,lgamma);
                
            }
            else{
                left->density = 0;
                left->Fcord = RT->Fcord - 0.5*n/Size;
            }
            //right branch
            //std::cout << "right: " << volume << " Im:" << Im << " prod: " << (1-(Im/ndivide)) << " In: " << volume*(1-(Im/ndivide)) << std::endl;
            vector<double> RD(2);
            for(i=0;i<2;i++){RD[i] = D[i];}
            RD[0] = Va;
            if(n_right > 0){
            right->Fcord = RT->Fcord - 0.5*n/Size + 1.0*n_left/Size + 0.5*n_right/Size;
            dat_right = CopulaForest(right, right_data, n_right, RD, length*(1-(Im/ndivide)), Depth,rgamma);
            }
            else{
                right->density = 0;
                right->Fcord = RT->Fcord - 0.5*n/Size + 1.0*n_left/Size;
            }
            
            int l1 = 0, l2 = 0;
            for(i = 0; i< n;i++){
                if(dat[i] < Va){
                    if(dat[i] == dat_left[l1]){
                        std::cout<<1;
                    }
                    dat[i] = dat_left[l1];
                    l1++;
                }
                else{
                    if(dat[i] == dat_right[l2]){
                        std::cout<<1;
                    }
                    dat[i] = dat_right[l2];
                    l2++;
                }
            }
            
        }
        
        //output density
        else{
            RT->density = gamma*n/(Size*length);
            //RT->bound[0] = D[2*Nm];
            for(i = 0;i < 2; i++){
                RT->bound[i] = D[i];
            }
            RT->volume = length;
            std::cout << "Leaf! Depth = " << Depth << "; " << "Density = " << RT->density << std::endl;
            for(i = 0; i < n; i++){
                dat[i] = RT->Fcord + RT->density * (dat[i] - (RT->bound[1] + RT->bound[0])*0.5);//+ 0*(((double) rand() / (RAND_MAX)) * 2 - 1)*0.5*RT->volume;
            }
        }
        
        if(Depth == 1){
            outfile3.close ();
        }

    /*
    for(i = 0; i < n; i++){
        std::cout << dat[i] <<std::endl;
    }
    std::cout << "#####################################" << Depth <<std::endl;
    */
    return dat;
    
}

//** construct the copula forest  **//
vector<double> Copula(TNode* Forest[NN], vector<double> dat, double mmax, double mmin){
    
    int i,j;
    int T = 0;
    vector<double> dat1(Size);
    vector<double> sample(Size);
    vector<double> sample1(Size);
    vector<double> D(2);
    D[0] = mmin;
    D[1] = mmax;
    vector<double> datt;
    
    //TNode* *TList = new TNode*[10];
    std::cout << "hehe" << dat[0] <<std::endl;
    for(i = 0; i < Size; i++){
        sample[i] = dat[i];
    }
   
    std::cout << dat[0] <<std::endl;
    //std::ofstream outfile6;
    //outfile6.open ("source1.txt");

    for(i = 0; i < NN; i++){
            //Forest[T] = new TNode;
        for(j = 0; j < Size; j++){
            sample1[j] = sample[j];
            //if(mmax == 1){
            //    outfile6 << sample1[i]
            //}
        }

        datt = CopulaForest(Forest[T],sample1, Size, D,(mmax - mmin), 0, 1);
        for(j = 0; j < Size; j++){
            dat[j] = (dat[j]*i*1.00 + datt[j])/(1.0*(i+1));
        }
        T++;
    }
    //outfile6.close();
    //*dat = *dat * 0.1;
    //Eliminate_Forest(TList);
    return dat;
}


//** get coordinate? **//
double get_cord(double dat, TNode* Forest[NN]){
    int i;
    double Fcord=0;
    double left = 0;
    double right = 1;
    for (i = 0;i < NN; i++){
        TNode* T = Forest[i];
        while(T->left != (TNode*)NULL && (T->density == -1)){
            if(dat < T->Division_position){
                T = T->left;
                left = T->bound[0];
                right = T->bound[1];
            }
            else{
                T = T->right;
                left = T->bound[0];
                right = T->bound[1];
            }
        }
        
        Fcord = Fcord + T->Fcord + T->density * T->volume * 1 * (dat - (right + left)*0.5)/(1 * (right - left));
    }
    return Fcord/(NN);
}

//** calculate the density from the forest given a univariate value dat **//
double get_density(double dat, TNode* Forest[NN]){
    int i;
    double density=0;
    
    for (i = 0;i < NN; i++){
        TNode* T = Forest[i];
        while(T->left != NULL && (T->density == -1)){
            if(dat < T->Division_position){  
                T = T->left;
            }
            else{
                T = T->right;
            }
        }
        density = density + T->density;
    }
    return density/(NN);

}

//** calculate th covariance matrix? **//
vector<int> covariance(vector< vector<double> > dat){
    vector< vector<double> > Cov_vector(dim, vector<double>(dim));
    vector< vector<double> > Cov(dim, vector<double>(dim));
    vector< vector< vector<double> > > MC(dim, vector< vector<double> >(dim, vector<double>(20)));
    vector< vector< vector<double> > > count(dim, vector< vector<double> >(dim, vector<double>(10)));
    vector< vector<double> > Test(dim, vector<double>(dim));
    vector< vector<double> > Ht(dim, vector<double>(dim));
    vector<double> mean(dim);
    int s;
    int aux;
    double max = 0;
    int maxj;
    int maxk;
    vector<int> J(dim);
    int i, j, k, g;
    
    for(j = 0; j < dim; j++){
        mean[j] = 0;
        for(k = j; k < dim; k++){
            //for(i = 0; i < 10; i++){
                //MC[j][k][i] = (rand()%999 + 1)*0.001;
                //MC[j][k][i+10] = (rand()%999 + 1)*0.001;
                //count[j][k][i] = 0;
                Ht[j][k] = 0;
            //}
        }
    }
    
    for(i = 0; i < Size; i++){
        s = 0;
        for(j = 0; j < dim; j++){
            mean[j] = mean[j] + dat[i][j];
            for(k = j; k < dim; k++){
                //Cov_vector[(dim + dim - (j-1))*(j)/2 + k-j+1] = Cov_vector[(dim + dim - (j-1))*(j)/2 + k-j+1] + dat[i][j]*dat[i][k];
                Cov_vector[j][k] = Cov_vector[j][k] + dat[i][j]*dat[i][k];
                //for(g = 0; g < 10; g++){
                //    if(dat[i][j] < MC[j][k][g] && dat[i][k] < MC[j][k][g+10]){count[j][k][g]++;}
                //}
                //new 12-7
                if(dat[i][j]/dat[i][k] < 0.5 || dat[i][j]/dat[i][k] > 2){
                    Ht[j][k]++;
                }
            }
        }
    }
    s = 1;
    for(j = 0; j < dim-1; j++){

        for(k = j+1; k < dim; k++){
            Cov[j][k] = fabs((mean[j]*mean[k]/(Size) - Cov_vector[j][k])/pow((Cov_vector[j][j] - mean[j]*mean[j]/Size)*(Cov_vector[k][k] - mean[k]*mean[k]/Size),0.5));
            if(Cov[j][k] > max){
                max = Cov[j][k];
                maxj = j;
                maxk = k;
            }
            
            //add j,k to J[]
            if(Cov[j][k] > Cov_threshold){
                aux = 1;
                while(J[aux] != j && aux < s){
                    aux++;
                }
                if(aux == s){
                    J[s] = j;
                    s++;
                }

                aux = 1;
                while(J[aux] != k && aux < s){
                    aux++;
                }
                if(aux == s){
                    J[s] = k;
                    s++;
                }
            }
            
            /*
            for(g = 0; g < 10; g++){
                if(count[j][k][g] == 0){Test[j][k] = Test[j][k] + 0.3;}
                else{
                Test[j][k] = Test[j][k] + fabs((count[j][k][g] - MC[j][k][g]*MC[j][k][g+10]*Size)/count[j][k][g]);
                }
            }
            
            if(Test[j][k] > 5){
                    aux = 1;
                    while(J[aux] != j && aux < s){
                        aux++;
                    }
                    if(aux == s){
                        J[s] = j;
                        s++;
                    }
                    aux = 1;
                    while(J[aux] != k && aux < s){
                        aux++;
                    }
                    if(aux == s){
                        J[s] = k;
                        s++;
                    }
            }*/
            
            //new 12-7
            if(Ht[j][k] > 0.9*Size || Ht[j][k] < 0.1*Size){
                aux = 1;
                while(J[aux] != j && aux < s){
                    aux++;
                }
                if(aux == s){
                    J[s] = j;
                    s++;
                }
                
                aux = 1;
                while(J[aux] != k && aux < s){
                    aux++;
                }
                if(aux == s){
                    J[s] = k;
                    s++;
                }
            }
        }
    }
    
    //new 12-7
    J[0] = s-1;
    return J;
}




//** not used **//
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}
