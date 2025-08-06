//
//  main.cpp
//  Random Forest
//
//  Created by zhang on 2017/7/16.
//  Copyright 20171. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector> 

#include "stdlib.h"
#include "RandomForest.hpp"
#include "RandomForest.cpp"
using namespace std;


int main(){
//** STEP 1 open data **//
    
    std::ifstream file;
    std::ofstream outfile6;
    bool real_data = false;

    /* load the dataset */
    if(real_data == true){
        file.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/real_data/Higgs/Higgs_train_1.txt");
    }
    else{
        //file.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/100dim/10j+90m/raw.txt");
        file.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/train_data.txt");
        
    }
    //OctTT!/testdata/dim3/cos+unif+beta/spar3.txt

    /*initialize the variable */
    //double data[Size][dim] = {-1}; // data matrix
    vector< vector<double> > data(Size, vector<double>(dim));
    vector<double> marginal_max(dim), marginal_min(dim); // boundaries of the space
    vector<double> margins(2*dim); // margins of each feature
    int i,j,k; // variables
    double volume = 1; // the volume
    vector<int> output_important_features; // a list of correlated features (important)
    for(i = 0; i < 2*dim;i++){
        if(i%2 == 0){margins[i] = 0;}
        else{margins[i] = 1;}
    }
    std::cout << "Finish initialization" <<std::endl;
    
    /////*****  generate data (from R)  *****/////
    //Modified at 10/11/2018 By Yiliang
    //Test the CDF of max & min point
    //double CDF[dim][2] = {-1};
    
    for(i = 0; i < Size; i++){
        for(j = 0; j< dim; j++){
            file >> data[i][j]; 
            if(i==0){marginal_min[j] = data[i][j];marginal_max[j] = data[i][j];}
            if(data[i][j] < marginal_min[j]){
                marginal_min[j] = data[i][j];
            }
            else if (data[i][j] > marginal_max[j]){
                marginal_max[j] = data[i][j];
            }
        }
    }
    std::cout << "Finish data loading" <<std::endl;
    file.close();
    
    /* calculate the volume of the space */
    for(i = 0; i < 2*dim;i++){
        if(i%2 == 0){margins[i] = marginal_min[i/2];}
        else{
            margins[i] = marginal_max[i/2];
            volume = volume*(margins[i]-margins[i-1]);
        }
    }
    
    /*J = covariance(data);
    K[0] = J[0];
    for(i = 1;i < K[0]+1;i++){
        K[i] = J[i];
    }*/
    
    
//** STEP 2 marginal estimation & Copula transformation **//
    // notice the random forest struct
    //TNode* RTree;
    //RTree = (TNode*)malloc(sizeof(TNode)*1);
    //double sample[Size][dim];
    //double sample_t[Size][dim];
    vector< vector<double> > sample(Size, vector<double>(dim));
    //vector< vector<double> > sample_ts(Size, vector<double>(dim));
    vector<double> Dat(Size);
    vector<double> datt;
    double mar;
    
    TNode* *TList = new TNode*[N_tree];
    //TNode* FList[dim][NN];  // = new (TNode* [10][dim]);
    vector< vector<TNode*> > FList(dim, vector<TNode*>(NN));
    
    for(i = 0;i < dim; i++){
        for(j = 0; j < NN; j++){
            FList[i][j] = new TNode;
        }
    }
    
    /* data preprocessing */
    TNode* Forest[NN];
    for(i = 0; i < dim; i++){
        for(j = 0; j < Size; j++){
            Dat[j] = data[j][i];
            //if(i==10){
            //    //outfile6 << Dat[j] << std::endl;
            //}
        }

        for(j = 0; j<NN;j++){
            Forest[j] = FList[i][j];
        }
        /* copula transformation */
        datt = Copula(Forest, Dat, marginal_max[i], marginal_min[i]);
        for(j = 0; j < Size; j++){
            sample[j][i] = datt[j];
            //datt++;
        }
        for(j = 0; j<NN;j++){
            FList[i][j] = Forest[j];
        }

        /* record the true density in the first dimension */
        //if(i == 1){
        //    for(k = 0;k<10000;k++){
        //        mar = (k*0.0001)*(marginal_max[1]-marginal_min[1]) + marginal_min[1];
        //        mar = get_density(mar, Forest);
        //        outfile6 << mar << std::endl;
        //    }
        //}
    }
    std::cout << "Finish copula transformation" <<std::endl;
    
/* select important features */
    output_important_features = covariance(sample);
    vector<int> important_features(dim+1);
    important_features[0] = output_important_features[0];
    for(i = 1;i < important_features[0]+1;i++){
        important_features[i] = output_important_features[i];
    }
    
/*
    for(i = 0; i < dim; i++){
        for(j = 0; j < Size; j++){
            Dat[j] = data[j][i];
        }
        //datt = Copula(Dat, marginal_max[i], marginal_min[i]);
        for(j = 0; j < Size; j++){
            sample[j][i] = datt[j];
            //datt++;
        }
    }
 */



//** STEP 3 Max Random Forest **//
    /* what is this? */
    for(i = 0; i < 2*dim; i++){
        if(i%2 == 0){margins[i] = 0;}
        else{
            margins[i] = 1;
            volume = 1;
        }
    }
    
    /* loop to grow each density trees */
    int initial_depth = 0, id_tree = 0;
    for(id_tree=0; id_tree<N_tree; id_tree++){
        TList[id_tree] = new TNode;
        RandomForest(TList[id_tree], sample, Size, margins, volume, initial_depth, 1, important_features);
    }
    
    /* update the volume */
    for(i=0; i<dim; i++){
        volume = volume*(marginal_max[i] - marginal_min[i]);
    }
    std::cout << "Finish model construction" <<std::endl;
    



//** STEP 4 calculate the error **//
    if (real_data == false){
        double H = Hellinger_Error(TList, FList, marginal_min, marginal_max, volume);
        std::cout << " H-error: " << H << std::endl;
        //Marginal_error(X,Y);
        //std::cout << "marginal bin: " << Cc << std::endl;
    }
    else{
        /* 1/30/2019 */
        std::ifstream outfile4;
        outfile4.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/real_data/Higgs/Higgs_test.txt");
        std::ofstream outfile5;
        outfile5.open("/Users/leoz/Desktop/Research/Max-Random Forest/Code/MRF/Data/real_data/Higgs/Higgs_result_1.txt");
        double g, gde;
        vector<double> dat(dim);
        
        for(i = 0; i < 1000; i++){
            g = 1;
            for(j = 0; j < dim; j++){
                //input the KL data
                outfile4 >> dat[j];
            }
            
            TNode* ForestM[NN];
            double marginal_pdf;
            //anti-copula
            for(j = 0; j<dim;j++){
                for(k = 0; k < NN; k++){
                    ForestM[k] = FList[j][k];
                }
                
                //to prevend ttt be 0
                marginal_pdf = (get_density(dat[j], ForestM) + 0.0000)/1.0000;
                if(marginal_pdf == 0){
                    marginal_pdf = 0.001;
                }
                
                //if KL data is outside [min,max]
                if(dat[j] > marginal_max[j] || dat[j] < marginal_min[j]){
                    g = 0;
                    break;
                }
                
                //calculate the coordinate after copula
                dat[j] = get_cord(dat[j], ForestM);
                g = g*marginal_pdf;
            }
            
            //take the copula density
            gde = Density_Estimate(dat, TList);
            
            //density equals to copula pdf multiplies n marginal pdfs
            g = g*gde;
            outfile5 << g << std::endl;
        }
        outfile4.close();
        outfile5.close();
    }
    std::cout << "Finish evaluation" <<std::endl;
    
//** STEP 5 delete the forest and release the memory **//
    Eliminate_Forest(TList);
    std::cout << "Finish model deletion" <<std::endl;

    return 0;
}
