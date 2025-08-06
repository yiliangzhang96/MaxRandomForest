
#ifndef LETTER_H
#define	LETTER_H

#include "general_utils.h"
#include "output.h"
#include "SISfunctions.h"
#include "buildpathcttable.h"
#include "sampling.h"
#include "data_store.h"

// end_col_normalization=FALSE for classification problem
inline vector<vector<double> > read_data(string filename, bool end_col_normalization) {
    int dim = 0;
   cerr<<"start read data!"<<endl;
    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

    while (!infile.eof()) {

        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
        data.push_back(d);
    }

    int k= dim;
    if(!end_col_normalization){ k=k-1;}

    vector<double> marginalmax(k, -100000000);
    vector<double> marginalmin(k,  1000000000);
    for(int i=0; i<k; i++){
        for (int j=0; j<(int)data.size(); j++){
            if(data[j][i]>marginalmax[i])  marginalmax[i]=data[j][i];
            if(data[j][i]<marginalmin[i])  marginalmin[i]=data[j][i];

        }
        cerr<<"max="<<marginalmax[i]<<"; min="<<marginalmin[i]<<endl;
    }
    for (int i = 0; i < k; i++) {
        if (marginalmin[i] == marginalmax[i]) {
            cout << "only one value in dim=" << i << endl;
            for (int j = 0; j < (int) data.size(); j++) {
                 data[j][i] = 0.45;
            }
            continue;
        }
        for (int j = 0; j < (int) data.size(); j++) {
            data[j][i] = (data[j][i] - marginalmin[i]) / (marginalmax[i] - marginalmin[i]);

        }
    }

    return data;
}


// end_col_normalization=FALSE for classification problem
inline vector<vector<double> > read_data(string filename,bool end_col_normalization, vector<double>& mmax, vector<double>& mmin) {
    int dim = 0;
   cerr<<"start read data!"<<endl;
    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();
    cerr<<"dim="<<dim<<endl;
    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

    while (!infile.eof()) {

        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
   //     cerr<<d[0]<<'\t';
        data.push_back(d);
    }

    int k= dim;
    if(!end_col_normalization){ k=k-1;}

    vector<double> marginalmax(k, -100000000);
    vector<double> marginalmin(k,  1000000000);
    for(int i=0; i<k; i++){
        for (int j=0; j<(int)data.size(); j++){
            if(data[j][i]>marginalmax[i])  marginalmax[i]=data[j][i];
            if(data[j][i]<marginalmin[i])  marginalmin[i]=data[j][i];

        }
        cerr<<"max="<<marginalmax[i]<<"; min="<<marginalmin[i]<<endl;
    }
    int sumsmall=0;
    for (int i = 0; i < k; i++) {
        if (marginalmin[i] == marginalmax[i]) {
            cout << "only one value in dim=" << i << endl;
            for (int j = 0; j < (int) data.size(); j++) {
                 data[j][i] = 0.45;
            }    
            continue;
        }
        for (int j = 0; j < (int) data.size(); j++) {
  //          data[j][i] = (data[j][i] - marginalmin[i] + 0.001) / (marginalmax[i] - marginalmin[i] + 0.002);
            data[j][i] = (data[j][i] - marginalmin[i]) / (marginalmax[i] - marginalmin[i]);
            if(data[j][i]<0) data[j][i]=0;
            if(data[j][i]>1) data[j][i]=1;

//            cerr<<data[j][i]<<'\t';
            if(data[j][i]<0.01) sumsmall+=1;

        }
    }
    cerr<<"sumsmall="<<sumsmall<<endl;
    mmax = marginalmax;
    mmin = marginalmin;

    return data;

}

inline vector<vector<double> > read_partition(string filename, bool end_line) {
    int dim = 0;

    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

    while (!infile.eof()) {

        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        bool good_data = true;
        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
            if (end_line) {
                if (i != (dim - 1) && (d[i] > 1.0 || d[i] < 0)) {   //parentheses by thchiu
                    good_data = false;
                    break;
                }
            } else {
                if (d[i] > (1.0 - (1.0 / 2147483648)) || d[i] < (0 + (1.0 / 2147483648))) {
                    good_data = false;
                    break;
                }
            }
        }
        if(good_data){
            data.push_back(d);
        }else{
            cerr << "Warning: Data out of range("<< dim <<"): " << line << '\n';
        }

    }
    return data;
}

inline vector<vector<double> > read_data(string filename) {
    int dim = 0;
   cerr<<"start read data!"<<endl;
    ifstream infile(filename.c_str());

    if (!infile.is_open()) {
        cerr << "ERROR: Could not open " << filename << '\n';
        exit(1);
    }

    string line;
    getline(infile, line);

    trim2(line);

    if(line.length() == 0){
        cerr << "ERROR: Empty file: " << filename << '\n';
        exit(1);
    }

    vector<string> line_list = split(line);

    dim = line_list.size();

    infile.close();

    infile.open(filename.c_str());

    vector<vector<double> > data;

    while (!infile.eof()) {

        getline(infile, line);
        trim2(line);

        if (line.length() == 0) continue;

        vector<string> ll = split(line);

        if ((int) ll.size() != dim) {
            cerr << "ERROR: Bad line dim: " << line << '\n';
            exit(2);
        }

        vector<double> d;
        d.resize(dim);

        for (int i = 0; i < dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
        data.push_back(d);
    }

    return data;
}



#endif







