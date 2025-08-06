/*
 Bayesian Sequential Partitioning
 2013-1  by LuoLu
 */

// hello
//#include <iostream>
#include "tree.h"
#include <math.h>
#include <cmath>
#include "readdata.h"
#include "examples.h"
#include "SISfunctions.h"
#include "timer_ed520.h"    //+ by ed520
TIMER timer;                //+ by ed520

int testF (vector< map< OnePartition , double, CompairSRegs> >& pathctmap); //don't been called anywhere, thchiu
int BSP_C(vector<string> params);
int BSP_NC(vector<string> params);
void print_usage_and_exit();
int density(vector<string> & params);
int BSP_C_forden(vector<string> params);
int BSP_tree_C(vector<string> params);     //+ by 520;
int BSP_tree_NC(vector<string> params);     //+ by 520;

void test_partition_generation( Partition &mytestpartition ,  Region_Node *  mytestregion_a    );

int countfrompartition(vector<string> params);

int main(int argc, char** argv) {

 //   cerr << "**" + c::PROG_NAME << '\n';
 //   cerr << "Version: " << c::BUILD << "\n";
    cerr << "Luo Lu\n";
    cerr << "http://www.stanford.edu/group/wonglab" << "\n";
    cerr << "\n";

    string mode = "";
    string temp = "";
    int error_num = 0;
    vector<string> params;

    if (argc < 2) {
        print_usage_and_exit();
    }

    mode = argv[1];
    for (int i = 2; i < argc; i++) {
                temp = argv[i];
                params.push_back(temp);
    }

    timer.Begin();
    if (mode == "BSP_NC") {
         error_num = BSP_NC(params);
    } else if(mode == "BSP_C"){
         error_num= BSP_C_forden(params);
//    } else if (mode == "hell_dist") {
 //        error_num = hell_dist(params);
    } else if (mode == "density") {
         error_num = density(params);
    } else if(mode == "count"){
        error_num = countfrompartition(params);
    } else if(mode == "BSP_tree_C"){
        error_num=BSP_tree_C(params);
    } else if(mode == "BSP_tree_NC"){
        error_num=BSP_tree_NC(params);
    }
    else{
        print_usage_and_exit();
    }
    timer.End();
    timer.RunTimeProfile();

    //densities_from_partition(<#vector<double> &testdata#>, <#OnePartition_data &P#>, <#int n#>);
    return error_num;

}

int BSP_C(vector<string> params){
    string usage_text = "Usage: BSP_C <data_file> <output_file> <level1=1000> <level2=200> \n data_file -- One sample each row; MAP partitions output to output_file. Log to STDOUT; level1 -- number of SIS levels for joint distribution; level2 -- number of SIS levels for marginal distribution";

    if (params.size() > 4 || params.size() < 2) {
        cerr << usage_text << endl;
        return 3;
    }
    vector<double> mmax, mmin;
    
    string ofilename = params[1];
    parameters  p;
    p.data = read_data(params[0], true, mmax, mmin);
    p.dim=(int)p.data[0].size();
 //   p.levels=strTo<int>(params[1]);
    if(params.size()==4) p.levels= strTo<int>(params[2]);
    else p.levels = 1000;
    p.maxpercentage = 0.9;
    p.n = 200;
    p.resampling=2;
    p.samplesize= (int)p.data.size();;
    p.steps=4;
    p.data1D = new double[p.dim*p.samplesize]; //+ by thchiu
    p.pt_start = new int [p.samplesize * p.n];//+ by thchiu
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

    int dim = p.dim;
    int samplesize = p.samplesize;
    parameters para;
    para.maxpercentage = 0.9;
    para.samplesize = samplesize;
    para.steps = p.steps;
    para.dim = 1;
    para.n = 200;
    para.data1D = new double[para.dim*para.samplesize]; //+ by thchiu
    para.pt_start = new int [para.samplesize * para.n];//+ by thchiu

    if(params.size()==4)    para.levels = strTo<int>(params[3]);
    else para.levels=200;
    para.resampling = 2;
    //vector<vector<double> >  dataF= p.data;   //commented by thchiu
    vector<OnePartition_data> marginalP;

     //vector<OnePartition_data>  bestlevelPs;  //commented by thchiu
     OnePartition_data bestlevelPs; //added by thchiu
     //int bestind; //commented by thchiu

    for(int d=0; d<dim; d++){
        cout<<"dim="<<d<<endl<<endl;
        para.data.clear();
        vector<double> onepiece(1,0);
        vector<double> onevariable(samplesize, 0);
        for(int i=0; i<samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);
            para.data1D[i*para.dim]=p.data[i][d];   //+ by thchiu
            for(int j=0; j<para.n; j++)             //+ by thchiu
                para.pt_start[j*para.samplesize+i] = i*para.dim;
        }


        /*SIS(para, true, bestlevelPs, bestind);
        marginalP.push_back( bestlevelPs[bestind]);*/ //commented by thchiu
        SIS(para, true, bestlevelPs);   //thchiu
        marginalP.push_back( bestlevelPs);//thchiu
        print_partition(marginalP[d]);

  //      print_partition(marginalP[d]);
        vector<double> newx = Ftransform(onevariable, marginalP[d]);

        for(int i=0; i<10; i++) cerr<<onevariable[i]<<" ";
        cerr<<endl;

        for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
        cerr<<endl;

        vector<double> recoveronevariable = inv_Ftransform(newx, marginalP[d]);

        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]<<" ";
        cerr<<endl;

        //for(int i=0; i<samplesize; i++)        //dataF[i][d]=newx[i];   //commented by thchiu
        //    p.data[i][d]=newx[i];

        for(int i=0; i<samplesize; i++){        //+ by thchiu
            p.data1D[i*p.dim+d] = newx[i];
        }
    }
     /*para.dim = dim;      //commented by thchiu
     para.data.clear();
     para.data=dataF;
     para.levels = p.levels;
     para.n = p.n;
     para.steps = p.steps;
     para.resampling = p.resampling;*/

     //SIS(para, true, bestlevelPs, bestind);   //commented by thchiu
     SIS(p, true, bestlevelPs);    //added by thchiu

     ///////// Output file is commented by ed520  ////////////
#define OUTPUT
#ifdef OUTPUT
    //     print_Finvtransform(P, marginalP);
    for (int d = 0; d < dim; d++) {
        string temp= ofilename+'_'+toStr<int>(d)+".txt";
        ofstream outfile (temp.c_str());
        // fprintf(pFile, "d= %d\n", d);
        print_partition(outfile, marginalP[d]);
        outfile.close();
    }
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());

    //OnePartition_data op = bestlevelPs[bestind];  //commented by thchiu
    OnePartition_data &op = bestlevelPs; //added by thchiu

     //print_partition(outfile, bestlevelPs[bestind]);  //commented by thchiu
    print_partition(outfile, bestlevelPs);  //added by thchiu

    outfile <<endl << "After transform back" << endl;
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < dim; d++) {
            range = convert_ranges(op.sregs[i].reg_code[d]);
            double ub,lb;
            ub = inv_Ftransform(range.second, marginalP[d]);
            ub = ub * (mmax[d] - mmin[d])+ mmin[d];
            lb =  inv_Ftransform(range.first, marginalP[d]);
            lb = lb * (mmax[d] - mmin[d])+ mmin[d];
            outfile <<lb << " " << ub <<" ";
        }
        outfile<< op.sregs[i].num<<endl;

    }

    outfile.close();

 // output the overall regions without copula

    temp= ofilename+"_Big_total.txt";
    ofstream outfile2 (temp.c_str());
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
          pair<double, double> range;
          vector<usint_mask> reg_codes= op.sregs[i].reg_code;
          vector<pair<double, double> > regs;
          for (int d = 0; d < dim; d++) {

              range = convert_ranges(reg_codes[d]);
              pair<double, double> bounds;
              bounds.first = inv_Ftransform(range.first, marginalP[d]);
              bounds.second =  inv_Ftransform(range.second, marginalP[d]);
              regs.push_back(bounds);
          }

          recover_fullregions(regs, marginalP,mmax, mmin, outfile2);

    }

    outfile2.close();
#endif

    return 0;

}

// output the partition after normalization
int BSP_NC(vector<string> params){

//    string usage_text = "Usage: BSP_NC <data_file> <output_file> <level1=1000>\n"
 //           + " data_file -- One sample each row\n"
//            + " MAP partitions output to output_file. Log to STDOUT \n"
//            + " level -- number of SIS levels for joint distribution";


    if (params.size() > 4 || params.size()<2) {
//        cerr << usage_text << endl;
        cerr<<"wrong input!\n";
        return 3;
    }
    string ofilename = params[1];
    parameters p;
    vector<double> mmax, mmin;
    p.data = read_data(params[0],true,mmax, mmin);
    //p.testdata = read_data(params[2]);
    cerr<<p.data[0][0]<<" "<<p.data[1][0]<<" "<<p.data[2][0]<<" "<<p.data[3][0]<<endl;

    p.dim=(int)p.data[0].size();
    if(params.size()>2)   p.levels=strTo<int>(params[2]);
    else
    {
     p.levels = 1000;
    }
    p.maxpercentage = 0.9;
    //if(params.size()>3)   p.n=strTo<int>(params[3]);   //+ by ed520 else
     p.n = 200;
    p.resampling=2;
    p.samplesize= (int)p.data.size();;
    p.steps=4;

    p.data1D = new double[p.dim*p.samplesize];  // new

    for(int d = 0;d<p.dim;d++){
        for(int i = 0;i<p.samplesize;i++){
            p.data1D[i*p.dim+d] = p.data[i][d];

        }
    }


    p.pt_start = new int [p.samplesize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

     /*vector<OnePartition_data>  Ps;    //commented by thchiu
    int bestind;        //index for choosing the best partition in Ps, thchiu
    */
    OnePartition_data Ps;   //added by thchiu
    vector<OnePartition_data> P;
    //mmax[0] = 1;mmax[1] = 1;mmin[0] = 0;mmin[1] = 0;  //added by Yiliang
    //SIS(p, true, Ps, bestind);    //commented by thchiu
    SIS(p, true, Ps, ofilename, mmax, mmin);   //added by thchiu
    
    
    /////////***********   12-21   ***********//////////////
    /*
    P.push_back(Ps);
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    Output_distance( p, P,NULL, 0, pathctmap, false);
*/
    
//    print_partition(Ps);
    //ofstream outfile (ofilename.c_str()); //- by ed520
    string temp= ofilename+"_Big.txt";  //+ by ed520
    ofstream outfile (temp.c_str()); //+ by ed520
 //   print_partition(outfile,Ps[bestind]);

  //  outfile <<endl<<endl;


    //print_partition(outfile,Ps[bestind], mmax, mmin); //commented by thchiu
    print_partition(outfile,Ps, mmax, mmin);    //added by thchiu

    outfile.close();
    //+ by ed520
    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<p.dim<<endl;
    cout<<"Data_num = "<<p.samplesize<<endl;
    cout<<"# of Partition = "<<p.n<<endl;
    //+ end by ed520
    return 0;

}


int density(vector<string> & params){

   /*string usage_text = "Usage: \n"
            + " density <-c/-n> <Partitions> <sample_data> \n"
            + "       -c/-n -- c=use copula, n=no copula\n"
            + "  partitions -- partitions of a distribution\n"
            + " sample_data -- Each row one data point\n "
            + " Output the density.\n";
   */
    cerr << params[0] << endl;
    cerr << params[1] << endl;
    cerr << params[2] << endl;
    
    
    if (params.size() != 3) {
        cerr<<params.size()<<endl;
        cerr << "usage_text" << endl;
        return 3;
    }
    vector<double> mmax, mmin;
//    cerr<<params[1]<<'\n'<<params[2]<<''
    vector<vector<double> > test_data = read_data(params[1],true, mmax, mmin);  //data (params[0], true, mmax, mmin);
    
    //vector<vector<double> > mmax_min = read_data(params[3]); // file size: dim*2

    int test_N = (int)test_data.size();
    int dim    = (int)test_data[0].size();
    vector<double> dens(test_N, 0);

    for(int i=0; i<test_N; i++){
        for(int d=0; d<dim; d++){
       test_data[i][d] = (test_data[i][d] - mmin[d]) / (mmax[d] - mmin[d]);
        }
    }
    cerr << test_N << " data points in " << dim << " dimensions.\n";
    cerr << "some data: "<<test_data[10][1]<<" "<<test_data[100][20]<<endl; //
    if (params[0]=="-c"){
         vector<vector<vector<double> > > MarginalPs;
         for(int d=0; d<dim; d++){
               MarginalPs.push_back(read_data( params[2]+'_'+toStr<int>(d)+".txt"));
         }
         vector<vector<double> > P = read_data( params[2]+"_Big.txt");//after copula transback, after normalize
         cerr<<"P="<<P[0][2]<<" "<<P[1][5]<<endl;

         for (int i = 0; i < test_N; i++) {
             if(max(test_data[i])>=1 || min(test_data[i])<=0){ // outside cubic
               //  cout<<i<<'\t'<<dens[i]<<'\n';
                 cerr<<"0"<<endl; //
                 continue;
             }
             dens[i] = densities_from_partition(test_data[i], P);  //the 2*dim col of P should be density

            if(i<10) cerr<<"\n i= \t"<<dens[i]<<endl;
            for (int d = 0; d < dim; d++) {
                vector<double> onepiece;
                onepiece.push_back(test_data[i][d]);
                dens[i] *= densities_from_partition(onepiece, MarginalPs[d]);
                 if(i<10) cerr<<dens[i]<<'\t';
            }

        }
        
    }
    else  if (params[0]=="-n"){
        vector<vector<double> > P = read_data( params[2]);
         for (int i = 0; i < test_N; i++) {
            dens[i] = densities_from_partition(test_data[i], P);
        }
    }
    else {
        cerr << "usage_text" << endl;
        return 1;
    }

    string tempout= params[2]+"_den.txt";
    ofstream outfileout (tempout.c_str());
    for(int i=0; i<test_N; i++){
        for(int d=0; d<dim; d++){
             dens[i]/= (mmax[d] - mmin[d]);

        }
        outfileout <<dens[i]<<'\n';
    }

    outfileout.close();


    return 0;

}

int countfrompartition(vector<string> params){

//    string usage_text = "Usage: count dim <data_file> <partition_file><output number file> \n"
 //           + " data_file -- One sample each row\n";


    if (params.size() > 4 || params.size()<4) {
//        cerr << usage_text << endl;
        cerr<<"wrong input!\n";
        return 3;
    }
    int dim = atoi(params[0].c_str());
    cerr<<"dim="<<dim<<endl;
    string datafilename = params[1];
    string partitionfilename = params[2];
    string ofilename = params[3];
    string line;
    ifstream infile;

    infile.open(partitionfilename.c_str());
    vector<vector<double> > partition;

    while (!infile.eof()) {
        getline(infile, line);
        trim2(line);
        if (line.length() == 0) continue;
        vector<string> ll = split(line);
        vector<double> d;
        d.resize(2*dim);
        for (int i = 0; i < 2*dim; i++) {
            d[i] = strTo<double>(ll[i]);
        }
  //      cerr<<d[7]<<" ";
        partition.push_back(d);
    }
//    cerr<<endl<<partition[3][7]<<endl;
    infile.close();
    int regnum = (int)partition.size();
    infile.open(datafilename.c_str());
    vector<int> nums(regnum, 0);

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
        for(int m=0; m<regnum; m++){
            bool inreg=true;
            for(int j=0; j<dim; j++){
                if(d[j]> partition[m][2*j+1] || d[j]<= partition[m][2*j]){
                    inreg=false;
                    break;
                }
            }
            if(inreg) {
                nums[m]++;
                break;
            }

        }
    }
    infile.close();

    ofstream outfile (ofilename.c_str());

    for(int i=0; i<regnum; i++) outfile<<nums[i]<<" ";
    outfile<<endl;
    cerr<< endl << sum(nums) <<endl;
    return 0;
}

int BSP_C_forden(vector<string> params){
    string usage_text = "Usage: BSP_C <data_file> <output_file> <level1=1000> <level2=200> \n data_file -- One sample each row; MAP partitions output to output_file. Log to STDOUT; level1 -- number of SIS levels for joint distribution; level2 -- number of SIS levels for marginal distribution";

//*** initialization procedure (start) + by FJKu
    if (params.size() > 10 || params.size() < 2) {
        cerr << usage_text << endl;
        return 3;
    }
    string ofilename = params[1];
    parameters  p;
    vector<double> mmax, mmin;
    
    p.data = read_data(params[0], true, mmax, mmin);

    p.testdata = read_data(params[2]);


    string tempmaxmin= ofilename+"_maxmin.txt";
    ofstream outfilemaxmin (tempmaxmin.c_str());
    for(int d=0; d<(int)p.data[0].size();d++){
        outfilemaxmin << mmax[d]<<'\t'<<mmin[d]<<'\n';
    }
    outfilemaxmin.close();

    p.dim=(int)p.data[0].size();
 //   p.levels=strTo<int>(params[1]);
    if(params.size()==4) p.levels= strTo<int>(params[2]);
    else p.levels = 1000;
    p.maxpercentage = 0.9;
    if(params.size()==4) p.n = strTo<int>(params[3]); //+ by ed520
    else p.n = 200;
    p.resampling=2;
    p.samplesize= (int)p.data.size();;
    p.steps=4;

    p.data1D = new double[p.dim*p.samplesize];  // new
    p.pt_start = new int [p.samplesize * p.n];
    for(int i=0; i<p.n; i++){
        for(int j=0; j<p.samplesize; j++){
            p.pt_start[i*p.samplesize+j] = j*p.dim;
        }
    }

    //    bool pathct=false;
    //    bool discrete=true;

    //    double smoothneighbordist=0.001;

    int dim = p.dim;
    int samplesize = p.samplesize;
    parameters para;
    para.maxpercentage = 0.9;
    para.samplesize = samplesize;
    para.steps = p.steps;
    para.dim = 1;
    //para.n = 200; //- by ed520
    para.n = p.n; //+ by ed520
    para.data1D = new double[para.dim*para.samplesize]; // new
    para.pt_start = new int[para.samplesize*para.n];

    //if(params.size()==4)    para.levels = strTo<int>(params[3]);  //- by ed520
    //else para.levels=200; //- by ed520
    para.levels = p.levels; //+ by ed520
    para.resampling = 2;
    //    vector<vector<double> >  dataF= p.data;
    vector<OnePartition_data> marginalP;

    //   vector<OnePartition_data>  bestlevelPs;
    OnePartition_data  bestlevelPs;
    //    int bestind;

//*** initialization procedure (end) + by FJKu

    string temptrans = ofilename+"_trans.txt";
    ofstream outfiletrans (temptrans.c_str());

//*** Copula procedure (start) + by FJKu

    for(int d=0; d<dim; d++){ // loop D-times, because we have to do 1-D BSP for D dimensions  + by FJKu
        cerr<<"dim="<<d<<endl<<endl;
        para.data.clear();
        vector<double> onepiece(1,0);
        vector<double> onevariable(samplesize, 0);
        for(int i=0; i<samplesize; i++)   {
            onepiece[0]=p.data[i][d];
            onevariable[i]= p.data[i][d];
            para.data.push_back(onepiece);

            para.data1D[i*para.dim]=p.data[i][d]; // new
            for(int j=0; j<para.n; j++)
                para.pt_start[j*para.samplesize+i] = i*para.dim;
        }

        SIS(para, true, bestlevelPs);

        marginalP.push_back( bestlevelPs);
        print_partition(marginalP[d]);

        //    print_partition(bestlevelPs);
        vector<double> newx = Ftransform(onevariable, marginalP[d]); // disperse data + by FJKu

        for(int i=0; i<10; i++) cerr<<onevariable[i]<<" ";
        cerr<<endl;

        for(int i=0; i<10; i++) cerr<<newx[i]<<" ";
        cerr<<endl;

        vector<double> recoveronevariable = inv_Ftransform(newx, marginalP[d]);

        for(int i=0; i<10; i++) cerr<<recoveronevariable[i]<<" ";
        cerr<<endl;


        for(int i=0; i<samplesize; i++){
            //dataF[i][d]=newx[i];            //new
            //p.data[i][d] = newx[i];       //- by ed520, we use below
            p.data1D[i*p.dim+d] = newx[i];  //+ by ed520
            outfiletrans << newx[i]<<'\t';
        }
        outfiletrans<<'\n';
    }
//*** Copula procedure (end) + by FJKu

    outfiletrans.close();
    /*    para.dim = dim;
          para.data.clear();
          para.data=dataF;
          para.levels = p.levels;
          para.n = p.n;
          para.steps = p.steps;
          para.resampling = p.resampling;*/           //new
    

//*** Main SIS procedure (start) + by FJKu
	// Start multi-dim BSP + by FJKu
    SIS(p, true, bestlevelPs);
//*** Main SIS procedure (end) + by FJKu
    

    /////////***********   12-21   ***********//////////////
    
    vector< map< OnePartition , double, CompairSRegs> > pathctmap;
    Output_distance( p, marginalP,NULL, 0, pathctmap, false);

//*** Below are not so important stuff for output + by FJKu
    for (int d = 0; d < dim; d++) {
        string temp= ofilename+'_'+toStr<int>(d)+".txt";
        ofstream outfile (temp.c_str());

        // fprintf(pFile, "d= %d\n", d);
        print_partition(outfile, marginalP[d]);
        outfile.close();
    }
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());
    //   print_partition(outfile, bestlevelPs[bestind]);
    //   OnePartition_data op = bestlevelPs[bestind];
    OnePartition_data &op = bestlevelPs;


    //    outfile <<endl << "After transform back" << endl;
    for (int i = 0; i < (int) (int) op.sregs.size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < dim; d++) {
            range = convert_ranges(op.sregs[i].reg_code[d]);
            double ub,lb;
            ub = inv_Ftransform(range.second, marginalP[d]);
            //     ub = ub * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            lb =  inv_Ftransform(range.first, marginalP[d]);
            //     lb = lb * (mmax[d] - mmin[d] + 0.002)+ mmin[d] - 0.001;
            //         outfile<<range.first<<" "<<range.second<<" ";
            outfile <<lb << " " << ub <<" ";
        }
        outfile<< (double)op.sregs[i].num/(double)samplesize/ exp(lprod_usint_mask(op.sregs[i].reg_code)) <<" "<< op.sregs[i].num<<endl;

    }

    outfile.close();
    //+ by ed520
    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<p.dim<<endl;
    cout<<"Data_num = "<<p.samplesize<<endl;
    cout<<"# of Partition = "<<p.n<<endl;
    //+ end by ed520
    return 0;

}
void print_usage_and_exit() {
    //    cerr << "Usage: " + c::PROG_NAME + " <BSP>" << "\n";
    cerr << "Options:" << "\n";
    cerr << "-== Density Estimation ==-" << '\n';
    cerr << "  BSP_C       -- MAP partitions for marginal distribution and MAP partition for joint distribution after copula transformation" << "\n";
    cerr << "  BSP_NC      -- MAP partition without copula" << "\n";
    cerr << "\n";
    //    cerr << "-== Other tools ==-" << '\n';
    //    cerr << "  hell_dist   -- Compute sample Hellinger distance from a known density" << "\n";
    //    cerr << "  density     -- Get the density at particular points" << "\n";
    exit(2);
}



int BSP_tree_C(vector<string> params ){
    string wrong_message = "ED520 tree version : wrong input";
    if(params.size()> 4 || params.size() < 2)
    {
        cout<<wrong_message<<endl;
        return 3;
    }
    else{
        cerr << "BSP_Tree is programmed by Parallel Computing System Lab (ED520) in NCTU\n";
        cerr << "https://sites.google.com/a/g2.nctu.edu.tw/parallel-computing-lab/" << "\n";
        cerr << "\n";
    }
    tree BspTree;
    BspTree.discrete=true;
    string ofilename = params[1] ;
    vector<double> mmax, mmin;    
    BspTree.Data = read_data(params[0], true, mmax, mmin);
    BspTree.Dim = BspTree.Data[0].size();
    if(params.size()==4)
        BspTree.Levels= strTo<int>(params[2]);
    else 
        BspTree.Levels = 1000;  // default # of cut
    if(params.size()==4)
        BspTree.n = strTo<int>(params[3]); 
    else 
        BspTree.n = 200;       // default # of partition
    
    string tempmaxmin= ofilename+"_maxmin.txt";
    ofstream outfilemaxmin (tempmaxmin.c_str());
    for(int d=0; d<BspTree.Dim   ;d++){
        outfilemaxmin << mmax[d]<<'\t'<<mmin[d]<<'\n';
    }
    outfilemaxmin.close();
    
    BspTree.tree_nonpara_init( );
    // copula initialization : 
    tree BspTree_para;
    BspTree_para.discrete=true;
    BspTree_para.maxpercentage=BspTree.maxpercentage; //some parameter for statistics
    BspTree_para.Levels=BspTree.Levels;     
    BspTree_para.n=BspTree.n;          
    BspTree_para.Resampling=BspTree.Resampling;
    BspTree_para.Samplesize=BspTree.Samplesize;   
    BspTree_para.Steps=BspTree.Steps;
    BspTree_para.beta=BspTree.beta;
    // end of copula parameter !!! 
    string temptrans = ofilename+"_trans.txt";
    ofstream outfiletrans (temptrans.c_str());
    // for final output file 
#ifdef OUTPUT_FILE
    vector<OnePartition_data> marginalP;
    vector<pair <double , double> > temp_intervals;
    vector <vector<pair <double , double> > > marginal_intervals;
#endif
    vector <double> onevariable(BspTree.Samplesize,0);
    BspTree_para.para_share_init();
    //start Copula
    for(int i=0;i<BspTree.Dim;i++){
        cout<<"dim is "<<i<<endl;
        BspTree_para.tree_para_init(  );
        
        for(int j=0 ; j < BspTree_para.Samplesize ; j++ ){
            BspTree_para.Data_array[j]= BspTree.Data_array[j*BspTree.Dim+i];
            onevariable[j]=BspTree.Data_array[j*BspTree.Dim+i]; //  Data[j][i]; 
        }
        for(int j =0 ; j <BspTree.Samplesize; ++j ){
            BspTree_para.Point_index[j]= BspTree_para.Dim*j;    //  create all index first 
        }

        BspTree_para.SIS_tree();
       
        BspTree_para.List_index=BspTree_para.List_index_best;
#ifdef DEBUG
        BspTree_para.print_partition1_t(  BspTree_para.P_bestlevel );
#endif
        vector<double> newx = BspTree_para.Ftransform_t(onevariable,  BspTree_para.P_bestlevel);

        for(int j=0; j<10; j++) cerr<<onevariable[j]<<" ";
        cerr<<endl;

        for(int j=0; j<10; j++) cerr<<newx[j]<<" ";
        cerr<<endl;
#ifdef DEBUG
        vector<double> recoveronevariable = BspTree_para.inv_Ftransform_t(newx, BspTree_para.P_bestlevel );

        for(int j=0; j<10; j++) cerr<<recoveronevariable[j]<<" ";
        cerr<<endl;
#endif
        for(int j=0; j<BspTree.Samplesize; j++){
            BspTree.Data_array[j*BspTree.Dim+i] = newx[j];  //+ by ed520  need to be replaced 
            outfiletrans << newx[j]<<'\t';
        }
        outfiletrans<<'\n';
#ifdef OUTPUT_FILE
        string temp= ofilename+'_'+toStr<int>(i)+".txt";
        ofstream outfile (temp.c_str());
        BspTree_para.print_partition2_t( outfile, BspTree_para.P_bestlevel );
        temp_intervals = BspTree_para.onedimpartitiontointerval_t_copy( BspTree_para.P_bestlevel );
        marginal_intervals.push_back(temp_intervals);
        outfile.close();
#endif
        for(unsigned del=0; del < BspTree_para.Bsp_Tree.size(); del++ )
        {
            BspTree_para.Bsp_Tree[del]->destroy(); 
        }
        cout<<" Finish on Dim of copula Dim is  "<<i<<"  \n\n";
    }
    //end Copula
    BspTree_para.tree_para_destroy();
    outfiletrans.close();
    
    //start Non-copula
    BspTree.SIS_tree();
    BspTree.List_index=BspTree.List_index_best;
    
    // for output file 
#ifdef OUTPUT_FILE
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());
    Partition &op = BspTree.P_bestlevel;

    for (unsigned i = 0; i < op.region[BspTree.List_index].size(); i++) {
        pair<double, double> range;
        for (int d = 0; d < BspTree.Dim; d++) {
            range = convert_ranges(op.region[BspTree.List_index][i]->Region_mask[d]);
            double ub,lb;
            ub = BspTree.inv_Ftransform_for_output(range.second, marginal_intervals[d] /* marginalP[d]*/  );
            lb = BspTree.inv_Ftransform_for_output(range.first,  marginal_intervals[d] /* marginalP[d] */ );
            outfile <<lb << " " << ub <<" ";
        }
        outfile<< (double)op.region[BspTree.List_index][i]->length/(double)BspTree.Samplesize/ exp(lprod_usint_mask(op.region[BspTree.List_index][i]->Region_mask)) <<" "<< op.region[BspTree.List_index][i]->length<<endl;
    }
    outfile.close();
#endif

    // end of output file 
    // end of process
    // end of the showing !! 

    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<BspTree.Dim<<endl;
    cout<<"Data_num = "<<BspTree.Samplesize<<endl;
    cout<<"# of Partition = "<<BspTree.n<<endl;   
#ifdef PROFILE
    cout<<"====== Profile Info ======"<<endl;
    cout<<"NC CPU COUNT = "<<timer.CPU_COUNT<<endl;
    cout<<"NC HIT COUNT = "<<timer.HIT_COUNT<<endl;
    cout<<"Memory Usage = "<<BspTree.Last_point_index*4/1024/1024<<" MB"<<endl;
    cout<<"Total Memory = "<<(long)(BspTree.n*BspTree.Samplesize*SCALE)*4/1024/1024<<" MB"<<endl;
#endif
    BspTree.SaveTree();
    return 0;

}

int BSP_tree_NC(vector<string> params ){
    string wrong_message = "ED520 tree version : wrong input";
    if(params.size()> 4 || params.size() < 2)
    {
        cout<<wrong_message<<endl;
        return 3;
    }
    else{
        cerr << "BSP_Tree is programmed by Parallel Computing System Lab (ED520) in NCTU\n";
        cerr << "https://sites.google.com/a/g2.nctu.edu.tw/parallel-computing-lab/" << "\n";
        cerr << "\n";
    }
    tree BspTree;
    BspTree.discrete=true;
    string ofilename = params[1] ;
    vector<double> mmax, mmin;
    BspTree.Data = read_data(params[0], true, mmax, mmin);
    BspTree.Dim = BspTree.Data[0].size();
    if(params.size()==4)
        BspTree.Levels= strTo<int>(params[2]);
    else
        BspTree.Levels = 1000;  // default # of cut
    if(params.size()==4)
        BspTree.n = strTo<int>(params[3]);
    else
        BspTree.n = 200;       // default # of partition
    
    BspTree.tree_nonpara_init( );

    BspTree.SIS_tree();
    BspTree.List_index=BspTree.List_index_best;

#ifdef OUTPUT_FILE
    string temp= ofilename+"_Big.txt";
    ofstream outfile (temp.c_str());

    BspTree.print_partition3_t(outfile, BspTree.P_bestlevel, mmax, mmin);

    outfile.close();
#endif

    cout<<"====== Bench Info ======"<<endl;
    cout<<"Dimension = "<<BspTree.Dim<<endl;
    cout<<"Data_num = "<<BspTree.Samplesize<<endl;
    cout<<"# of Partition = "<<BspTree.n<<endl;
#ifdef PROFILE
    cout<<"====== Profile Info ======"<<endl;
    cout<<"NC CPU COUNT = "<<timer.CPU_COUNT<<endl;
    cout<<"NC HIT COUNT = "<<timer.HIT_COUNT<<endl;
    cout<<"Memory Usage = "<<BspTree.Last_point_index*4/1024/1024<<" MB"<<endl;
    cout<<"Total Memory = "<<(long)(BspTree.n*BspTree.Samplesize*SCALE)*4/1024/1024<<" MB"<<endl;
#endif
    BspTree.SaveTree();
    return 0;
}

void test_partition_generation( Partition &mytestpartition ,  Region_Node *  mytestregion_a    ){

    int paradim=1;
    for( int ii=0 ; ii <10 ; ii ++){
        mytestregion_a[ii].length = ii;
        //mytestregion_b[i].length = 0;
        if(ii==0){
            mytestregion_a[ii].um.mask=0;
            mytestregion_a[ii].um.x=0;
            for(int j= 0 ; j <paradim ; j++)
                mytestregion_a[ii].Region_mask.push_back(mytestregion_a[ii].um);
        }
        else{
            for(int j = 0 ; j < paradim ; j ++){
                mytestregion_a[ii].Region_mask=mytestregion_a[ii-1].Region_mask;
                mytestregion_a[ii].Region_mask[j].x = mytestregion_a[ii-1].Region_mask[j].x<<1;
                mytestregion_a[ii].Region_mask[j].mask = (mytestregion_a[ii-1].Region_mask[j].mask<<1)+1 ;
            }
        }
    }
    for(int ii=0 ; ii <10 ; ii++)
    {
        mytestpartition.region_1.push_back(& mytestregion_a[ii]);
    }
    mytestpartition.region_2.push_back(NULL);
    mytestpartition.region.push_back(mytestpartition.region_1);
    mytestpartition.region.push_back(mytestpartition.region_2);


}

