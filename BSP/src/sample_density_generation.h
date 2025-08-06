


#ifndef SAMPLE_DENSITY_GENERATION_H
#define SAMPLE_DENSITY_GENERATION_H

#include "sampling.h"
#include "beta.h"


static const double L_PI=3.141596;

inline void sample_quarter_circle_2(int n, vector<vector<double> >& data, int d=2) {
	int dim = d;
	data.resize(n);
	int i = 0;
	while (i < n) {
		vector<double> temp;
		temp.resize(dim);
		for (int j = 0; j < dim; j++) {
			temp[j] = rand_double();
		}
		if (temp[d-1] * temp[d-1] + temp[d-2] * temp[d-2] > 1) {
			data[i] = temp;
			i++;
		} else {
			if (rand_double() < .2) {
				data[i] = temp;
				i++;
			}
		}
	}
}

inline double sample_quarter_circle_2_density(const vector<double> &x) {
    double area = L_PI / 4;
    double density = 1/(1-0.8*area);
    if (x[1] * x[1] + x[0] * x[0] <= 1) density /= 5;
    return density;
}


inline void sample_circle_2(int n, vector<vector<double> >& data, int d=2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = rand_double();
        }
        if ((temp[1] - 0.5) * (temp[1] - 0.5) + (temp[0] - 0.5) * (temp[0] - 0.5) < 0.3 * 0.3) {
            data[i] = temp;
            i++;
        }
    }
}


inline double sample_circle_2_density(const vector<double> &x) {
    double area = L_PI * 0.3 * 0.3;
    double density = 1e-10;
    if ((x[1] - 0.5) * (x[1] - 0.5) + (x[0] - 0.5) * (x[0] - 0.5) < 0.3 * 0.3) density = 1/area;
    return density;
}

inline void sample_normal(int n, vector<vector<double> >& data, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = -1;
            while(temp[j] < 0 || temp[j] > 1) {
                temp[j] = rnorm(0.5, 0.1);
            }
        }
        data[i] = temp;
        i++;
    }
}

inline double sample_normal_density(const vector<double> &x) {
    double density = 1.0;
    for (int i = 0; i < (int)x.size(); i++) {
        density *= pnorm(x[i], 0.5, 0.1);
    }
    return density;
}



inline void sample_5normal(int n, vector<vector<double> >& data, int d=1) {
    double scale = 1;
    vector<double> sd(7,0);
    sd[0]= 0.1;
    sd[1]= 0.06;
    sd[2]= 0.044;
    sd[3]= 0.03;
    sd[4]= 0.02;
    sd[5]= 0.013;
    sd[6]= 0.01;
    double shift = 0;
    vector<double> mu(7,0);
    mu[0]= 0.5/scale+ shift;
    mu[1]= 0.4/scale+ shift;
    mu[2]= 0.33/scale+ shift;
    mu[3]= 0.28/scale+ shift;
    mu[4]= 0.26/ scale + shift;
    mu[5]= 0.24/scale+ shift;
    mu[6]= 0.21/scale+ shift;

    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp(1, -1);
        double r = rand_double();
        int choose = (int) floor(r * 7);
        while (temp[0] < 0 || temp[0] > 1) {
            temp[0] = rnorm(mu[choose], sd[choose]);
        }
        data[i] = temp;
        i++;
    }
}


inline double sample_5normal_density(const vector<double> &x) {
     double scale = 1;
    vector<double> sd(7,0);
    sd[0]= 0.1;
    sd[1]= 0.06;
    sd[2]= 0.044;
    sd[3]= 0.03;
    sd[4]= 0.02;
    sd[5]= 0.013;
    sd[6]= 0.01;
    double shift = 0;
    vector<double> mu(7,0);
    mu[0]= 0.5/scale+ shift;
    mu[1]= 0.4/scale+ shift;
    mu[2]= 0.33/scale+ shift;
    mu[3]= 0.28/scale+ shift;
    mu[4]= 0.26/ scale + shift;
    mu[5]= 0.24/scale+ shift;
    mu[6]= 0.21/scale+ shift;
    double density = 0;
    for(int k=0; k<7; k++)    density += pnorm(x[0], mu[k], sd[k])/7;
    return density;
}


inline void sample_mix2normal(int n, vector<vector<double> >& data, mixnormalparas & mixp, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < 2; j++) {
            double u = rand_double();
            if(u<mixp.ratio){
                temp[0] = -1;
                while(temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu1, mixp.sd1);
                }

                temp[1] = -1;
                while(temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm( mixp.mu2, mixp.sd2);
                }
            }
            else{
                temp[0] = -1;
                while(temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu3, mixp.sd1);
                }
                temp[1] = -1;
                while(temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu4, mixp.sd2);
                }

            }
        }
        for(int j=2; j<dim; j++)  temp[j]= rand_double();
        data[i] = temp;
        i++;
    }
}

inline double sample_mix2normal_density(const vector<double> &x, mixnormalparas & mixp) {
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    return density;
}

//using in 'mixnormal'
inline void sample_mix2normalpro(int n, vector<vector<double> >& data, mixnormalparas & mixp, int d = 2) {
    int dim = d;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < 2; j++) {
            double u = rand_double();
            if (u < mixp.ratio) {
                temp[0] = -1;
                while (temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu1, mixp.sd1);
                }

                temp[1] = -1;
                while (temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu2, mixp.sd2);
                }
            } else {
                temp[0] = -1;
                while (temp[0] < 0 || temp[0] > 1) {
                    temp[0] = rnorm(mixp.mu3, mixp.sd1);
                }
                temp[1] = -1;
                while (temp[1] < 0 || temp[1] > 1) {
                    temp[1] = rnorm(mixp.mu4, mixp.sd2);
                }

            }
        }
        for (int j = 2; j < 4 && j < dim; j++) {
            temp[j] = -1;
            while (temp[j] < 0 || temp[j] > 1) {
                temp[j] = rnorm(0.5, 0.1);
            }
        }
        for (int j = 4; j < dim; j++) {
            if (rand_double() > 0.5) {
                temp[j] = -1;
                while (temp[j] < 0 || temp[j] > 1) {
                    temp[j] = rnorm(0.35, 0.1);
                }
            } else {
                temp[j] = -1;
                while (temp[j] < 0 || temp[j] > 1) {
                    temp[j] = rnorm(0.6, 0.05);
                }
            }
        }


        // the >30 dimensional case, replace the last two dimensions temp[dim-2] and temp[dim-1]
        if (dim > 31) {
            double ruf = rand_double();
            if (ruf < 0.2) {
                double a=0;
                double b=-1;
                while(a>=b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]=a;
                temp[dim-1]=b;
            }
            else if (ruf < 0.7) {
                double a=rand_double();
                double b=rand_double();
                if(a+b>=1){
                    temp[dim-2] = 0.5+0.5*a;
                    temp[dim-1] = 0.5*b;
                }
                else{
                    temp[dim-2]= 1-0.5*b;
                    temp[dim-1]= 0.5+0.5*a;
                }
            }
            else if (ruf < 0.8) {
                double a=-2;
                double b=-1;
                while(a<b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]= 1-0.5*a;
                temp[dim-1]= 0.5*b;

            }
            else {
                double a=-2;
                double b=-1;
                while(a<b){
                    a=rand_double();
                    b=rand_double();
                }
                temp[dim-2]= 0.5*a;
                temp[dim-1]= 0.5*b;
            }

        }
        data[i] = temp;
        i++;

    }
    return;

}

//using
inline double sample_mix2normalpro_density(const vector<double> &x, mixnormalparas & mixp) {
    int dim = (int)x.size();
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    for(int j=2; j<4 && j<dim; j++)  density *= pnorm(x[j], 0.5, 0.1);
    if(dim<=31){
        for(int j=4; j<dim; j++)  density *= (0.5*pnorm(x[j],0.35,0.1)+0.5*pnorm(x[j], 0.6,0.05));
        return density;
    }
    for(int j=4; j<dim-2; j++)  density *= (0.5*pnorm(x[j],0.35,0.1)+0.5*pnorm(x[j], 0.6,0.05));
    if(x[dim-2]<x[dim-1]) density *= 0.4;
    else if((x[dim-2]+x[dim-1])>1) density *= 2;
    else if(x[dim-2]>0.5)  density *= 0.8;
    else density *= 1.6;
    return density;

}


inline double sample_mix2normalpro_sep_density(const vector<double> &x, mixnormalparas & mixp){
    double density = mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)*pnorm(x[1],mixp.mu2,mixp.sd2)+(1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1)*pnorm(x[1], mixp.mu4,mixp.sd2);
    density = density/ (mixp.ratio*pnorm(x[0],mixp.mu1,mixp.sd1)+ (1-mixp.ratio)*pnorm(x[0], mixp.mu3,mixp.sd1));
    density = density/ (mixp.ratio*pnorm(x[1],mixp.mu2,mixp.sd2)+ (1-mixp.ratio)*pnorm(x[1], mixp.mu4,mixp.sd2));
    return density;
}


//Confidence interval: [-5,5]            temp = (temp +5)/10;
inline void sample_skewedmixnormal(int n, vector<vector<double> >& data) {
    data.resize(n);
    for (int j = 0; j < n; j++) {
        int i = (int)floor(rand_double()*8);
        if (i == 8) i = 7;
        double temp = rnorm(3 * (pow((double)2.0 / 3.0, (int)i) - 1), pow((double)2.0 / 3.0, (int)i));
        temp = (temp +5.0)/10.0;
        data[j].push_back(temp);
        if(temp<0 && temp>1) cout<<"out of range!\n";
    }
}

inline double sample_skewedmixnormal_density(const vector<double> &x) {
    double t = x[0]*10.0-5.0;
    double density = 0;
    for(int i=0; i<8; i++)  {
        density +=   10.0/8.0 * pnorm(t, 3 * (pow(2.0 / 3.0, i) - 1),pow(2.0 / 3.0, i));         // jacobian: 10
    }

    return density;

}

inline void sample_mixstepnormal(int n, vector<vector<double> >& data, int d=6) {

    double temp;
    data.resize(n);
    for (int i = 0; i < n/2; i++) {
        data[i].push_back(rand_double() * .25 + .25);
        data[i].push_back(rand_double() * .125 + .75);
        data[i].push_back(rand_double() * .5);
        data[i].push_back(rand_double() * .5 + .5);
        temp = -1;
        while(temp < 0 || temp > 1) {
            temp = rnorm(0.3, 0.05);
        }
        data[i].push_back(temp);

        temp = -1;
        while(temp < 0 || temp > 1) {
            temp = rnorm(0.6, 0.1);
        }
        data[i].push_back(temp);
    }
    for (int i = n/2; i < n; i++) {
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
    }
}

inline double sample_mixstepnormal_density(const vector<double> &x) {
    double density1 = 0;

    if (x[0] >= 0.25 && x[0] <= 0.5 && x[1] >= 0.75 && x[1] <= 0.875 && x[2] <= 0.5 && x[3] >= 0.5) density1 += 1 / 0.25 / 0.125 / 0.5 / 0.5;
    double density2 = pnorm(x[4],0.3,0.05)*pnorm(x[5],0.6,0.1);
    return (0.5 + 0.5* density1*density2);
}

inline void sample_beta15(int n,vector<vector<double> >& data) {
    data.resize(n);
    for (int i = 0; i < n ; i++) {
        double ind = rand_double();
 //       if(ind>0.7)     data[i].push_back(rbeta(rand_double(),300,100));
 //       else  data[i].push_back(rbeta(rand_double(), 3,12));
        if(ind>0.7)     data[i].push_back(rbeta(rand_double(),300,10));
        else  data[i].push_back(rbeta(rand_double(), 3,12));
    }
}

inline double sample_beta15_density(const vector<double> &x) {
    double density = 0.3*beta(x[0],300,10) + 0.7*beta(x[0], 3,12);
    return density;
}



inline void sample_1(int n,vector<vector<double> >& data) {

    data.resize(n);
    for (int i = 0; i < n / 2; i++) {
        data[i].push_back(rand_double() * .125 + 0.5);
    }
    for (int i = n / 2; i < n; i++) {
        data[i].push_back(rand_double());
    }
}

inline double sample_1_density(const vector<double> &x) {
    double density = 0.5;
    if (x[0] >= 0.5 && x[0] <= 0.625) density += 0.5 / 0.125;
    return density;
}


inline void sample_5(int n, vector<vector<double> >& data) {

    data.resize(n);
    for (int i = 0; i < n / 2; i++) {
        data[i].push_back(rand_double() * .25 + .25);
        data[i].push_back(rand_double() * .125 + .75);
        data[i].push_back(rand_double() * .5);
        data[i].push_back(rand_double() * .5 + .5);
        data[i].push_back(rand_double());
    }
    for (int i = n / 2; i < n; i++) {
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
        data[i].push_back(rand_double());
    }
}

inline double sample_5_density(const vector<double> &x) {
    double density = 0.5;
    if (x[0] >= 0.25 && x[0] <= 0.5 && x[1] >= 0.75 && x[1] <= 0.875 && x[2] <= 0.5 && x[3] >= 0.5) density += 0.5 / 0.25 / 0.125 / 0.5 / 0.5;
    return density;
}

inline void sample_10(int n,  vector<vector<double> >& data) {
    int dim = 10;
    data.resize(n);
    int i = 0;
    while (i < n) {
        vector<double> temp;
        temp.resize(dim);
        for (int j = 0; j < dim; j++) {
            temp[j] = rand_double();
        }
        if (temp[1] > 0.8 * sin(temp[0] * M_PI * 0.8)) {
            data[i] = temp;
            i++;
        } else {
            if (rand_double() < .2) {
                data[i] = temp;
                i++;
            }
        }
    }
}

inline double sample_10_density(const vector<double> &x) {
    double area = (1-cos(0.8*M_PI))/M_PI;
    double density = 1/(1-0.8*area);
    if (x[1] <= 0.8 * sin(x[0] * M_PI * 0.8)) density /= 5;
    return density;
}



#endif //





