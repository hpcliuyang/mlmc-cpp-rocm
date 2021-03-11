#include <iostream>
#include <math.h>  
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;

double eval_P0(int nsamples, int nsteps, double paras[])
{
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0.0,1.0);

    double T = paras[0];
    double r = paras[1];
    double sig = paras[2];
    double S0 = paras[3];
    double K = paras[4];

    double hf = T/nsteps;
    double Sf[nsamples];
    double Pf[nsamples];
    double dWf[nsamples];

    for (int i = 0; i < nsamples; i ++)
    {
        Sf[i] = S0;
    }
    for (int j = 0; j < nsteps; j ++){
    for (int i = 0; i < nsamples; i ++){
        dWf[i] = distribution(generator);
        dWf[i] = sqrt(hf)*dWf[i];
        Sf[i] = Sf[i] + r*hf*Sf[i] + sig*Sf[i]*dWf[i];
    }
    }

    double Psum = 0.0;
    for (int i = 0; i < nsamples; i ++)
    {
        Pf[i] = max(0.0,Sf[i]-K);
        Pf[i] = exp(-r*T)*Pf[i];
        Psum = Psum+Pf[i];
    }
    double Pavg = Psum/nsamples;
    return Pavg;
}



double eval_correction(int nsamples, int nsteps, double paras[], int Mfactor)
{
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(0.0,1.0);

    double T = paras[0];
    double r = paras[1];
    double sig = paras[2];
    double S0 = paras[3];
    double K = paras[4];

    double hc = T/nsteps;
    double hf = hc/Mfactor;
    double Sf[nsamples];
    double Sc[nsamples];
    double Pf[nsamples];
    double Pc[nsamples];
    double dWf[nsamples];
    double dWc[nsamples];

    for (int i = 0; i < nsamples; i ++){
        Sf[i] = S0;
        Sc[i] = S0;
    }
    for (int j = 0; j < nsteps; j ++){
    for (int i = 0; i < nsamples; i ++){
        dWc[i] = 0.0;
    }
    for (int m = 0; m < Mfactor; m ++){
        for (int i = 0; i < nsamples; i ++){
        dWf[i] = distribution(generator);
        dWf[i] = sqrt(hf)*dWf[i];
        Sf[i] = Sf[i] + r*hf*Sf[i] + sig*Sf[i]*dWf[i];
        dWc[i] = dWc[i]+dWf[i];
        }
    }
        for (int i = 0; i < nsamples; i ++){
        Sc[i] = Sc[i] + r*hc*Sc[i] + sig*Sc[i]*dWc[i];
        }    
    }

    double Pfsum = 0.0;
    double Pcsum = 0.0;

    for (int i = 0; i < nsamples; i ++){
        Pf[i] = max(0.0,Sf[i]-K);
        Pf[i] = exp(-r*T)*Pf[i];
        Pfsum = Pfsum+Pf[i];
        Pc[i] = max(0.0,Sc[i]-K);
        Pc[i] = exp(-r*T)*Pc[i];
        Pcsum = Pcsum+Pc[i];
    }

    cout  << "Pf  " << Pfsum/nsamples << endl;
    cout  << "Pc  " << Pcsum/nsamples << endl;        

    double Pcor = (Pfsum-Pcsum)/nsamples;

    return Pcor;
}




int main()
{
    int Msamples[3];
    Msamples[0] = 100000;
    Msamples[1] = 25000;
    Msamples[2] = 5000;

    int Mfactor = 4;
    int Msteps[3];
    Msteps[0] = 1;
    Msteps[1] = 1;
    Msteps[2] = Msteps[1]*Mfactor;

    double Pc[2]; 
    
    double paras[5];
//    double T = 1.0;
//    double r = 0.05;
//    double sig = 0.2;
//    double S0 = 100.0;
//    double K = 100.0;
    paras[0] = 1.0; 
    paras[1] = 0.05; 
    paras[2] = 0.2;
    paras[3] = 100.0;
    paras[4] = 100.0;

    int Msmp = Msamples[0];
    int Mstp = Msteps[0];
    double P0  = eval_P0(Msmp,Mstp,paras);

    Msmp = Msamples[1];
    Mstp = Msteps[1];    
    Pc[0]  = eval_correction(Msmp,Mstp,paras,Mfactor);

    Msmp = Msamples[2];
    Mstp = Msteps[2];
    Pc[1]  = eval_correction(Msmp,Mstp,paras,Mfactor);

    double P = P0+Pc[0]+Pc[1];

    cout  << "Evaluate P0  " << P0 << endl;
    cout  << "Evaluate correction Pc  " << Pc[0] << endl;
    cout  << "Evaluate correction Pc  " << Pc[1] << endl;
    cout  << "Final P  " << P << endl;

    return 0;
}

