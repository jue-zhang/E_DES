//
//  main.cpp
//  E_DES
//
//  Created by Jue on 6/24/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "multinest.h"

#include "E_DES.hpp"

using namespace std;

//void MultiNest_LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context){
//    
//    E_DES* modP = static_cast<E_DES *>(context);
//    
////    // for k1:k8, k11, k12, sigma, scatter pts within [10^-5, 10^2]
////    modP->k1 = pow(10, -5 + Cube[0] * 7);
////    modP->k2 = pow(10, -5 + Cube[1] * 7);
////    modP->k3 = pow(10, -5 + Cube[2] * 7);
////    modP->k4 = pow(10, -5 + Cube[3] * 7);
////    modP->k5 = pow(10, -5 + Cube[4] * 7);
////    modP->k6 = pow(10, -5 + Cube[5] * 7);
////    modP->k7 = pow(10, -5 + Cube[6] * 7);
////    modP->k8 = pow(10, -5 + Cube[7] * 7);
////    //k9
////    //k10
////    modP->k11 = pow(10, -5 + Cube[8] * 7);
////    modP->k12 = pow(10, -5 + Cube[9] * 7);
////    modP->sigma = pow(10, -5 + Cube[10] * 7);
////    // for KM, scatter pts within [1, 10^3]
////    modP->KM = pow(10, Cube[11]*3);
//    
//    // for k1:k8, k11, k12, sigma, scatter pts within [10^-5, 10^2]
//    modP->k1 = 0.1 * Cube[0];
//    modP->k2 = Cube[1];
//    modP->k3 = 0.01 * Cube[2];
//    modP->k4 = 0.001 * Cube[3];
//    modP->k5 = Cube[4];
//    modP->k6 = Cube[5];
//    modP->k7 = 10 * Cube[6];
//    modP->k8 = 10 * Cube[7];
//    //k9
//    //k10
//    modP->k11 = 0.1 * Cube[8];
//    modP->k12 = Cube[9];
//    modP->sigma = 10 * Cube[10];
//    // for KM, scatter pts within [1, 10^3]
//    modP->KM = 100 * Cube[11];
//    
//    const int fp_sz_1 = 8; // k1 : k8
//    const int fp_sz_2 = 3; // k11, k12, sigma
//
//    
////    std::vector<double> fitted_params;
////    const int fp_sz_1 = 8; // k1 : k8
////    for (int i = 0; i < fp_sz_1; ++i)
////        fitted_params.push_back(pow(10, -5+Cube[i]*7)); // for k1:k8, scatter pts within [10^-5, 10^2]
////    fitted_params.push_back(0); // k9
////    fitted_params.push_back(0); // k10
////    const int fp_sz_2 = 3; // k11, k12, sigma
////    for (int i = 0; i < fp_sz_2; ++i)
////        fitted_params.push_back(pow(10, -5+Cube[fp_sz_1+i]*7)); // for k11, k12, sigma, scatter pts within [10^-5, 10^2]
////    fitted_params.push_back(pow(10, Cube[fp_sz_1+fp_sz_2]*3)); // KM within [1, 10^3]
////    modP->SetFittedParams(fitted_params);
//    
//    double SSR_tot = 0.;
//    for (const auto &param_est_data_set: modP->param_est_data_sets) {
//        // set the input params and check_pts, according to the data used for parameter estimation
//        modP->check_pts = param_est_data_set[0];
//        // evolve the system
//        modP->Solver_gsl();
//        // calculate SSR
//        double SSR_tmp = 0.;
//        vector<double> glucoses = modP->glucoses;
//        vector<double> insulins = modP->insulins;
//        vector<double> glucoses_t = param_est_data_set[1];
//        vector<double> insulins_t = param_est_data_set[2];
//        for (int i = 0; i < glucoses.size(); ++i) {
//            SSR_tmp += pow(glucoses[i] - glucoses_t[i], 2);
//            SSR_tmp += pow((insulins[i] - insulins_t[i])/10., 2);
////            cout << SSR_tmp << endl;
//        }
//        SSR_tot += SSR_tmp;
//    }
//    
//    cout << "------------" << SSR_tot << endl;
//    lnew = - SSR_tot;
//    
//    // store the scattered points in terms of k1:k8, k11, k12, sigma, KM
//    vector<double> fitted_params = modP->GetFittedParams();
//    for (int i = 0; i < fp_sz_1; ++i){
//        Cube[i] = fitted_params[i];
//    }
//    for (int i = 0; i < fp_sz_2; ++i){
//        Cube[fp_sz_1+i] = fitted_params[fp_sz_1+2+i]; // skip k9 and k10 in fitted_params
//    }
//    Cube[fp_sz_1+fp_sz_2] = fitted_params[fp_sz_1+fp_sz_2+2];
//
//
//}
//
//
//void MultiNest_Dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context) {}


int main(int argc, const char * argv[]) {
    
//    string root_dir = "/Users/Jue/Desktop/precision_health/";
//    
//    ofstream check_output;
//    string output_file_path = root_dir + "models/E_DES/E_DES/output.dat";
//    check_output.open(output_file_path, std::ofstream::out);
//    
//    double tI = 0.;
//    double tF = 240.;
//    int steps = 10;
//
//    E_DES eDES;
//    eDES.SetCheckPts(tI, tF, steps);
//    //    eDES.Solver_gsl();
//    eDES.Solver_gsl(check_output);
//    check_output.close();
    
    E_DES eDES;
    double foodIntake = 75E3; // unit: mg
    double bodyMass = 75.; // unit: kg
    double Gpl_init = 5.; // unit: mmol/L
    double Ipl_init = 8.; // unit: mU/L
    cout << eDES.TwoHourGlucose(foodIntake, bodyMass, Gpl_init, Ipl_init) << endl;
    
//    ifstream data_para_est_file;
//    string dpe_file_dir = root_dir + "/data/E_DES/data_parameter_estimation/clean_data/";
//    vector<string> dpe_file_names = {"34_female", "34_male", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44"};
//    string dpe_glucose_file_suffix = "_plasma_glucose_cleaned.dat";
//    string dpe_insulin_file_suffix = "_plasma_insulin_cleaned.dat";
//    vector<string> dpe_glucose_files, dpe_insulin_files;
//    for (auto dpe_file_name: dpe_file_names){
//        dpe_glucose_files.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_file_suffix);
//        dpe_insulin_files.push_back(dpe_file_dir + dpe_file_name + dpe_insulin_file_suffix);
//    }
//    
//    E_DES eDES;
//    eDES.SetDataForParameterEstimation(dpe_glucose_files, dpe_insulin_files);
//
//    
//    // set the MultiNest sampling parameters
//    
//    int IS = 0;					// do Nested Importance Sampling?
//    
//    int mmodal = 1;					// do mode separation?
//    
//    int ceff = 0;					// run in constant efficiency mode?
//    
//    int nlive = 5000;				// number of live points
//    
//    double efr = 0.5;				// set the required efficiency
//    
//    double tol = 0.05;				// tol, defines the stopping criteria
//    
//    int ndims = 12;					// dimensionality (no. of free parameters)
//    
//    int nPar = 12;					// total no. of parameters including free & derived parameters
//    
//    int nClsPar = 12;				// no. of parameters to do mode separation on
//    
//    int updInt = 10000;				// after how many iterations feedback is required & the output files should be updated
//    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
//    
//    double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
//    
//    int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
//    
//    int pWrap[ndims];				// which parameters to have periodic boundary conditions?
//    //for(int i = 0; i < ndims; i++) pWrap[i] = 0;
//    
//    char root[100] = "chains/test_";			// root for output files
//    
//    int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
//    
//    int fb = 1;					// need feedback on standard output?
//    
//    int resume = 0;					// resume from a previous job?
//    
//    int outfile = 1;				// write output files?
//    
//    int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
//    // set it to F if you want your main program to handle MPI initialization
//    
//    double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
//    
//    int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
//    // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
//    
//    void *context = 0;				// not required by MultiNest, any additional information user wants to pass
//    
//    context = &eDES;
//    
//    // calling MultiNest
//    nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, MultiNest_LogLike, MultiNest_Dumper, context);

    
    
    return 0;
}







