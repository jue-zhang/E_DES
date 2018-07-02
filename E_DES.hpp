//
//  E_DES.hpp
//  E_DES
//
//  Created by Jue on 6/24/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#ifndef E_DES_hpp
#define E_DES_hpp

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//#include "multinest.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

class E_DES{
    friend void MultiNest_LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
public:
    E_DES() = default;
    // Interfaces for modifing model paramers
    void SetInputParams(const std::vector<double> &input_params);
    void SetFittedParams(const std::vector<double> &fitted_params);
    void SetInitConditions(const std::vector<double> &init_conditions);
    // Interfaces for specifying the check pts (i.e., the time instants that one wants to inspect)
    void SetCheckPts(double tI, double tF, int steps); // tI -- initial time, tF -- final time, steps -- num of pts to be inspected
    void SetCheckPts(const std::vector<double> &check_pts_input); // directly specify the to-be-inspected time instants
    // Obtain model parameters
    std::vector<double> GetInputParams() ;
    std::vector<double> GetFittedParams() ;
    std::vector<double> GetConstParams() ;
    std::vector<double> GetCheckPts() {return check_pts;};
    std::vector<double> GetInitConditions() ;
    
    // Set the subject type and the corresponding fitted parameters
    //      Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    void SetSubjectTypeFittedParams(const int &type);
    
    // clear the results from the previous run
    void ClearPreRuns();
    
    // solve the ODEs
    int Solver_gsl(); // store the results in 'glucoses' and 'insulins' internally
    int Solver_gsl(std::ofstream &output); //store the results in 'glucoses' and 'insulins' internally, and output them in terminal and the specified output file
    
    // Obtain the glucose two hours after food intake
    //      input: foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose after two hours of food intake, in unit of 'mmol/L'
    double TwoHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    
    // Obtain the glucose levels in four hours after food intake
    //      input: foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose levels in four hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> FourHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    // Obtain the glucose levels in eight hours after food intake
    //      input: foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose levels in eight hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> EightHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    
    
    void SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_files,
                                       const std::vector<std::string> &dpe_insulin_files);
//    void EstimateParameters();
    
private:
    // internal usage for using GSL to solve ODEs
    static int gsl_ODEs (double t, const double y[], double f[], void *paramsP);
//    // for MultiNest
//    static void MultiNest_LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
//    static void MultiNest_Dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context) {}
    
private:
    // type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    static int type;
    // input params: Dmeal -- amount of food intake (mg); Mb -- body mass (Kg)
    static double Dmeal, Mb;
    // fitted params: params that are to be obtained through fitting
    static double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, sigma, KM;
    // constant params in the model
    static double gbliv, Gthpl, vG, vI, beta, fC, tau_i, t_int, tau_d, c1;
    // initial_conditions: status before eating anything (Note: Jpl is the first order derivative of Ipl )
    static double MGgut_init, Gpl_init, Ipl_init, Jpl_init, Iif_init;
    // time instants that are to be checked (Note: the first check_pt is always the initial time instant, mostly, t = 0.)
    std::vector<double> check_pts;
    // store the obtained glucoses and insulins at the specified time instants in check_pts
    std::vector<double> glucoses;
    std::vector<double> insulins;
    
    // store the data-sets for parameter estitmation
    using data_set = std::vector<std::vector<double>>; // data_set = time, glucose, insulin
    std::vector<data_set> param_est_data_sets;

};



#endif /* E_DES_hpp */
