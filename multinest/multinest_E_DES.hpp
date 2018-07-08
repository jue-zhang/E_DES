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
#include <gsl/gsl_multimin.h>

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
    static void SetInputParams(const std::vector<double> &input_params);
    static void SetFittedParams(const std::vector<double> &fitted_params);
    static void SetInitConditions(const std::vector<double> &init_conditions);
    // Interfaces for specifying the check pts (i.e., the time instants that one wants to inspect)
    void SetCheckPts(double tI, double tF, int steps); // tI -- initial time, tF -- final time, steps -- num of pts to be inspected
    void SetCheckPts(const std::vector<double> &check_pts_input); // directly specify the to-be-inspected time instants
    // Obtain model parameters
    std::vector<double> GetInputParams() ;
    std::vector<double> GetFittedParams() ;
    std::vector<double> GetConstParams() ;
    std::vector<double> GetCheckPts() {return check_pts;};
    static std::vector<double> GetInitConditions() ;
    std::vector<double> GetCurrentEvolvedParams() ;
    std::vector<double> GetMinFittedParamsSSR() { return minParams_fval_curr;}
    
    // Set the subject type and the corresponding fitted parameters
    //      Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    void SetSubjectTypeFittedParams(const int &type);
    
    // clear the results from the previous run
    static void ClearPreRuns();
    
    // solve the ODEs
    static int Solver_gsl(); // store the results in 'glucoses' and 'insulins' internally
    int Solver_gsl(std::ofstream &output); //store the results in 'glucoses' and 'insulins' internally, and output them in terminal and the specified output file
    
    // Obtain the glucose two hours after food intake
    //      input:
    //          foodIntake -- the amount of intaked food, in unit of 'mg'
    //          bodyMass -- the body mass, in unit of 'kg'
    //          Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //          Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose after two hours of food intake, in unit of 'mmol/L'
    double TwoHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    
    // Obtain the glucose levels in four hours after food intake
    //      input:
    //          foodIntake -- the amount of intaked food, in unit of 'mg'
    //          bodyMass -- the body mass, in unit of 'kg'
    //          Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //          Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose levels in four hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> FourHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    // Obtain the glucose levels in eight hours after food intake
    //      input:
    //          foodIntake -- the amount of intaked food, in unit of 'mg'
    //          bodyMass -- the body mass, in unit of 'kg'
    //          Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //          Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      output: the glucose levels in eight hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> EightHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    // Obtain the glucose levels in 10-min intervals under a set of consective food intake events
    //      input:
    //          bodyMass -- the body mass, in unit of 'kg'
    //          Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //          Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //          foodIntakeEvents -- a list of food intake events in the form of <time, foodIntake>, where "time" is the time
    //                              instant starting from the initial time, and 'foodIntake' is the amount of intaked food,
    //                              in unit of 'mg'. The first element of 'foodIntakeEvents' is the initial time, i.e.,
    //                              foodIntakeEvents[0] = {0., initial_food_intake}, while the last element is the final time
    //                              instant. Note: If we want results at 8h after the final food intake, we should set the last
    //                              element in 'foodIntakeEvents' as {tF + 8*60, 0.}, where 'tF' is the time instant of the final
    //                              non-zero food intake event.
    //      output: the glucose levels at the time instants that are integers of 10's, in the form of <time, glucose>
    std::vector<std::pair<double, double>> GlucoseUnderFoodIntakeEvents(double bodyMass, double Gpl_init_input, double Ipl_init_input, std::vector<std::pair<double, double>> foodIntakeEvents);
    
    
    //  Store the data used for parameter estimation
    //      1. separate glucose and insulin files:
    //          file format: for every row "i" in the glucose file, we have "t[i] glucose[i] glu_err[i]",
    //                       where "glu_err[i]" is the error of "glucose[i]"; same for the insulin files
    void SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_files,
                                       const std::vector<std::string> &dpe_insulin_files);
    //      2. a single file with both glucose and insulin information:
    //          file format: every row "i", "t[i] glucose[i] glu_err[i] insulin[i] ins_err[i]"
    void SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_insulin_files,
                                       const std::vector<std::vector<double>> &input_parameter_sets);
    
    
    static void EstimateFittedParameters_gsl();
    
private:
    // internal usage for using GSL to solve ODEs
    static int gsl_ODEs (double t, const double y[], double f[], void *paramsP);
    
    // internal usage for minimizing the fitted parameters
    static double gsl_min_fitted_params_SSR_func (const gsl_vector *v, void *params);
    static double ComputeSSR(const std::vector<std::vector<double>> &param_est_data_set,
                             const std::vector<double> &glucoses,
                             const std::vector<double> &insulins);
    
    
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
    static double t_init, MGgut_init, Gpl_init, Ipl_init, Jpl_init, Iif_init;
    // evolved parameters at the current time instant (updated after each evolution)
    static double t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Jpl_curr, Iif_curr;
    // time_offset: the acutal time instant when one starts performing evolution (we need this b/c in the actual evolution the starting time instant must be 0!)
    static double time_offset;
    
    // time instants that are to be checked (Note: the first check_pt is always the initial time instant, mostly, t = 0.)
    static std::vector<double> check_pts;
    
    // ouput the obtained glucoses and insulins at the specified time instants in 'time_instants'
    static std::vector<double> glucoses;
    static std::vector<double> insulins;
    static std::vector<double> time_instants;
    
    // store the data-sets for parameter estitmation
    using data_set = std::vector<std::vector<double>>; // data_set[i] = {t[i], glucose[i], glu_err[i], insulin[i], ins_err[i]}
    static std::vector<data_set> param_est_data_sets;
    static std::vector<std::vector<double>> input_param_sets; // store the input parameters, <foodIntake, bodyMass>, for each data set
    static std::vector<double> minParams_init; // the initial values of the parameters to be minimized
    static std::vector<double> minParams_fval_curr; // include the SSR_curr

};



#endif /* E_DES_hpp */
