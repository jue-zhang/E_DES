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


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

class E_DES{
    
    // ++++++++++++++++++++++++++++++   High level methods  ++++++++++++++++++++++++++++++
public:
    // Method:
    //      Obtain the glucose two hours after food intake
    // Input:
    //      foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    // Output:
    //      the glucose after two hours of food intake, in unit of 'mmol/L'
    double TwoHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    // Method:
    //      Obtain the glucose levels in four hours after food intake
    // Input:
    //      foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    // Output:
    //      the glucose levels in four hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> FourHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    
    // Method:
    //      Obtain the glucose levels in eight hours after food intake
    // Input:
    //      foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    // Output:
    //      the glucose levels in eight hours (10 min intervals), in unit of 'mmol/L'
    std::vector<double> EightHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    
    // Method:
    //      Obtain the glucose levels in eight hours after food intake (EVERY MINIMUTE)
    // Input:
    //      foodIntake -- the amount of intaked food, in unit of 'mg'
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    // Output:
    //      the glucose levels in eight hours (per min), in unit of 'mmol/L'
    std::vector<double> EightHourGlucosePerMin(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input);
    
    // Method:
    //      Obtain the glucose levels in 10-min intervals under a set of consective food intake events
    // Input:
    //      bodyMass -- the body mass, in unit of 'kg'
    //      Gpl_init_input -- the initial glucose in the plasma before eating, in unit of 'mmol/L'; ref. value = 5 mmol/L
    //      Ipl_init_input -- the initial insulin in the plasma before eating, in unit of 'mU/L'; ref. value = 8 mU/L
    //      foodIntakeEvents -- a list of food intake events in the form of <time, foodIntake>, where "time" is the time
    //              instant starting from the initial time, and 'foodIntake' is the amount of intaked food, in unit of 'mg'.
    //              The first element of 'foodIntakeEvents' is the initial time, i.e.,
    //              foodIntakeEvents[0] = {0., initial_food_intake},  while the last element is the final time instant.
    //              Note: If we want results at 8h after the final food intake, we should set the last element in
    //              'foodIntakeEvents' as {tF + 8*60, 0.}, where 'tF' is the time instant of the final non-zero food intake event.
    // Output:
    //      the glucose levels at the time instants that are integers of 10's, in the form of <time, glucose>
    std::vector<std::pair<double, double>> GlucoseUnderFoodIntakeEvents(double bodyMass, double Gpl_init_input, double Ipl_init_input, std::vector<std::pair<double, double>> foodIntakeEvents);
    

    
    // ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++
public:
    // -------------  Modifying model paramers
    
    static void SetInputParams(const std::vector<double> &input_params);
    static void SetFittedParams(const std::vector<double> &fitted_params);
    static void SetInitConditions(const std::vector<double> &init_conditions);
    // Specifying the check pts (i.e., the time instants that one wants to inspect)
    void SetCheckPts(double tI, double tF, int steps); // tI -- initial time, tF -- final time, steps -- num of pts to be inspected
    void SetCheckPts(const std::vector<double> &check_pts_input); // directly specify the to-be-inspected time instants
    // Set the subject type and the corresponding fitted parameters
    //      Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    void SetSubjectTypeFittedParams(const int &type); // use pre-seted fitted-params
    // Load the optimized fitted-params
    void LoadFittedParams(std::ifstream &param_file);
    
    // -------------  Extracting model paramers
    
    std::vector<double> GetInputParams() ;
    std::vector<double> GetFittedParams() ;
    std::vector<double> GetConstParams() ;
    std::vector<double> GetCheckPts() {return check_pts;};
    static std::vector<double> GetInitConditions() ;
    std::vector<double> GetCurrentEvolvedParams() ;
    std::vector<double> GetMinFittedParamsSSR() { return minParams_fval_curr;} // get minimized fitted-params + SSR_min
    std::vector<std::pair<double, double>> GetGlucoses(); // return the obtained glucoses at the specied time instants
    std::vector<std::pair<double, double>> GetInsulins(); // return the obtained insulins at the specied time instants
    
    // ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++
public:
    // clear the results from the previous run
    static void ClearPreRuns();
    // solve the ODEs
    static int Solver_gsl(); // store the results in 'glucoses' and 'insulins' internally
private:
    // internal usage for using GSL to solve ODEs
    static int gsl_ODEs (double t, const double y[], double f[], void *paramsP);
    
    
    // ++++++++++++++++++++++++++++++   Methods for estimating fitted-params  ++++++++++++++++++++++++++++++
public:
    
    // Method:
    //      Set up the data sets used for estimating fitted-params (separate glucose and insulin files)
    // Input:
    //      dpe_glucose_files -- paths of glucose files
    //      dpe_insulin_files -- paths of insulin files
    //      file format: for every row "i" in the glucose file, we have "t[i] glucose[i] glu_err[i]",
    //                       where "glu_err[i]" is the error of "glucose[i]"; same for the insulin files
    // Output:
    //      'param_est_data_sets' stores all the data sets used for estimating fitted-params
    void SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_files,
                                       const std::vector<std::string> &dpe_insulin_files);
    // Method:
    //      Set up the data sets used for estimating fitted-params (a single file with both glucose and insulin information)
    // Input:
    //      dpe_glucose_insulin_files -- paths of files
    //      file format: every row "i", "t[i] glucose[i] glu_err[i] insulin[i] ins_err[i]"
    // Output:
    //      'param_est_data_sets' stores all the data sets used for estimating fitted-params
    void SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_insulin_files,
                                       const std::vector<std::vector<double>> &input_parameter_sets);
    
    // Estimate fitted-params
    static void EstimateFittedParameters_gsl();
    
private:
    // internal usage for minimizing the fitted parameters
    static double gsl_min_fitted_params_SSR_func (const gsl_vector *v, void *params);
    static double ComputeSSR(const std::vector<std::vector<double>> &param_est_data_set,
                             const std::vector<double> &glucoses,
                             const std::vector<double> &insulins);
    
    
    // ++++++++++++++++++++++++++++++   Methods for class itself  ++++++++++++++++++++++++++++++
public:
    E_DES() = default;
    
    // ++++++++++++++++++++++++++++++   Model Parameters  ++++++++++++++++++++++++++++++
private:
    
    // ---------------  Subject type and the corresponding fitted parameters
    
    // type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    static int type;
    // pre_setted fitted params
    static double k1_H, k2_H, k3_H, k4_H, k5_H, k6_H, k7_H, k8_H, k9_H, k10_H, k11_H, k12_H, sigma_H, KM_H; // healthy
    static double k1_D1, k2_D1, k3_D1, k4_D1, k5_D1, k6_D1, k7_D1, k8_D1, k9_D1, k10_D1, k11_D1, k12_D1, sigma_D1, KM_D1; // D1
    static double k1_D2, k2_D2, k3_D2, k4_D2, k5_D2, k6_D2, k7_D2, k8_D2, k9_D2, k10_D2, k11_D2, k12_D2, sigma_D2, KM_D2; // D2
    
    // ---------------- Parameters necessary for evolution
    
    // input params: Dmeal -- amount of food intake (mg); Mb -- body mass (Kg)
    static double Dmeal, Mb;
    // fitted params: params that are to be obtained through fitting
    static double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, sigma, KM;
    // constant params in the model
    static double gbliv, Gthpl, vG, vI, beta, fC, tau_i, t_int, tau_d, c1;
    // initial_conditions: status before eating anything (Note: Jpl is the first order derivative of Ipl )
    static double t_init, MGgut_init, Gpl_init, Ipl_init, Jpl_init, Iif_init;
    // time_offset: the acutal time instant when one starts performing evolution. We need this b/c
    //              in the actual evolution the starting time instant must be 0!
    static double time_offset;
    
    // ---------------- Vars used for checking evolution
    
    // time instants that are to be checked (Note: the first check_pt is always the initial time instant, t = 0.)
    static std::vector<double> check_pts;
    // evolved parameters at the current time instant (updated after each evolution)
    static double t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Jpl_curr, Iif_curr;
    // output the obtained glucoses and insulins at the specified time instants in 'time_instants'
    static std::vector<double> time_instants;
    static std::vector<double> glucoses;
    static std::vector<double> insulins;
    
    // ---------------- Vars used for estimating fitted-params
    
    // store the input parameters, <foodIntake, bodyMass>, for each data set
    static std::vector<std::vector<double>> input_param_sets;
    // store the correponding set of data points corresponding to each set of input parameters.
    //      format of each data set: data_set[i] = {t[i], glucose[i], glu_err[i], insulin[i], ins_err[i]}
    using data_set = std::vector<std::vector<double>>;
    static std::vector<data_set> param_est_data_sets;
    // initial values of the parameters to be minimized, NOT include the SSR_init
    static std::vector<double> minParams_init;
    // current values of the parameters to be minimized, include the SSR_curr
    static std::vector<double> minParams_fval_curr;

};



#endif /* E_DES_hpp */
