//
//  E_DES_Ex.hpp
//  E_DES
//
//  Created by Jue on 7/13/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#ifndef E_DES_Ex_hpp
#define E_DES_Ex_hpp

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multimin.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <map>


class E_DES_Ex{
    friend class SubjectGlucose;
    
    // ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++
public:
    // -------------  Modifying model paramers
    
    // set the vars. specific to the evolution: <type, bodyMass, Dmeal, Gpl_basl, Ipl_basal, ex_intensity, ex_type>
    void SetEvolutionSpecifics(const std::tuple<int, double, double, double, double, double, double, int> &specifics);
    // set the initial conditions
    void SetInitConditions(const std::vector<double> &init_conditions);
    // set fitted-params
    void SetFittedParams(const std::vector<double> &fitted_params);
    void SetFittedParams(const int &type); // use pre-seted fitted-params
//    void SetFittedParams(std::ifstream &param_file); // Load the optimized fitted-params from file
    void SetFittedParams(const std::string &param_file_path); // Load the optimized fitted-params from file
    void SetFittedParams(std::map<std::string, double> dict_fitted_params);
    // Specifying the check pts (i.e., the time instants that one wants to inspect)
    void SetCheckPts(double tI, double tF, int steps); // tI -- initial time, tF -- final time, steps -- num of pts to be inspected
    void SetCheckPts(const std::vector<double> &check_pts_input) { check_pts = check_pts_input;}
    
    // -------------  Extracting model paramers
    std::tuple<int, double, double, double, double, double, double, double, int> GetEvolutionSpecifics() {
        return {type, Mb, Gpl_basal, Ipl_basal, this_meal, latest_meal, delta_eat, ex_type};
    }
    std::vector<double> GetInitConditions() {
        return {t_init, MGgut_init, Gpl_init, Ipl_init, Jpl_init, Iif_init, ExPre_init, Ex_init};
    }
    std::vector<double> GetFittedParams() {
        return {k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, sigma, KM, k4e, k5e, k6e, k8e, lam};
    }
    std::vector<double> GetConstParams() { return {gbliv, Gthpl, vG, vI, beta, fC, tau_i, t_int, tau_d, c1}; }
    std::vector<double> GetCurrentEvolvedParams() {
        return {t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Jpl_curr, Iif_curr, ExPre_curr, Ex_curr};
    }
    
    std::map<std::string, double> GetFittedParamsDict();
    
    // ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++
public:
    // clear the results from the previous run
    void ClearPreRunOutput();
    // solve the ODEs
    int Solver_gsl(); // store the results in 'glucoses' and 'insulins' internally
private:
    // internal usage for using GSL to solve ODEs
    static int gsl_ODEs (double t, const double y[], double f[], void *paramsP);
    
    
    // ++++++++++++++++++++++++++++++   Methods for class itself  ++++++++++++++++++++++++++++++
public:
    E_DES_Ex() = default;
    
    // ++++++++++++++++++++++++++++++   Model Parameters  ++++++++++++++++++++++++++++++
private:
    
    // ---------------- Parameters necessary for evolution
    
    // Evolution specifics: overall setting
    int type = 0; // type of the subject, which determines the set of 'fitted-params'
    double Mb = 75.; // kg
    double Gpl_basal = 5.; // basal glucose in plasma, mmol/L
    double Ipl_basal = 8.; // basal insulin in plasma, U/L
    double this_meal = 75E3; // mg
    double latest_meal = 0.; // the amount of intaked food in the latest food intake event
    double delta_eat = 0.; // how much time has passed since the latest food intake
    int ex_type = 0;
    
    // Initial_conditions
    double t_init = 0.;      // min
    double MGgut_init= 0.;  // mg
    double Gpl_init = 5.;   // mmol/L
    double Ipl_init = 8.;   // U/L
    double Jpl_init = 0.;   // mU/L/min
    double Iif_init = 0.;   // mU/L
    double ExPre_init = 0.; // pre-existing exercise intensity
    double Ex_init = 0.;    // total exercise intensity
    
    // fitted params: params that can obtained through fitting
    double k1 = 0.015262;
    double k2 = 0.304526;
    double k3 = 0.00808413;
    double k4 = 0.000177229;
    double k5 = 0.0798414;
    double k6 = 0.28092;
    double k7 = 0.0321147;
    double k8 = 6.86228;
    double k9 = 0.;         // short-acting insulin
    double k10 = 0.;        // short-acting insulin
    double k11 = 0.0271874;
    double k12 = 0.294954;
    double sigma = 1.57349;
    double KM = 19.3604;
    double k4e = 1.;        // exercise
    double k5e = 1.5;        // exercise
    double k6e = 0.5;        // exercise
    double k8e = 0.5;        // exercise
    double lam = 1/60.;        // exercise
    
    
    // pre_setted fitted params
    static double k1_H, k2_H, k3_H, k4_H, k5_H, k6_H, k7_H, k8_H, k9_H, k10_H, k11_H, k12_H, sigma_H, KM_H; // healthy
    static double k1_D1, k2_D1, k3_D1, k4_D1, k5_D1, k6_D1, k7_D1, k8_D1, k9_D1, k10_D1, k11_D1, k12_D1, sigma_D1, KM_D1; // D1
    static double k1_D2, k2_D2, k3_D2, k4_D2, k5_D2, k6_D2, k7_D2, k8_D2, k9_D2, k10_D2, k11_D2, k12_D2, sigma_D2, KM_D2; // D2
    static double k4e_H, k5e_H, k6e_H, k8e_H, lam_H;   // H, ex.
    static double k4e_D1, k5e_D1, k6e_D1, k8e_D1, lam_D1;   // D1, ex.
    static double k4e_D2, k5e_D2, k6e_D2, k8e_D2, lam_D2;   // D2, ex.
    
    // constant params in the model
    double gbliv = 0.043;
    double Gthpl = 9.;
    double vG = 17/70.;
    double vI = 13/70.;
    double beta = 1.;
    double fC = 1/180.16;
    double tau_i = 31.;
    double t_int = 30.;
    double tau_d = 3.;
    double c1 = 0.1;
    
    
    // ---------------- Vars used for checking evolution
    
    // time instants that are to be checked (Note: the first check_pt is always the initial time instant, t = 0.)
    std::vector<double> check_pts;
    // evolved parameters at the current time instant (updated after each evolution)
    double t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Jpl_curr, Iif_curr, Ex_curr, ExPre_curr;
    // output the obtained glucoses and insulins at the specified time instants in 'time_instants'
    std::vector<double> time_instants;
    std::vector<double> glucoses;
    std::vector<double> insulins;
    
    // ---------------- Vars used for parameter estimation
    //    double k1_f, k2_f, k3_f, k4_f, k5_f, k6_f, k7_f, k8_f, k9_f, k10_f, k11_f, k12_f, sigma_f, KM_f;
    std::vector<double> fittedParamsCurr;
    
};

#endif /* E_DES_Ex_hpp */
