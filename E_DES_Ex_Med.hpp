//
//  E_DES_Ex_Med.hpp
//  E_DES
//
//  Created by Jue on 7/16/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#ifndef E_DES_Ex_Med_hpp
#define E_DES_Ex_Med_hpp

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


class E_DES_Ex_Med{
    friend class SubjectGlucose;
    
    // ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++
public:
    // -------------  Modifying model paramers
    
    // set the vars. specific to the evolution:
    //      <type, bodyMass, Gpl_basl, Ipl_basal, this_meal, latest_meal, delta_eat, ex_type, insulin_sa_type, insulin_la_type, insulin_sa, insulin_la_curr,  insulin_la_prev, delta_insulin, oralMed_type, oralMed>
    void SetEvolutionSpecifics(const std::tuple<
                               int, double, double, double, double, double, double, int,
                               int, int, double, double, double, double, // insulin
                               int, double  // oralMed
                               > &specifics);
    
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
    std::tuple<int, double, double, double, double, double, double, int, int, int, int, double, double, double, int, double > GetEvolutionSpecifics() {
        return std::make_tuple(type, Mb, Gpl_basal, Ipl_basal, this_meal, latest_meal, delta_eat, ex_type,
                               insulin_sa_type, insulin_la_type_curr, insulin_la_type_latest, insulin_la_curr, insulin_la_latest, delta_insulin,
                               oralMed_type, oralMed);
    }
    
    std::vector<double> GetInitConditions() {
        return {t_init, MGgut_init, Gpl_init, Ipl_init, Int_init, ExPre_init, Ex_init, UI_sc1_init, UI_sc2_init};
    }
    std::vector<double> GetFittedParams() {
        return {k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, sigma, KM, k4e, k5e, k6e, k8e, lam};
    }
    std::vector<double> GetConstParams() { return {gbliv, Gthpl, vG, vI, beta, fC, tau_i, t_int, tau_d, c1, c2, c3, h_prev, a_prev, b_prev, kd_prev, h_curr, a_curr, b_curr, kd_curr, r1, r2, r3}; }
    std::vector<double> GetCurrentEvolvedParams() {
        return {t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Int_curr, ExPre_curr, Ex_curr, UI_sc1_curr, UI_sc2_curr};
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
    E_DES_Ex_Med() = default;
    
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
    int ex_type = 0;        // type of exercise
    int insulin_sa_type = 0;   // type of sa. insulin
    int insulin_la_type_curr = 0;   // type of current la. insulin
    int insulin_la_type_latest = 0; // type of the latest la. insulin
    double insulin_la_curr = 0.; // current long-acting insulin injection
    double insulin_la_latest = 0.; // latest long-acting insulin injection
    double delta_insulin = 0.; // the passed time since the latest insulin injection
    int oralMed_type = 0;  // type of oral
    double oralMed = 0.; // the intaked oral medicine
    
    // Initial_conditions
    double t_init = 0.;      // min
    double MGgut_init= 0.;  // mg
    double Gpl_init = 5.;   // mmol/L
    double Ipl_init = 8.;   // mU/L
    double Int_init = 0.;   // mU/L/min
    double ExPre_init = 0.; // pre-existing exercise intensity
    double Ex_init = 0.;    // total exercise intensity
    double UI_sc1_init = 0.; // mU
    double UI_sc2_init = 0.; // mU
    
    // fitted params: params that can obtained through fitting
    double k1 = 0.;
    double k2 = 0.;
    double k3 = 0.;
    double k4 = 0.;
    double k5 = 0.;
    double k6 = 0.;
    double k7 = 0.;
    double k8 = 0.;
    double k9 = 0.;
    double k10 = 0.;
    double k11 = 0.;
    double sigma = 1.57349;
    double KM = 19.3604;
    double k4e = 1.;        // exercise
    double k5e = 1.5;        // exercise
    double k6e = 0.5;        // exercise
    double k8e = 0.5;        // exercise
    double lam = 1/60.;        // exercise
    
    
    // pre_setted fitted params
    static double k1_H, k2_H, k3_H, k4_H, k5_H, k6_H, k7_H, k8_H, k9_H, k10_H, k11_H, sigma_H, KM_H; // healthy
    static double k1_D1, k2_D1, k3_D1, k4_D1, k5_D1, k6_D1, k7_D1, k8_D1, k9_D1, k10_D1, k11_D1, sigma_D1, KM_D1; // D1
    static double k1_D2, k2_D2, k3_D2, k4_D2, k5_D2, k6_D2, k7_D2, k8_D2, k9_D2, k10_D2, k11_D2, sigma_D2, KM_D2; // D2
    static double k4e_H, k5e_H, k6e_H, k8e_H, lam_H;   // H, ex.
    static double k4e_D1, k5e_D1, k6e_D1, k8e_D1, lam_D1;   // D1, ex.
    static double k4e_D2, k5e_D2, k6e_D2, k8e_D2, lam_D2;   // D2, ex.
    static std::vector<double> insulin_la_specifics;
    static std::vector<double> insulin_sa_specifics;
    
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
//    double c2 = 0.198788; // (optimized 7/19/2018, based on Fig.14 of Mass' thesis)
//    double c3 = 0.000198892; // (optimized 7/19/2018, based on Fig.14 of Mass' thesis)
    double c2 = 0.156152; // (optimized 8/2/2018, based on the literature data in EDES paper, using multinest)
    double c3 = 0.00885685; // (optimized 8/2/2018, based on the literature data in EDES paper, using multinest)
    //      long-acting insulin (depending on the brand)
    double h_prev = 1.79;   // Levemir
    double a_prev = 2.88E-3;
    double b_prev = 5.61E2;
    double kd_prev = 8.49E-1;
    double h_curr = 1.79;   // Levemir
    double a_curr = 2.88E-3;
    double b_curr = 5.61E3;
    double kd_curr = 8.49E-1;
    //      short-acting insulin
    double r1 = 0.2483;
    double r2 = 8.04891E-3;
    double r3 = 0.0764;
    
    
    // ---------------- Vars used for checking evolution
    
    // time instants that are to be checked (Note: the first check_pt is always the initial time instant, t = 0.)
    std::vector<double> check_pts;
    // evolved parameters at the current time instant (updated after each evolution)
    double t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Int_curr, Ex_curr, ExPre_curr, UI_sc1_curr, UI_sc2_curr;
    // output the obtained glucoses and insulins at the specified time instants in 'time_instants'
    std::vector<double> time_instants;
    std::vector<double> glucoses;
    std::vector<double> insulins;
    
    
    // ---------------- Vars used for parameter estimation
    std::vector<double> fittedParamsCurr;
    
};

#endif /* E_DES_Ex_Med_hpp */
