//
//  E_DES_Ex.cpp
//  E_DES
//
//  Created by Jue on 7/13/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "E_DES_Ex.hpp"


// ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++

void E_DES_Ex::SetEvolutionSpecifics(const std::tuple<int, double, double, double, double, double, double, int> &specifics){
    type        = std::get<0>(specifics);
    Mb          = std::get<1>(specifics);
    Gpl_basal   = std::get<2>(specifics);
    Ipl_basal   = std::get<3>(specifics);
    this_meal       = std::get<4>(specifics);
    latest_meal   = std::get<5>(specifics);
    delta_eat = std::get<6>(specifics);
    ex_type   = std::get<7>(specifics);
}


void E_DES_Ex::SetInitConditions(const std::vector<double> &init_conditions){
    if (init_conditions.size() != 8) {
        std::cout << "ERROR(E_DES_Ex::SetInitConditions): the size of the init_conditions params does NOT match!" << std::endl;
        return;
    }
    t_init      = init_conditions[0];
    MGgut_init  = init_conditions[1];
    Gpl_init    = init_conditions[2];
    Ipl_init    = init_conditions[3];
    Jpl_init    = init_conditions[4];
    Iif_init    = init_conditions[5];
    ExPre_init  = init_conditions[6];
    Ex_init     = init_conditions[7];
    
    t_curr      = t_init;
    MGgut_curr  = MGgut_init;
    Gpl_curr    = Gpl_init;
    Ipl_curr    = Ipl_init;
    Jpl_curr    = Jpl_init;
    Iif_curr    = Iif_init;
    ExPre_curr  = ExPre_init;
    Ex_curr     = Ex_init;
    
}

void E_DES_Ex::SetFittedParams(const std::vector<double> &fitted_params){
    if (fitted_params.size() != 18) {
        std::cout << "ERROR(E_DES_Ex::SetFittedParams): the size of the fitted params does NOT match!" << std::endl;
        return;
    }
    k1      = fitted_params[0];
    k2      = fitted_params[1];
    k3      = fitted_params[2];
    k4      = fitted_params[3];
    k5      = fitted_params[4];
    k6      = fitted_params[5];
    k7      = fitted_params[6];
    k8      = fitted_params[7];
    k9      = fitted_params[8];
    k10     = fitted_params[9];
    k11     = fitted_params[10];
    k12     = fitted_params[11];
    sigma   = fitted_params[12];
    KM      = fitted_params[13];
    
    k4e     = fitted_params[14];
    k5e     = fitted_params[15];
    k6e     = fitted_params[16];
    k8e     = fitted_params[17];
    lam     = fitted_params[18];
}

void E_DES_Ex::SetFittedParams(const int &type_input){
    switch (type_input) {
        case 0: // healthy person
            k1 = k1_H;
            k2 = k2_H;
            k3 = k3_H;
            k4 = k4_H;
            k5 = k5_H;
            k6 = k6_H;
            k7 = k7_H;
            k8 = k8_H;
            k9 = k9_H;         // short-acting insulin
            k10 = k10_H;        // short-acting insulin
            k11 = k11_H;
            k12 = k12_H;
            sigma = sigma_H;
            KM = KM_H;
            k4e = k4e_H;        //exercise
            k5e = k5e_H;        //exercise
            k6e = k6e_H;        //exercise
            k8e = k8e_H;        //exercise
            lam = lam_H;        //exercise
            break;
        case 1: // type-I
            k1 = k1_D1;
            k2 = k2_D1;
            k3 = k3_D1;
            k4 = k4_D1;
            k5 = k5_D1;
            k6 = k6_D1;
            k7 = k7_D1;
            k8 = k8_D1;
            k9 = k9_D1;         // short-acting insulin
            k10 = k10_D1;        // short-acting insulin
            k11 = k11_D1;
            k12 = k12_D1;
            sigma = sigma_D1;
            KM = KM_D1;
            k4e = k4e_D1;        //exercise
            k5e = k5e_D1;        //exercise
            k6e = k6e_D1;        //exercise
            k8e = k8e_D1;        //exercise
            lam = lam_D1;        //exercise
            break;
        case 2: // type-II
            k1 = k1_D2;
            k2 = k2_D2;
            k3 = k3_D2;
            k4 = k4_D2;
            k5 = k5_D2;
            k6 = k6_D2;
            k7 = k7_D2;
            k8 = k8_D2;
            k9 = k9_D2;         // short-acting insulin
            k10 = k10_D2;        // short-acting insulin
            k11 = k11_D2;
            k12 = k12_D2;
            sigma = sigma_D2;
            KM = KM_D2;
            k4e = k4e_D2;        //exercise
            k5e = k5e_D2;        //exercise
            k6e = k6e_D2;        //exercise
            k8e = k8e_D2;        //exercise
            lam = lam_D2;        //exercise
            break;
        default:// healthy person
            std::cout << "ERROR(E_DES::SetObjectTypeFittedParams): invalid input of subject type! Using the healthy person instead!" << std::endl;
            k1 = k1_H;
            k2 = k2_H;
            k3 = k3_H;
            k4 = k4_H;
            k5 = k5_H;
            k6 = k6_H;
            k7 = k7_H;
            k8 = k8_H;
            k9 = k9_H;         // short-acting insulin
            k10 = k10_H;        // short-acting insulin
            k11 = k11_H;
            k12 = k12_H;
            sigma = sigma_H;
            KM = KM_H;
            k4e = k4e_H;        //exercise
            k5e = k5e_H;        //exercise
            k6e = k6e_H;        //exercise
            k8e = k8e_H;        //exercise
            lam = lam_H;        //exercise
            break;
    }
}

void E_DES_Ex::SetFittedParams(const std::string &param_file_path){
    
    std::ifstream param_file;
    param_file.open(param_file_path, std::ifstream::in);
    
    double tmp = 0.;
    std::vector<double> v_tmp;
    
    while (param_file >> tmp) {
        v_tmp.push_back(tmp);
    }
    for (int i = 0; i < v_tmp.size(); ++i) {
        k1 = v_tmp[0];
        k2 = v_tmp[1];
        k3 = v_tmp[2];
        k4 = v_tmp[3];
        k5 = v_tmp[4];
        k6 = v_tmp[5];
        k7 = v_tmp[6];
        k8 = v_tmp[7];
        k9 = v_tmp[8];         // short-acting insulin
        k10 = v_tmp[9];        // short-acting insulin
        k11 = v_tmp[10];
        k12 = v_tmp[11];
        sigma = v_tmp[12];
        KM = v_tmp[13];
        k4e = v_tmp[14];        //exercise
        k5e = v_tmp[15];        //exercise
        k6e = v_tmp[16];        //exercise
        k8e = v_tmp[17];        //exercise
        lam = v_tmp[18];        //exercise
    }
    param_file.close();
    
}


void E_DES_Ex::SetFittedParams(std::map<std::string, double> dict_fp){
    
    if (dict_fp.find("k1") != dict_fp.end())    k1 = dict_fp["k1"];
    if (dict_fp.find("k2") != dict_fp.end())    k2 = dict_fp["k2"];
    if (dict_fp.find("k3") != dict_fp.end())    k3 = dict_fp["k3"];
    if (dict_fp.find("k4") != dict_fp.end())    k4 = dict_fp["k4"];
    if (dict_fp.find("k5") != dict_fp.end())    k5 = dict_fp["k5"];
    if (dict_fp.find("k6") != dict_fp.end())    k6 = dict_fp["k6"];
    if (dict_fp.find("k7") != dict_fp.end())    k7 = dict_fp["k7"];
    if (dict_fp.find("k8") != dict_fp.end())    k8 = dict_fp["k8"];
    if (dict_fp.find("k9") != dict_fp.end())    k9 = dict_fp["k9"];
    if (dict_fp.find("k10") != dict_fp.end())   k10 = dict_fp["k10"];
    if (dict_fp.find("k11") != dict_fp.end())   k11 = dict_fp["k11"];
    if (dict_fp.find("k12") != dict_fp.end())   k12 = dict_fp["k12"];
    if (dict_fp.find("sigma") != dict_fp.end()) sigma = dict_fp["sigma"];
    if (dict_fp.find("KM") != dict_fp.end())    KM = dict_fp["KM"];
    // exercise
    if (dict_fp.find("k4e") != dict_fp.end())    KM = dict_fp["k4e"];
    if (dict_fp.find("k5e") != dict_fp.end())    KM = dict_fp["k5e"];
    if (dict_fp.find("k6e") != dict_fp.end())    KM = dict_fp["k6e"];
    if (dict_fp.find("k8e") != dict_fp.end())    KM = dict_fp["k8e"];
    if (dict_fp.find("lam") != dict_fp.end())    KM = dict_fp["lam"];
    
}

void E_DES_Ex::SetCheckPts(double tI, double tF, int steps){
    if (tF < tI || steps < 0) {
        std::cout << "ERROR(E_DES_params::SetCheckPts)" << std::endl;
        return;
    }
    check_pts.clear();
    for (int i = 0; i <= steps; ++i) {
        check_pts.push_back(tI + i * (tF - tI) /steps);
    }
}

std::map<std::string, double> E_DES_Ex::GetFittedParamsDict(){
    std::map<std::string, double> ret = {
        {"k1", k1},
        {"k2", k2},
        {"k3", k3},
        {"k4", k4},
        {"k5", k5},
        {"k6", k6},
        {"k7", k7},
        {"k8", k8},
        {"k9", k9},
        {"k10", k10},
        {"k11", k11},
        {"k12", k12},
        {"sigma", sigma},
        {"KM", KM},
        // exercise
        {"k4e", k4e},
        {"k5e", k5e},
        {"k6e", k6e},
        {"k8e", k8e},
        {"lam", lam}
    };
    return ret;
}


// ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++

void E_DES_Ex::ClearPreRunOutput(){
    time_instants.clear();
    glucoses.clear();
    insulins.clear();
}


int E_DES_Ex::Solver_gsl() {
    // ref. https://www.gnu.org/software/gsl/doc/html/ode-initval.html
    
    // -------- set-up the ODE solver
    const int para_nums = 5 + 2; // number of dynamic parames in ODEs,                  //  add ex.
    
    // -------- set-up initial conditions
    double time_offset = check_pts[0]; // store the acutual initial time, as in performing the evolution t = 0
    double t = 0;
    double dynamic_vars[para_nums];
    std::vector<double> init_conditions = GetInitConditions();
    for(int i = 0; i < para_nums; ++i)
        dynamic_vars[i] = init_conditions[i+1];
    
    // -------- pass the instance to the static gsl_ODEs function
    gsl_odeiv2_system ode_sys = {gsl_ODEs, nullptr, para_nums, this}; // set the Jacobian to nullptr
    
    // -------- define the high-level wrapper, "driver", to solve ODEs
    // step function: gsl_odeiv2_step_rkf45
    double hstart = 1e-6; // the initial step size
    // error control: for each dynamic variable y, the desired error level
    //                  D_i = epsabs + epsrel * (a_y |y_i| + a_dydt h |y_i^\prime|)
    double epsabs = 1e-6; // the desired absolute error
    double epsrel = 0.; // the desired relative error
    gsl_odeiv2_driver * ode_driver = gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
    
    
    // -------- evolution; obtain results at the specified check_pts
    for (auto check_pt: check_pts){
        int status = gsl_odeiv2_driver_apply(ode_driver, &t, check_pt - time_offset, dynamic_vars);
        if (status != GSL_SUCCESS) {
            std::cout << "ERROR(E_DES_Solver): return value = " << status << std::endl;
            return status;
        }
        // store the evolved parameters at the current time
        t_curr = t + time_offset;
        MGgut_curr = dynamic_vars[0];
        Gpl_curr = dynamic_vars[1];
        Ipl_curr = dynamic_vars[2];
        Jpl_curr = dynamic_vars[3];
        Iif_curr = dynamic_vars[4];
        ExPre_curr = dynamic_vars[5];
        Ex_curr = dynamic_vars[6];
        // store the evolved glucoses and insulins for output
        time_instants.push_back(t_curr);
        glucoses.push_back(Gpl_curr);
        insulins.push_back(Ipl_curr);
//        // command line
//        std::cout << (t + time_offset) / 60.;
//        for(int i = 0; i < para_nums; ++i) std::cout << " " << dynamic_vars[i];
//        std::cout << std::endl;
    }
    
    gsl_odeiv2_driver_free(ode_driver);
    return GSL_SUCCESS;
    
}

int E_DES_Ex::gsl_ODEs (double t, const double y[], double f[], void *paramsP){
    (void)(t); /* avoid unused parameter warning */
    
    // passing the necessary parameters
    
    auto eDES_tmp = *static_cast<E_DES_Ex *>(paramsP);
    auto specifics_tmp = eDES_tmp.GetEvolutionSpecifics();
    auto fittedParams_tmp = eDES_tmp.GetFittedParams();
    auto constParams_tmp = eDES_tmp.GetConstParams();
    
    auto Mb     = std::get<1>(specifics_tmp);
    auto Gpl_basal  = std::get<2>(specifics_tmp);
    auto Ipl_basal = std::get<3>(specifics_tmp);
    auto this_meal = std::get<4>(specifics_tmp);
    auto latest_meal = std::get<5>(specifics_tmp);
    auto delta_eat = std::get<6>(specifics_tmp);
    auto ex_type = std::get<7>(specifics_tmp);
    
    auto k1     = fittedParams_tmp[0];
    auto k2     = fittedParams_tmp[1];
    auto k3     = fittedParams_tmp[2];
    auto k4     = fittedParams_tmp[3];
    auto k5     = fittedParams_tmp[4];
    auto k6     = fittedParams_tmp[5];
    auto k7     = fittedParams_tmp[6];
    auto k8     = fittedParams_tmp[7];
    //    auto k9     = fittedParams_tmp[8];
    //    auto k10    = fittedParams_tmp[9];
    auto k11    = fittedParams_tmp[10];
    auto k12    = fittedParams_tmp[11];
    auto sigma  = fittedParams_tmp[12];
    auto KM     = fittedParams_tmp[13];
    auto k4e    = fittedParams_tmp[14];     // ex
    auto k5e    = fittedParams_tmp[15];     // ex
    auto k6e    = fittedParams_tmp[16];     // ex
    auto k8e    = fittedParams_tmp[17];     // ex
    auto lam    = fittedParams_tmp[18];     // ex
    
    
    auto gbliv  = constParams_tmp[0];
    auto Gthpl  = constParams_tmp[1];
    auto vG     = constParams_tmp[2];
    //    auto vI     = constParams_tmp[3];
    auto beta   = constParams_tmp[4];
    auto fC     = constParams_tmp[5];
    auto tau_i  = constParams_tmp[6];
    auto t_int  = constParams_tmp[7];
    auto tau_d  = constParams_tmp[8];
    auto c1     = constParams_tmp[9];
    
    double Gbpl = Gpl_basal;
    double Ibpl = Ipl_basal;
    
    // variable mapping:
    //      y[0]: MGgut
    //      y[1]: Gpl
    //      y[2]: Ipl
    //      y[3]: Jpl
    //      y[4]: Iif
    //      y[5]: ExPre                                                                 // add ex.
    //      y[6]: Ex                                                                    // add ex.
    //
    // Glucose in the gut
    // Note: careful treatment on the previous meal
    double mGmeal_this = sigma * pow(k1, sigma) * pow(t, sigma-1) * exp(-pow(k1*t, sigma)) * this_meal;
    double mGmeal_latest = sigma * pow(k1, sigma) * pow(t+delta_eat, sigma-1) * exp(-pow(k1*(t+delta_eat), sigma)) * latest_meal;
    double mGpl = k2 * y[0];
    f[0] = mGmeal_this + mGmeal_latest - mGpl;
    // Glucose in the plasma
    double gliv = gbliv - k3 * (y[1] - Gbpl) - k4 * beta * y[4] * (1 + k4e * y[6]);     // add ex.
    double ggut = k2 * (fC/vG/Mb) * y[0];
    double gnonit = gbliv * (KM+Gbpl)/Gbpl * y[1] / (KM+y[1]);
    double git = k5 * beta * y[4] * y[1] / (KM + y[1]) * (1 + k5e * y[6]);              // add ex.
    double gren = (y[1] >= Gthpl)? (c1/vG/Mb) * (y[1] - Gthpl) : 0;
    f[1] = gliv + ggut - gnonit - git - gren;
    // Insulin in the plasma
    //    double integral_part = (t >= t_int)? (k7/tau_i)*(y[1]-Gbpl) : 0.;
    double integral_part = (t <= t_int)? (k7/tau_i)*(y[1]-Gbpl) : 0.;
    double ilivDiff = (k7/tau_i) * (Gbpl/beta) * y[3] / Ibpl;
    double iif = k11 * (y[2] - Ibpl);
    double iifDiff = k11 * y[3];
    f[2] = y[3];
    // Insulin in the interstitial fluid
    f[4] = iif - k12 * y[4];
    // evolution of pre-existing ex.                                                    // add ex.
    f[5] = - lam * y[5];                                                                // add ex.
    // evolution of total ex.                                                           // add ex.
    f[6] = - lam * y[5];                                                                // add ex.
    // Jp1 (Note: f[3] needs to put at last, as it involves f[4])
    double glivDiff = -k3 * f[1] - k4 * beta * f[4] * (1 + k4e * y[6]) - k4 * beta * y[4] * k4e * f[6];          // add ex.
    double ggutDiff = k2 * (fC/vG/Mb) * f[0];
    double gnonitDiff = gbliv * (KM+Gbpl)/Gbpl * (KM*f[1]) / pow(KM+y[1], 2);
    double gitDiff = k5 * beta * (KM*y[4]*f[1] + y[1]*f[4]*(KM+y[1])) / pow(KM+y[1], 2) * (1 + k5e * y[6])
                    + k5 * beta * y[4] * y[1] / (KM + y[1]) * k5e * f[6] ;              // add ex.
    double grenDiff = (y[1] >= Gthpl)? (c1/vG/Mb) * f[1] : 0;
    double GplDiff2 = glivDiff + ggutDiff - gnonitDiff - gitDiff - grenDiff;
    double ipncDiff = 1/beta * (k6 * f[1] * (1 - k6e*y[6] ) + k6 * (y[1]-Gbpl) * (-k6e*f[6]) + integral_part
                                + k8*tau_d*GplDiff2 * (1 - k8e * y[6])
                                + k8 * (-k8e * f[6]) * tau_d * f[1]);                   // add ex.
    f[3] = ipncDiff - ilivDiff - iifDiff;
    
//    std::cout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " " << y[6] << std::endl;
    
    return GSL_SUCCESS;
}


// ++++++++++++++++++++++++++++++   Initialization of static memembers  ++++++++++++++++++++++++++++++

// pre-setted fitted-params: healthy person (optimized: 7/7/2018) (add ex. 7/13/2018)
double E_DES_Ex::k1_H = 0.015262;
double E_DES_Ex::k2_H = 0.304526;
double E_DES_Ex::k3_H = 0.00808413;
double E_DES_Ex::k4_H = 0.000177229;
double E_DES_Ex::k5_H = 0.0798414;
double E_DES_Ex::k6_H = 0.28092;
double E_DES_Ex::k7_H = 0.0321147;
double E_DES_Ex::k8_H = 6.86228;
double E_DES_Ex::k9_H = 0.;         // short-acting insulin
double E_DES_Ex::k10_H = 0.;        // short-acting insulin
double E_DES_Ex::k11_H = 0.0271874;
double E_DES_Ex::k12_H = 0.294954;
double E_DES_Ex::sigma_H = 1.57349;
double E_DES_Ex::KM_H = 19.3604;
double E_DES_Ex::k4e_H = 0.5;        // exercise
double E_DES_Ex::k5e_H = 0.5;        // exercise
double E_DES_Ex::k6e_H = 0.5;        // exercise
double E_DES_Ex::k8e_H = 0.5;        // exercise
double E_DES_Ex::lam_H = 1/120.;        // exercise

// pre-setted fitted-params: D1 (currently, same as D2) (optimized: 7/7/2018) (add ex. 7/13/2018)
double E_DES_Ex::k1_D1 = 0.0157906;
double E_DES_Ex::k2_D1 = 0.105023;
double E_DES_Ex::k3_D1 = 0.00603478;
double E_DES_Ex::k4_D1 = 0.000202276;
double E_DES_Ex::k5_D1 = 0.00727093;
double E_DES_Ex::k6_D1 = 0.0826157;
double E_DES_Ex::k7_D1 = 0.00439569;
double E_DES_Ex::k8_D1 = 2.5473;
double E_DES_Ex::k9_D1 = 0.;         // short-acting insulin
double E_DES_Ex::k10_D1 = 0.;       // short-acting insulin
double E_DES_Ex::k11_D1 = 0.00877222;
double E_DES_Ex::k12_D1 = 0.0180231;
double E_DES_Ex::sigma_D1 = 1.4483;
double E_DES_Ex::KM_D1 = 23.9015;
double E_DES_Ex::k4e_D1 = 0.5;        // exercise
double E_DES_Ex::k5e_D1 = 0.5;        // exercise
double E_DES_Ex::k6e_D1 = 0.5;        // exercise
double E_DES_Ex::k8e_D1 = 0.5;        // exercise
double E_DES_Ex::lam_D1 = 1/120.;        // exercise

// pre-setted fitted-params: D2 (optimized: 7/7/2018) (add ex. 7/13/2018)
double E_DES_Ex::k1_D2 = 0.0157906;
double E_DES_Ex::k2_D2 = 0.105023;
double E_DES_Ex::k3_D2 = 0.00603478;
double E_DES_Ex::k4_D2 = 0.000202276;
double E_DES_Ex::k5_D2 = 0.00727093;
double E_DES_Ex::k6_D2 = 0.0826157;
double E_DES_Ex::k7_D2 = 0.00439569;
double E_DES_Ex::k8_D2 = 2.5473;
double E_DES_Ex::k9_D2 = 0.;         // short-acting insulin
double E_DES_Ex::k10_D2 = 0.;       // short-acting insulin
double E_DES_Ex::k11_D2 = 0.00877222;
double E_DES_Ex::k12_D2 = 0.0180231;
double E_DES_Ex::sigma_D2 = 1.4483;
double E_DES_Ex::KM_D2 = 23.9015;
double E_DES_Ex::k4e_D2 = 0.0;        // exercise
double E_DES_Ex::k5e_D2 = 0.5;        // exercise
double E_DES_Ex::k6e_D2 = 0.05;        // exercise
double E_DES_Ex::k8e_D2 = 0.0;        // exercise
double E_DES_Ex::lam_D2 = 1/120.;        // exercise

//// fitted-params for healthy person in original E_DES paper http://journals.sagepub.com/doi/abs/10.1177/1932296814562607
//double E_DES::k1_H = 1.45E-2;
//double E_DES::k2_H = 2.76E-1;
//double E_DES::k3_H = 6.07E-3;
//double E_DES::k4_H = 2.35E-4;
//double E_DES::k5_H = 9.49E-2;
//double E_DES::k6_H = 1.93E-1;
//double E_DES::k7_H = 1.15;
//double E_DES::k8_H = 7.27;
//double E_DES::k9_H = 0.;         // short-acting insulin
//double E_DES::k10_H = 0.;        // short-acting insulin
//double E_DES::k11_H = 3.83E-2;
//double E_DES::k12_H = 2.84E-1;
//double E_DES::sigma_H = 1.34;
//double E_DES::KM_H = 13.0995;

//// Manually tweaked fitted-params for healty person, comparing with those in original E_DES paper (6/25/2018)
//double E_DES::k1_H = 1.45E-2;
//double E_DES::k2_H = 2.76E-1;
//double E_DES::k3_H = 0.015; // tweaked params
//double E_DES::k4_H = 2.35E-4;
//double E_DES::k5_H = 0.06; // tweaked params
//double E_DES::k6_H = 1.; // tweaked params
//double E_DES::k7_H = 0.5; // tweaked params
//double E_DES::k8_H = 7.27;
//double E_DES::k9_H = 0.;         // short-acting insulin
//double E_DES::k10_H = 0.;        // short-acting insulin
//double E_DES::k11_H = 0.05; // tweaked params
//double E_DES::k12_H = 2.84E-1;
//double E_DES::sigma_H = 1.34;
//double E_DES::KM_H = 13.0995;

//// Manually tweaked fitted-params for D2, comparing with those in original E_DES paper (7/1/2018)
//double E_DES::k1_D2 = 1.45E-2;
//double E_DES::k2_D2 = 2.76E-1;
//double E_DES::k3_D2 = 0.005; // tweaked params
//double E_DES::k4_D2 = 2.35E-4;
//double E_DES::k5_D2 = 0.01; // tweaked params
//double E_DES::k6_D2 = 0.4; // tweaked params
//double E_DES::k7_D2 = 0.3; // tweaked params
//double E_DES::k8_D2 = 2.; // tweaked params
//double E_DES::k9_D2 = 0.;         // short-acting insulin
//double E_DES::k10_D2 = 0.;        // short-acting insulin
//double E_DES::k11_D2 = 0.05; // tweaked params
//double E_DES::k12_D2 = 2.84E-1;
//double E_DES::sigma_D2 = 1.34;
//double E_DES::KM_D2 = 16.65; // tweaked params


