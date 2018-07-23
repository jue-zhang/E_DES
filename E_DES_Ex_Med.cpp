//
//  E_DES_Ex_Med.cpp
//  E_DES
//
//  Created by Jue on 7/16/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "E_DES_Ex_Med.hpp"


// ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++

void E_DES_Ex_Med::SetEvolutionSpecifics(const std::tuple<int, double, double, double, double, double, double, int, int, int, double, double, double, double, int, double> &specifics){
    
    type            = std::get<0>(specifics);
    Mb              = std::get<1>(specifics);
    Gpl_basal       = std::get<2>(specifics);
    Ipl_basal       = std::get<3>(specifics);
    this_meal       = std::get<4>(specifics);
    latest_meal     = std::get<5>(specifics);
    delta_eat       = std::get<6>(specifics);
    ex_type         = std::get<7>(specifics);
    insulin_sa_type = std::get<8>(specifics);
    insulin_la_type_curr = std::get<9>(specifics);
    insulin_la_type_latest      = std::get<10>(specifics);
    insulin_la_curr = std::get<11>(specifics);
    insulin_la_latest = std::get<12>(specifics);
    delta_insulin   = std::get<13>(specifics);
    oralMed_type    = std::get<14>(specifics);
    oralMed         = std::get<15>(specifics);
}


void E_DES_Ex_Med::SetInitConditions(const std::vector<double> &init_conditions){
    if (init_conditions.size() != 9) {
        std::cout << "ERROR(E_DES_Ex::SetInitConditions): the size of the init_conditions params does NOT match!" << std::endl;
        return;
    }
    t_init      = init_conditions[0];
    MGgut_init  = init_conditions[1];
    Gpl_init    = init_conditions[2];
    Ipl_init    = init_conditions[3];
    Int_init    = init_conditions[4];
    ExPre_init  = init_conditions[5];
    Ex_init     = init_conditions[6];
    UI_sc1_init = init_conditions[7];
    UI_sc2_init = init_conditions[8];
    
    t_curr      = t_init;
    MGgut_curr  = MGgut_init;
    Gpl_curr    = Gpl_init;
    Ipl_curr    = Ipl_init;
    Int_curr    = Int_init;
    ExPre_curr  = ExPre_init;
    Ex_curr     = Ex_init;
    UI_sc1_curr = UI_sc1_init;
    UI_sc2_curr = UI_sc2_init;
    
}


void E_DES_Ex_Med::SetFittedParams(const int &type_input){
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
            k9 = k9_H;
            k10 = k10_H;
            k11 = k11_H;
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
            k9 = k9_D1;
            k10 = k10_D1;
            k11 = k11_D1;
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
            k9 = k9_D2;
            k10 = k10_D2;
            k11 = k11_D2;
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
            k9 = k9_H;
            k10 = k10_H;
            k11 = k11_H;
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

void E_DES_Ex_Med::SetFittedParams(const std::string &param_file_path){
    
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
        k9 = v_tmp[8];
        k10 = v_tmp[9];
        k11 = v_tmp[10];
        sigma = v_tmp[11];
        KM = v_tmp[12];
        k4e = v_tmp[13];        //exercise
        k5e = v_tmp[14];        //exercise
        k6e = v_tmp[15];        //exercise
        k8e = v_tmp[16];        //exercise
        lam = v_tmp[17];        //exercise
    }
    param_file.close();
    
}


void E_DES_Ex_Med::SetFittedParams(std::map<std::string, double> dict_fp){
    
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
    if (dict_fp.find("sigma") != dict_fp.end()) sigma = dict_fp["sigma"];
    if (dict_fp.find("KM") != dict_fp.end())    KM = dict_fp["KM"];
    // exercise
    if (dict_fp.find("k4e") != dict_fp.end())    KM = dict_fp["k4e"];
    if (dict_fp.find("k5e") != dict_fp.end())    KM = dict_fp["k5e"];
    if (dict_fp.find("k6e") != dict_fp.end())    KM = dict_fp["k6e"];
    if (dict_fp.find("k8e") != dict_fp.end())    KM = dict_fp["k8e"];
    if (dict_fp.find("lam") != dict_fp.end())    KM = dict_fp["lam"];
    
}

void E_DES_Ex_Med::SetCheckPts(double tI, double tF, int steps){
    if (tF < tI || steps < 0) {
        std::cout << "ERROR(E_DES_params::SetCheckPts)" << std::endl;
        return;
    }
    check_pts.clear();
    for (int i = 0; i <= steps; ++i) {
        check_pts.push_back(tI + i * (tF - tI) /steps);
    }
}

std::map<std::string, double> E_DES_Ex_Med::GetFittedParamsDict(){
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

void E_DES_Ex_Med::ClearPreRunOutput(){
    time_instants.clear();
    glucoses.clear();
    insulins.clear();
}


int E_DES_Ex_Med::Solver_gsl() {
    // ref. https://www.gnu.org/software/gsl/doc/html/ode-initval.html
    
    // -------- set-up the ODE solver
    const int para_nums = 8; // number of dynamic parames in ODEs,
    
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
        Int_curr = dynamic_vars[3];
        ExPre_curr = dynamic_vars[4];
        Ex_curr = dynamic_vars[5];
        UI_sc1_curr = dynamic_vars[6];
        UI_sc2_curr = dynamic_vars[7];
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

int E_DES_Ex_Med::gsl_ODEs (double t, const double y[], double f[], void *paramsP){
    (void)(t); /* avoid unused parameter warning */
    
    // passing the necessary parameters
    
    auto eDES_tmp = *static_cast<E_DES_Ex_Med *>(paramsP);
    auto specifics_tmp = eDES_tmp.GetEvolutionSpecifics();
    auto fittedParams_tmp = eDES_tmp.GetFittedParams();
    auto constParams_tmp = eDES_tmp.GetConstParams();
    //      <type, bodyMass, Gpl_basl, Ipl_basal, this_meal, latest_meal, delta_eat, ex_type, insulin_sa_type, insulin_la_type, insulin_sa, insulin_la_curr,  insulin_la_prev, delta_insulin, oralMed_type, oralMed>
    
    auto Mb                 = std::get<1>(specifics_tmp);
    auto Gpl_basal          = std::get<2>(specifics_tmp);
    auto Ipl_basal          = std::get<3>(specifics_tmp);
    auto this_meal          = std::get<4>(specifics_tmp);
    auto latest_meal        = std::get<5>(specifics_tmp);
    auto delta_eat          = std::get<6>(specifics_tmp);
    auto ex_type            = std::get<7>(specifics_tmp);
    auto insulin_sa_type    = std::get<8>(specifics_tmp);
    auto insulin_la_type_curr    = std::get<9>(specifics_tmp);
    auto insulin_la_type_latest        = std::get<10>(specifics_tmp);
    auto insulin_la_curr    = std::get<11>(specifics_tmp);
    auto insulin_la_prev    = std::get<12>(specifics_tmp);
    auto delta_insulin      = std::get<13>(specifics_tmp);
    auto oralMed_type       = std::get<14>(specifics_tmp);
    auto oralMed            = std::get<15>(specifics_tmp);
    
    auto k1     = fittedParams_tmp[0];
    auto k2     = fittedParams_tmp[1];
    auto k3     = fittedParams_tmp[2];
    auto k4     = fittedParams_tmp[3];
    auto k5     = fittedParams_tmp[4];
    auto k6     = fittedParams_tmp[5];
    auto k7     = fittedParams_tmp[6];
    auto k8     = fittedParams_tmp[7];
    auto k9     = fittedParams_tmp[8];
    auto k10    = fittedParams_tmp[9];
    auto k11    = fittedParams_tmp[10];
    
    auto sigma  = fittedParams_tmp[11];
    auto KM     = fittedParams_tmp[12];
    
    auto k4e    = fittedParams_tmp[13];     // ex
    auto k5e    = fittedParams_tmp[14];     // ex
    auto k6e    = fittedParams_tmp[15];     // ex
    auto k8e    = fittedParams_tmp[16];     // ex
    auto lam    = fittedParams_tmp[17];     // ex
    
    auto gbliv  = constParams_tmp[0];
    auto Gthpl  = constParams_tmp[1];
    auto vG     = constParams_tmp[2];
    auto vI     = constParams_tmp[3];
    auto beta   = constParams_tmp[4];
    auto fC     = constParams_tmp[5];
    auto tau_i  = constParams_tmp[6];
    auto t_int  = constParams_tmp[7];
    auto tau_d  = constParams_tmp[8];
    auto c1     = constParams_tmp[9];
    auto c2     = constParams_tmp[10];
    auto c3     = constParams_tmp[11];
    auto h_prev      = constParams_tmp[12];
    auto a_prev      = constParams_tmp[13];
    auto b_prev      = constParams_tmp[14];
    auto kd_prev     = constParams_tmp[15];
    auto h_curr      = constParams_tmp[16];
    auto a_curr      = constParams_tmp[17];
    auto b_curr      = constParams_tmp[18];
    auto kd_curr     = constParams_tmp[19];
    auto r1    = constParams_tmp[20];
    auto r2    = constParams_tmp[21];
    auto r3    = constParams_tmp[22];
    
    double Gbpl = Gpl_basal;
    double Ibpl = Ipl_basal;
    
    // variable mapping:
    //      y[0]: MGgut
    //      y[1]: Gpl
    //      y[2]: Ipl
    //      y[3]: Int
    //      y[4]: ExPre
    //      y[5]: Ex
    //      y[6]: UI_sc1
    //      y[7]: UI_sc2
    
//    std::cout << "y: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " " << y[6] << " " << y[7] << std::endl;
    
    // Glucose in the gut
    // Note: careful treatment on the previous meal
    double mGmeal_this = sigma * pow(k1, sigma) * pow(t, sigma-1) * exp(-pow(k1*t, sigma)) * exp(-k11*t) * this_meal;
    double mGmeal_latest = sigma * pow(k1, sigma) * pow(t+delta_eat, sigma-1) * exp(-pow(k1*(t+delta_eat), sigma)) * exp(-k11*(t+delta_eat)) * latest_meal;
    double mGpl = k2 * y[0];
    f[0] = mGmeal_this + mGmeal_latest - mGpl;
    
    // Glucose in the plasma
    double ggut = k2 * (fC/vG/Mb) * y[0];
    double gliv_up = gbliv - k3 * (y[1] - Gbpl) - k4 * (1 + k4e * y[5]) * beta * (y[2] - Ibpl);
    double gliv_down = gbliv - k10 * (y[1] - Gbpl) - k4 * (1 + k4e * y[5]) * beta * (y[2] - Ibpl);
    double gliv = (y[1] >= Gbpl)? gliv_up : gliv_down;
//    double gnonit_h = gbliv * (KM+Gbpl)/Gbpl * y[1] / (KM+y[1]);
//    double gnonit = gnonit_h;
    double gnonit = c2 * y[1] / (KM+y[1]);
    double git = k5 * (1 + k5e * y[5]) * beta * y[2] * y[1] / (KM + y[1]);
    double gren = (y[1] >= Gthpl)? (c1/vG/Mb) * (y[1] - Gthpl) : 0;
    f[1] = ggut + gliv - gnonit - git - gren;
    
    // Insulin in the plasma
    double ipnc = 1/beta* (k6*(1-k6e*y[5])*(y[1]-Gbpl) + k7/tau_i*y[3] + k7/tau_i*Gbpl + k8*(1-k8e*y[5])*tau_d*f[1]);
//    double iliv_h = k7*Gbpl*y[2]/(beta*tau_i*Ibpl);
//    double iliv = iliv_h;
    double iliv = c3 * y[2];
    double iif  = k9*(y[2] - Ibpl);
    double isa  = r3/vI/Mb*y[7];
    //      careful treatment on long-acting insulin
    double tHalf_curr = a_curr * insulin_la_curr + b_curr;
    double ila_curr = insulin_la_curr * (1-kd_curr) * h_curr*pow(tHalf_curr, h_curr)*pow(t, h_curr-1)/pow(pow(tHalf_curr, h_curr)+pow(t, h_curr), 2) / vI / Mb;
    double tHalf_prev = a_prev * insulin_la_prev + b_prev;
    double ila_prev = insulin_la_prev * (1-kd_prev) * h_prev*pow(tHalf_prev, h_prev)*pow(t+delta_insulin, h_prev-1)/pow(pow(tHalf_prev, h_prev)+pow(t+delta_insulin, h_prev), 2) / vI / Mb;
    double ila = ila_curr + ila_prev;
    f[2] = ipnc + isa + ila - iliv - iif;
    
    // Integrate part
    f[3] = (t < t_int)? y[1] - Gbpl : 0 ;
    
    // evolution of pre-existing ex.
    f[4] = - lam * y[4];
    // evolution of total ex.
    f[5] = - lam * y[4];
    
    // short-acting insulin
    f[6] = - r2 * y[6];
    f[7] = r2 * y[6] - r1 * y[7];
    
//    std::cout << "f:" << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << f[4] << " " << f[5] << " " << f[6] << " " << f[7] << std::endl;
    
    return GSL_SUCCESS;
}


// ++++++++++++++++++++++++++++++   Initialization of static memembers  ++++++++++++++++++++++++++++++

//// pre-setted fitted-params: healthy person (optimized 7/19/2018, based on Fig.14 of Mass' thesis)
//double E_DES_Ex_Med::k1_H = 0.0143161;
//double E_DES_Ex_Med::k2_H = 0.309044;
//double E_DES_Ex_Med::k3_H = 0.0171316;
//double E_DES_Ex_Med::k4_H = 0.00107114;
//double E_DES_Ex_Med::k5_H = 0.000537166;
//double E_DES_Ex_Med::k6_H = 0.251337;
//double E_DES_Ex_Med::k7_H = 0.0216814;
//double E_DES_Ex_Med::k8_H = 6.74236;
//double E_DES_Ex_Med::k9_H = 0.024347;
//double E_DES_Ex_Med::k10_H = 0.05;
//double E_DES_Ex_Med::k11_H = 0.;
//double E_DES_Ex_Med::sigma_H = 1.31931;
//double E_DES_Ex_Med::KM_H = 16.4845;
//double E_DES_Ex_Med::k4e_H = 0.;        // exercise
//double E_DES_Ex_Med::k5e_H = 10.;        // exercise
//double E_DES_Ex_Med::k6e_H = 0.05;        // exercise
//double E_DES_Ex_Med::k8e_H = 0.5;        // exercise
//double E_DES_Ex_Med::lam_H = 1/120.;        // exercise

// pre-setted fitted-params: healthy person (optimized 7/22/2018, based on Fig.17 of Mass' thesis)
double E_DES_Ex_Med::k1_H = 0.0163619;
double E_DES_Ex_Med::k2_H = 0.18896;
double E_DES_Ex_Med::k3_H = 0.0566132;
double E_DES_Ex_Med::k4_H = 0.00143478;
double E_DES_Ex_Med::k5_H = 0.000120797;
double E_DES_Ex_Med::k6_H = 0.550253;
double E_DES_Ex_Med::k7_H = 0.00749822;
double E_DES_Ex_Med::k8_H = 2.78914;
double E_DES_Ex_Med::k9_H = 0.0261024;
double E_DES_Ex_Med::k10_H = 0.05;
double E_DES_Ex_Med::k11_H = 0.;
double E_DES_Ex_Med::sigma_H = 1.22447;
double E_DES_Ex_Med::KM_H = 15.532;
double E_DES_Ex_Med::k4e_H = 0.;        // exercise
double E_DES_Ex_Med::k5e_H = 10.;        // exercise
double E_DES_Ex_Med::k6e_H = 0.05;        // exercise
double E_DES_Ex_Med::k8e_H = 0.5;        // exercise
double E_DES_Ex_Med::lam_H = 1/120.;        // exercise

// pre-setted fitted-params: D1
double E_DES_Ex_Med::k1_D1 = 0.0143161;
double E_DES_Ex_Med::k2_D1 = 0.309044;
double E_DES_Ex_Med::k3_D1 = 0.0171316;
double E_DES_Ex_Med::k4_D1 = 0.00107114;
double E_DES_Ex_Med::k5_D1 = 5.37166e-04;
double E_DES_Ex_Med::k6_D1 = 0.251337;
double E_DES_Ex_Med::k7_D1 = 2.16814e-05;
double E_DES_Ex_Med::k8_D1 = 1.0832;
double E_DES_Ex_Med::k9_D1 = 0.024347;
double E_DES_Ex_Med::k10_D1 = 0.05;
double E_DES_Ex_Med::k11_D1 = 0.;
double E_DES_Ex_Med::sigma_D1 = 1.31931;
double E_DES_Ex_Med::KM_D1 = 23.8407;
double E_DES_Ex_Med::k4e_D1 = 0.0;        // exercise
double E_DES_Ex_Med::k5e_D1 = 5.;        // exercise
double E_DES_Ex_Med::k6e_D1 = 0.1;        // exercise
double E_DES_Ex_Med::k8e_D1 = 0.0;        // exercise
double E_DES_Ex_Med::lam_D1 = 1/120.;        // exercise

// pre-setted fitted-params: D2   (optimized 7/19/2018, based on Fig.14 of Mass' thesis)
double E_DES_Ex_Med::k1_D2 = 0.0143161;
double E_DES_Ex_Med::k2_D2 = 0.309044;
double E_DES_Ex_Med::k3_D2 = 0.0171316;
double E_DES_Ex_Med::k4_D2 = 0.00107114;
double E_DES_Ex_Med::k5_D2 = 5.37166e-04;
double E_DES_Ex_Med::k6_D2 = 0.251337;
double E_DES_Ex_Med::k7_D2 = 2.16814e-05;
double E_DES_Ex_Med::k8_D2 = 1.0832;
double E_DES_Ex_Med::k9_D2 = 0.024347;
double E_DES_Ex_Med::k10_D2 = 0.05;
double E_DES_Ex_Med::k11_D2 = 0.;
double E_DES_Ex_Med::sigma_D2 = 1.31931;
double E_DES_Ex_Med::KM_D2 = 23.8407;
double E_DES_Ex_Med::k4e_D2 = 0.0;        // exercise
double E_DES_Ex_Med::k5e_D2 = 5.;        // exercise
double E_DES_Ex_Med::k6e_D2 = 0.1;        // exercise
double E_DES_Ex_Med::k8e_D2 = 0.0;        // exercise
double E_DES_Ex_Med::lam_D2 = 1/120.;        // exercise

// the specifices of long-acting insulin: <h, a, b, ke>
std::vector<double> E_DES_Ex_Med::insulin_la_specifics = {1.79, 2.88E-3, 5.61E2, 8.49E-1}; // Levemir
std::vector<double> E_DES_Ex_Med::insulin_sa_specifics = {0.248, 0.00805, 0.172};

