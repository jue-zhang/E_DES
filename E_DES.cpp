//
//  E_DES.cpp
//  E_DES
//
//  Created by Jue on 6/24/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "E_DES.hpp"

// input params
double E_DES::Dmeal = 75E3;
double E_DES::Mb = 75.;

// fitted params
////      healthy person (original)
//double E_DES::k1 = 1.45E-2;
//double E_DES::k2 = 2.76E-1;
//double E_DES::k3 = 6.07E-3;
//double E_DES::k4 = 2.35E-4;
//double E_DES::k5 = 9.49E-2;
//double E_DES::k6 = 1.93E-1;
//double E_DES::k7 = 1.15;
//double E_DES::k8 = 7.27;
//double E_DES::k9 = 0.;         // short-acting insulin
//double E_DES::k10 = 0.;        // short-acting insulin
//double E_DES::k11 = 3.83E-2;
//double E_DES::k12 = 2.84E-1;
//double E_DES::sigma = 1.34;
//double E_DES::KM = 13.0995;
////      healthy person (tweaked)
double E_DES::k1 = 1.45E-2;
double E_DES::k2 = 2.76E-1;
double E_DES::k3 = 0.015; // tweaked params
double E_DES::k4 = 2.35E-4;
double E_DES::k5 = 0.06; // tweaked params
double E_DES::k6 = 1.; // tweaked params
double E_DES::k7 = 0.5; // tweaked params
double E_DES::k8 = 7.27;
double E_DES::k9 = 0.;         // short-acting insulin
double E_DES::k10 = 0.;        // short-acting insulin
double E_DES::k11 = 0.05; // tweaked params
double E_DES::k12 = 2.84E-1;
double E_DES::sigma = 1.34;
double E_DES::KM = 13.0995;

// const params
double E_DES::gbliv = 0.043;
double E_DES::Gthpl = 9.;
double E_DES::vG = 17/70.;
double E_DES::vI = 13/70.;
double E_DES::beta = 1.;
double E_DES::fC = 1/180.16;
double E_DES::tau_i = 31.;
double E_DES::t_int = 30.;
double E_DES::tau_d = 3.;
double E_DES::c1 = 0.1;

// initial_conditions
double E_DES::t_init = 0.;      // min
double E_DES::MGgut_init= 0.;  // mg
double E_DES::Gpl_init = 5.;   // mmol/L
double E_DES::Ipl_init = 8.;   // U/L
double E_DES::Jpl_init = 0.;   // mU/L/min
double E_DES::Iif_init = 0.;   // mU/L

// evolved parameters at the current time
double E_DES::t_curr = 0.;      // min
double E_DES::MGgut_curr = 0.;  // mg
double E_DES::Gpl_curr = 5.;   // mmol/L
double E_DES::Ipl_curr = 8.;   // U/L
double E_DES::Jpl_curr = 0.;   // mU/L/min
double E_DES::Iif_curr = 0.;   // mU/L

// time_offset
double E_DES::time_offset = 0.;

void E_DES::SetInputParams(const std::vector<double> &input_params){
    if (input_params.size() != 2) {
        std::cout << "ERROR(E_DES_params::SetInputParams): the size of the input params does NOT match!" << std::endl;
        return;
    }
    Dmeal = input_params[0];
    Mb = input_params[1];
}

void E_DES::SetFittedParams(const std::vector<double> &fitted_params){
    if (fitted_params.size() != 14) {
        std::cout << "ERROR(E_DES_params::SetFittedParams): the size of the fitted params does NOT match!" << std::endl;
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
}

void E_DES::SetCheckPts(double tI, double tF, int steps){
    if (tF < tI || steps < 0) {
        std::cout << "ERROR(E_DES_params::SetCheckPts)" << std::endl;
        return;
    }
    for (int i = 0; i <= steps; ++i) {
        check_pts.push_back(tI + i * tF/steps);
    }
}

void E_DES::SetCheckPts(const std::vector<double> &check_pts_input){
    check_pts = check_pts_input;
}

void E_DES::SetInitConditions(const std::vector<double> &init_conditions){
    if (init_conditions.size() != 6) {
        std::cout << "ERROR(E_DES_params::SetInitConditions): the size of the init_conditions params does NOT match!" << std::endl;
        return;
    }
    t_init      = init_conditions[0];
    MGgut_init  = init_conditions[1];
    Gpl_init    = init_conditions[2];
    Ipl_init    = init_conditions[3];
    Jpl_init    = init_conditions[4];
    Iif_init    = init_conditions[5];
    
    t_curr      = t_init;
    MGgut_curr  = MGgut_init;
    Gpl_curr    = Gpl_init;
    Ipl_curr    = Ipl_init;
    Jpl_curr    = Jpl_init;
    Iif_curr    = Iif_init;
    
}

// Set the subject type and the corresponding fitted parameters
// Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
void E_DES::SetSubjectTypeFittedParams(const int &type){
    switch (type) {
        case 0: // healthy person
            k1 = 1.45E-2;
            k2 = 2.76E-1;
            k3 = 0.015; // tweaked params
            k4 = 2.35E-4;
            k5 = 0.06; // tweaked params
            k6 = 1.; // tweaked params
            k7 = 0.5; // tweaked params
            k8 = 7.27;
            k9 = 0.;         // short-acting insulin
            k10 = 0.;        // short-acting insulin
            k11 = 0.05; // tweaked params
            k12 = 2.84E-1;
            sigma = 1.34;
            KM = 13.0995;
            break;
        case 1: // type-I (Note: now it is the same as type-II)
            k1 = 1.45E-2;
            k2 = 2.76E-1;
            k3 = 0.005; // tweaked params
            k4 = 2.35E-4;
            k5 = 0.01; // tweaked params
            k6 = 0.4; // tweaked params
            k7 = 0.3; // tweaked params
            k8 = 2.;
            k9 = 0.;         // short-acting insulin
            k10 = 0.;        // short-acting insulin
            k11 = 0.05; // tweaked params
            k12 = 2.84E-1;
            sigma = 1.34;
            KM = 16.65; // tweaked params
            break;
        case 2: // type-II
            k1 = 1.45E-2;
            k2 = 2.76E-1;
            k3 = 0.005; // tweaked params
            k4 = 2.35E-4;
            k5 = 0.01; // tweaked params
            k6 = 0.4; // tweaked params
            k7 = 0.3; // tweaked params
            k8 = 2.;
            k9 = 0.;         // short-acting insulin
            k10 = 0.;        // short-acting insulin
            k11 = 0.05; // tweaked params
            k12 = 2.84E-1;
            sigma = 1.34;
            KM = 16.65; // tweaked params
            break;
        default:// healthy person
            std::cout << "ERROR(E_DES::SetObjectTypeFittedParams): invalid input of subject type! Using the healthy person instead!" << std::endl;
            k1 = 1.45E-2;
            k2 = 2.76E-1;
            k3 = 0.015; // tweaked params
            k4 = 2.35E-4;
            k5 = 0.06; // tweaked params
            k6 = 1.; // tweaked params
            k7 = 0.5; // tweaked params
            k8 = 7.27;
            k9 = 0.;         // short-acting insulin
            k10 = 0.;        // short-acting insulin
            k11 = 0.05; // tweaked params
            k12 = 2.84E-1;
            sigma = 1.34;
            KM = 13.0995;
            break;
    }
}


void E_DES::ClearPreRuns(){
    check_pts.clear();
    time_instants.clear();
    glucoses.clear();
    insulins.clear();
}
std::vector<double> E_DES::GetInputParams() {
    std::vector<double> inputParams = {Dmeal, Mb};
    return inputParams;
}

std::vector<double> E_DES::GetFittedParams() {
    std::vector<double> fittedParams = {k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, sigma, KM};
    return fittedParams;
}

std::vector<double> E_DES::GetConstParams() {
    std::vector<double> constParams = {gbliv, Gthpl, vG, vI, beta, fC, tau_i, t_int, tau_d, c1};
    return constParams;
}

std::vector<double> E_DES::GetInitConditions() {
    std::vector<double> initConditions = {t_init, MGgut_init, Gpl_init, Ipl_init, Jpl_init, Iif_init};
    return initConditions;
}

std::vector<double> E_DES::GetCurrentEvolvedParams() {
    std::vector<double> currentEvolvedParams = {t_curr, MGgut_curr, Gpl_curr, Ipl_curr, Jpl_curr, Iif_curr};
    return currentEvolvedParams;
}

int E_DES::gsl_ODEs (double t, const double y[], double f[], void *paramsP){
    (void)(t); /* avoid unused parameter warning */

    double Gbpl = Gpl_init;
    double Ibpl = Ipl_init;
    
    // -------- set-up ODEs
    // y[0]: MGgut
    // y[1]: Gpl
    // y[2]: Ipl
    // y[3]: Jpl
    // y[4]: Iif
    //
    // Glucose in the gut
    double mGmeal = sigma * pow(k1, sigma) * pow(t, sigma-1) * exp(-pow(k1*t, sigma)) * Dmeal;
    double mGpl = k2 * y[0];
    f[0] = mGmeal - mGpl;
    // Glucose in the plasma
    double gliv = gbliv - k3 * (y[1] - Gbpl) - k4 * beta * y[4];
    double ggut = k2 * (fC/vG/Mb) * y[0];
    double gnonit = gbliv * (KM+Gbpl)/Gbpl * y[1] / (KM+y[1]);
    double git = k5 * beta * y[4] * y[1] / (KM + y[1]);
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
    // Jp1 (Note: f[3] needs to put at last, as it involves f[4])
    double glivDiff = -k3 * f[1] - k4 * beta * f[4];
    double ggutDiff = k2 * (fC/vG/Mb) * f[0];
    double gnonitDiff = gbliv * (KM+Gbpl)/Gbpl * (KM*f[1]) / pow(KM+y[1], 2);
    double gitDiff = k5 * beta * (KM*y[4]*f[1] + y[1]*f[4]*(KM+y[1])) / pow(KM+y[1], 2);
    double grenDiff = (y[1] >= Gthpl)? (c1/vG/Mb) * f[1] : 0;
    double GplDiff2 = glivDiff + ggutDiff - gnonitDiff - gitDiff - grenDiff;
    double ipncDiff = 1/beta * (k6*f[1] + integral_part + k8*tau_d*GplDiff2);
    f[3] = ipncDiff - ilivDiff - iifDiff;
    
    return GSL_SUCCESS;
}


int E_DES::Solver_gsl() {
    
    // -------- set-up the ODE solver
    const int para_nums = 5; // number of dynamic parames in ODEs
    gsl_odeiv2_system ode_sys = {gsl_ODEs, nullptr, para_nums, nullptr}; // set the Jacobian to nullptr
    
    // -------- define the high-level wrapper, "driver", to solve ODEs
    // step function: gsl_odeiv2_step_rkf45
    double hstart = 1e-6; // the initial step size
    // error control: for each dynamic variable y, the desired error level
    //                  D_i = epsabs + epsrel * (a_y |y_i| + a_dydt h |y_i^\prime|)
    double epsabs = 1e-6; // the desired absolute error
    double epsrel = 0.; // the desired relative error
    gsl_odeiv2_driver * ode_driver = gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
    
    // -------- set-up initial conditions
    double evol_var = check_pts[0];
    double dynamic_vars[para_nums];
    std::vector<double> init_conditions = GetInitConditions();
    for(int i = 0; i < para_nums; ++i)
        dynamic_vars[i] = init_conditions[i+1];
    
    // -------- evolution; obtain results at the specified check_pts
    for (auto check_pt: check_pts){
        int status = gsl_odeiv2_driver_apply(ode_driver, &evol_var, check_pt, dynamic_vars);
        if (status != GSL_SUCCESS) {
            std::cout << "ERROR(E_DES_Solver): return value = " << status << std::endl;
            return status;
        }
        // store the evolved parameters at the current time
        t_curr = evol_var + time_offset;
        MGgut_curr = dynamic_vars[0];
        Gpl_curr = dynamic_vars[1];
        Ipl_curr = dynamic_vars[2];
        Jpl_curr = dynamic_vars[3];
        Iif_curr = dynamic_vars[4];
        // store the evolved glucoses and insulins for output
        time_instants.push_back(t_curr);
        glucoses.push_back(Gpl_curr);
        insulins.push_back(Ipl_curr);
//        // command line
//        std::cout << evol_var;
//        for(int i = 0; i < para_nums; ++i)
//            std::cout << " " << dynamic_vars[i];
//        std::cout << std::endl;
    }
    
    gsl_odeiv2_driver_free(ode_driver);
    return GSL_SUCCESS;
    
}

int E_DES::Solver_gsl(std::ofstream &check_output) {
    
    // -------- set-up the ODE solver
    const int para_nums = 5; // number of dynamic parames in ODEs
    gsl_odeiv2_system ode_sys = {gsl_ODEs, nullptr, para_nums, nullptr}; // set the Jacobian to nullptr
    
    // -------- define the high-level wrapper, "driver", to solve ODEs
    // step function: gsl_odeiv2_step_rkf45
    double hstart = 1e-6; // the initial step size
    // error control: for each dynamic variable y, the desired error level
    //                  D_i = epsabs + epsrel * (a_y |y_i| + a_dydt h |y_i^\prime|)
    double epsabs = 1e-6; // the desired absolute error
    double epsrel = 0.; // the desired relative error
    gsl_odeiv2_driver * ode_driver = gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
    
    // -------- set-up initial conditions
    double evol_var = check_pts[0];
    double dynamic_vars[para_nums];
    std::vector<double> init_conditions = GetInitConditions();
    for(int i = 0; i < para_nums; ++i)
        dynamic_vars[i] = init_conditions[i+1];
    
    // -------- evolution; obtain results at the specified check_pts
    for (auto check_pt: check_pts){
        int status = gsl_odeiv2_driver_apply(ode_driver, &evol_var, check_pt, dynamic_vars);
        if (status != GSL_SUCCESS) {
            std::cout << "ERROR(E_DES_Solver): return value = " << status << std::endl;
            return status;
        }
        // store the evolved parameters at the current time
        t_curr = evol_var + time_offset;
        MGgut_curr = dynamic_vars[0];
        Gpl_curr = dynamic_vars[1];
        Ipl_curr = dynamic_vars[2];
        Jpl_curr = dynamic_vars[3];
        Iif_curr = dynamic_vars[4];
        // store the evolved glucoses and insulins
        time_instants.push_back(t_curr);
        glucoses.push_back(Gpl_curr);
        insulins.push_back(Ipl_curr);
        // output to the command line and file
        // command line
        std::cout << evol_var;
        for(int i = 0; i < para_nums; ++i)
            std::cout << " " << dynamic_vars[i];
        std::cout << std::endl;
        // file
        check_output << evol_var;
        for(int i = 0; i < para_nums; ++i)
            check_output << " " << dynamic_vars[i];
        check_output << std::endl;
    }
    
    gsl_odeiv2_driver_free(ode_driver);
    
    return GSL_SUCCESS;
    
}

double E_DES::TwoHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input){
    ClearPreRuns();
    std::vector<double> inputParams = {foodIntake, bodyMass};
    SetInputParams(inputParams);
    std::vector<double> initialConditions = {0., 0., Gpl_init_input, Ipl_init_input, 0., 0.};
    SetInitConditions(initialConditions);
    std::vector<double> checkPtsInput = {0., 120.};
    SetCheckPts(checkPtsInput);
    Solver_gsl();
    return glucoses[1];
}

std::vector<double> E_DES::FourHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input){
    ClearPreRuns();
    std::vector<double> inputParams = {foodIntake, bodyMass};
    SetInputParams(inputParams);
    std::vector<double> initialConditions = {0., 0., Gpl_init_input, Ipl_init_input, 0., 0.};
    SetInitConditions(initialConditions);
    SetCheckPts(0, 240., 24);
    Solver_gsl();
    std::vector<double> ret;
    for (auto iter = glucoses.begin() + 1; iter != glucoses.end(); ++iter) { // skip glucoses[0] (initial value)
        if (*iter < Gpl_init_input) ret.push_back(Gpl_init_input);
        else ret.push_back(*iter);
    }
    return ret;
}

std::vector<double> E_DES::EightHourGlucose(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input){
    ClearPreRuns();
    std::vector<double> inputParams = {foodIntake, bodyMass};
    SetInputParams(inputParams);
    std::vector<double> initialConditions = {0., 0., Gpl_init_input, Ipl_init_input, 0., 0.};
    SetInitConditions(initialConditions);
    SetCheckPts(0, 480., 48);
    Solver_gsl();
    std::vector<double> ret;
    for (auto iter = glucoses.begin() + 1; iter != glucoses.end(); ++iter) { // skip glucoses[0] (initial value)
        ret.push_back(*iter);
//        if (*iter < Gpl_init_input) ret.push_back(Gpl_init_input);
//        else ret.push_back(*iter);
    }
    return ret;
}

std::vector<double> E_DES::EightHourGlucosePerMin(double foodIntake, double bodyMass, double Gpl_init_input, double Ipl_init_input){
    ClearPreRuns();
    std::vector<double> inputParams = {foodIntake, bodyMass};
    SetInputParams(inputParams);
    std::vector<double> initialConditions = {0., Gpl_init_input, Ipl_init_input, 0., 0.};
    SetInitConditions(initialConditions);
    SetCheckPts(0, 480., 480);
    Solver_gsl();
    std::vector<double> ret;
    for (auto iter = glucoses.begin() + 1; iter != glucoses.end(); ++iter) { // skip glucoses[0] (initial value)
        ret.push_back(*iter);
//        if (*iter < Gpl_init_input) ret.push_back(Gpl_init_input);
//        else ret.push_back(*iter);
    }
    return ret;
}


std::vector<std::pair<double, double>> E_DES::GlucoseUnderFoodIntakeEvents(double bodyMass, double Gpl_init_input, double Ipl_init_input, std::vector<std::pair<double, double>> foodIntakeEvents){
    // check if there are at least two elements in 'foodIntakeEvents', one for the initial time, and the last one for
    //  the final time corresponding to the time that is 8 hours after the final non-zero food intake.
    if (foodIntakeEvents.size() < 2) {
        std::cout << "ERROR(E_DES::GlucoseUnderFoodIntakeEvents): the input size of 'foodIntakeEvents' must be larger than 1" << std::endl;
    }
    std::vector<double> initialConditions = {foodIntakeEvents[0].first, 0., Gpl_init_input, Ipl_init_input, 0., 0.};
    SetInitConditions(initialConditions);
    std::vector<std::pair<double, double>> ret;
    for (auto event = foodIntakeEvents.begin(); event != foodIntakeEvents.end() - 1; ++event) {
        // set initial conditions of the current evolution as the evolved parameters from the last evolution
        ClearPreRuns();
        double foodIntake_tmp = event->second;
        std::vector<double> inputParams_tmp = {foodIntake_tmp, bodyMass};
        SetInputParams(inputParams_tmp);
        std::vector<double> initialConditions_tmp = GetCurrentEvolvedParams();
        SetInitConditions(initialConditions_tmp);
        // set check_pts in the current evolution
        double tI_tmp = t_curr;
        double tF_tmp = (event+1)->first; // the time instant of the next food intake event
        time_offset = tI_tmp; // set time_offset
        std::vector<double> check_pts_tmp;
        // careful treatment when 'tI_tmp' and 'tF_tmp' are not integers of 10's
        double tI_tmp_10 = floor(tI_tmp/10.) * 10.;
        if ( fabs(tI_tmp - tI_tmp_10) > 1E-5 ) // manually include tI_tmp if it is not 10's
            check_pts_tmp.push_back(tI_tmp - time_offset);
        double tF_tmp_10 = floor(tF_tmp/10.) * 10.;
        int steps_tmp = static_cast<int>( (tF_tmp_10 - tI_tmp_10) / 10.);
        for (int i = 0; i <= steps_tmp; ++i) check_pts_tmp.push_back(tI_tmp_10 + i * 10. - time_offset);
        if ( fabs(tF_tmp - tF_tmp_10) > 1E-5 ) // manually include tF_tmp if it is not 10's
            check_pts_tmp.push_back(tF_tmp - time_offset);
        SetCheckPts(check_pts_tmp);
        // evolve
        Solver_gsl();
        // export the glucose levels at the time instants that are multiple of 10's
        for (std::size_t i = 0; i < time_instants.size(); ++i) {
            double time_res = fmod(time_instants[i], 10.);
            if ( fabs(time_res) < 1E-5 ) ret.push_back({time_instants[i], glucoses[i]});
        }
    }
    return ret;
}




void E_DES::SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_files,
                                          const std::vector<std::string> &dpe_insulin_files){
    
    std::ifstream ifile_g, ifile_i;
    if (dpe_glucose_files.size() != dpe_insulin_files.size()) {
        std::cout << "ERROR(E_DES::SetDataForParameterEstimation): the number of glucose and insulin files does NOT match!" << std::endl;
        return;
    }
    std::vector<double> times;
    std::vector<double> glucoses;
    std::vector<double> insulins;
    data_set param_est_data_set;
    for (auto i = 0; i < dpe_glucose_files.size(); ++i) {
        ifile_g.open(dpe_glucose_files[i], std::ifstream::in);
        ifile_i.open(dpe_insulin_files[i], std::ifstream::in);
        std::string head;
        getline(ifile_g, head);
        getline(ifile_i, head);
        times.clear();
        glucoses.clear();
        insulins.clear();
        param_est_data_set.clear();
        double time, glucose, insulin;
        while (ifile_g >> time >> glucose) {
            times.push_back(time);
            glucoses.push_back(glucose);
        }
        while (ifile_i >> time >> insulin) {
            insulins.push_back(insulin);
        }
        param_est_data_set.push_back(times);
        param_est_data_set.push_back(glucoses);
        param_est_data_set.push_back(insulins);
        param_est_data_sets.push_back(param_est_data_set);
        ifile_g.close();
        ifile_i.close();
    }
    
}

//void E_DES::MultiNest_LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context){
//    
//    // for k1:k8, k11, k12, sigma, scatter pts within [10^-5, 10^2]
//    k1 = pow(10, -5 + Cube[0] * 7);
//    k2 = pow(10, -5 + Cube[1] * 7);
//    k3 = pow(10, -5 + Cube[2] * 7);
//    k4 = pow(10, -5 + Cube[3] * 7);
//    k5 = pow(10, -5 + Cube[4] * 7);
//    k6 = pow(10, -5 + Cube[5] * 7);
//    k7 = pow(10, -5 + Cube[6] * 7);
//    k8 = pow(10, -5 + Cube[7] * 7);
//    //k9
//    //k10
//    k11 = pow(10, -5 + Cube[8] * 7);
//    k12 = pow(10, -5 + Cube[9] * 7);
//    sigma = pow(10, -5 + Cube[10] * 7);
//    // for KM, scatter pts within [1, 10^3]
//    KM = pow(10, Cube[11]*3);
//    
//    std::vector<double> fitted_params;
//    const int fp_sz_1 = 8; // k1 : k8
//    for (int i = 0; i < fp_sz_1; ++i)
//        fitted_params.push_back(pow(10, -5+Cube[i]*7)); // for k1:k8, scatter pts within [10^-5, 10^2]
//    fitted_params.push_back(0); // k9
//    fitted_params.push_back(0); // k10
//    const int fp_sz_2 = 3; // k11, k12, sigma
//    for (int i = 0; i < fp_sz_2; ++i)
//        fitted_params.push_back(pow(10, -5+Cube[fp_sz_1-1+i]*7)); // for k11, k12, sigma, scatter pts within [10^-5, 10^2]
//    fitted_params.push_back(pow(10, Cube[fp_sz_1+fp_sz_2]*3)); // KM within [1, 10^3]
//    
//    SetFittedParams(fitted_params);
//    
//    double SSR_tot = 0.;
//    for (const auto &param_est_data_set: param_est_data_sets) {
//        check_pts = param_est_data_set[0];
//        Solver_gsl();
//        // calculate SSR
//        double SSR_tmp = 0.;
//        
//    }
//    
//}

//void E_DES::EstimateParameters(){
//    
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
//    int ndims = 3;					// dimensionality (no. of free parameters)
//    
//    int nPar = 4;					// total no. of parameters including free & derived parameters
//    
//    int nClsPar = 3;				// no. of parameters to do mode separation on
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
//    char root[100] = "chains/r_0p05_";			// root for output files
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
//    
//    // calling MultiNest
//    nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, MultiNest_LogLike, MultiNest_Dumper, context);
//
//}



