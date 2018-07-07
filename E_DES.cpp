//
//  E_DES.cpp
//  E_DES
//
//  Created by Jue on 6/24/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "E_DES.hpp"

// ++++++++++++++++++++++++++++++   High level methods  ++++++++++++++++++++++++++++++

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

std::vector<double> E_DES::EightHourGlucosePerMin(double foodIntake, double bodyMass,
                                                  double Gpl_init_input, double Ipl_init_input){
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

// ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++


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

void E_DES::SetCheckPts(const std::vector<double> &check_pts_input){
    check_pts = check_pts_input;
}

// Set the subject type and the corresponding fitted parameters
// Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
void E_DES::SetSubjectTypeFittedParams(const int &type_input){
    type = type_input;
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
            break;
    }
}

// load the optimized fitted parameters
void E_DES::LoadFittedParams(std::ifstream &param_file){
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
    }
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

std::vector<std::pair<double, double>> E_DES::GetGlucoses(){
    std::vector<std::pair<double, double>> ret;
    for (int i = 0; i < time_instants.size(); ++i) {
        ret.push_back({time_instants[i], glucoses[i]});
    }
    return ret;
}

// ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++

void E_DES::ClearPreRuns(){
    check_pts.clear();
    time_instants.clear();
    glucoses.clear();
    insulins.clear();
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
//        for(int i = 0; i < para_nums; ++i) std::cout << " " << dynamic_vars[i];
//        std::cout << std::endl;
    }
    
    gsl_odeiv2_driver_free(ode_driver);
    return GSL_SUCCESS;
    
}

int E_DES::gsl_ODEs (double t, const double y[], double f[], void *paramsP){
    (void)(t); /* avoid unused parameter warning */

    double Gbpl = Gpl_init;
    double Ibpl = Ipl_init;
    
    // variable mapping:
    //      y[0]: MGgut
    //      y[1]: Gpl
    //      y[2]: Ipl
    //      y[3]: Jpl
    //      y[4]: Iif
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

// ++++++++++++++++++++++++++++++   Methods for estimating fitted-params  ++++++++++++++++++++++++++++++

void E_DES::SetDataForParameterEstimation(const std::vector<std::string> &dpe_glucose_insulin_files,
                                   const std::vector<std::vector<double>> &input_parameter_sets){
    input_param_sets = input_parameter_sets;
    std::ifstream ifile_gi;
    data_set param_est_data_set;
    std::vector<double> param_est_data_set_row;
    double ti, glu_i, glu_err_i, ins_i, ins_err_i;
    for (auto i = 0; i < dpe_glucose_insulin_files.size(); ++i) {
        param_est_data_set.clear();
        ifile_gi.open(dpe_glucose_insulin_files[i], std::ifstream::in);
        while (ifile_gi >> ti >> glu_i >> glu_err_i >> ins_i >> ins_err_i){
            param_est_data_set_row.clear();
            param_est_data_set_row.push_back(ti);
            param_est_data_set_row.push_back(glu_i);
            param_est_data_set_row.push_back(glu_err_i);
            param_est_data_set_row.push_back(ins_i);
            param_est_data_set_row.push_back(ins_err_i);
            param_est_data_set.push_back(param_est_data_set_row);
        }
        param_est_data_sets.push_back(param_est_data_set);
        ifile_gi.close();
    }
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


void E_DES::EstimateFittedParameters_gsl(){ // ref. of using GSL: https://www.gnu.org/software/gsl/doc/html/multimin.html
    
    // set up the initial point and initial step size
    minParams_init = {k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, sigma, KM};
    const int num_params = 14;
    
    // search in the parameter space of log10(ki)
    double interval = 2.;
    std::vector<std::pair<double, double>> log10_min_params_range;
    // default treatment on the lower and higher bounds of the minimized params
    for (int i = 0; i < minParams_init.size(); ++i) {
        double tmp = log10(minParams_init[i]);
        log10_min_params_range.push_back( {tmp-interval, tmp+interval} );
    }
    // special treatment on some parameters
    log10_min_params_range[0] = { log10(0.5*k1), log10(2.*k1) }; // k1
    log10_min_params_range[1] = { log10(0.5*k2), log10(2.*k2) }; // k2
    log10_min_params_range[8] = {log10(0.99*1E-10), log10(1E-10)}; // k9
    log10_min_params_range[9] = {log10(0.99*1E-10), log10(1E-10)}; // k10
    log10_min_params_range[12] = {log10(1.), log10(2.)}; // sigma
    log10_min_params_range[13] = {log10(5.), log10(30.)}; // KM
    
    // additional type-II constraints on parameters obtained from healthy persons
    // Note: only used when obtaining the fitted parameter of D2 from those of H
    log10_min_params_range[4] = { log10(0.01*k5), log10(k5) }; // k5
    log10_min_params_range[5] = { log10(0.01*k6), log10(k6) }; // k6
    log10_min_params_range[6] = { log10(0.01*k7), log10(k7) }; // k7
    log10_min_params_range[7] = { log10(0.01*k8), log10(k8) }; // k8
    log10_min_params_range[13] = {log10(KM), log10(2*KM)}; // KM
    
    // set up the minimized params (gsl_vector) used in the gsl subroutine of minimization
    gsl_vector *initial_step_size, *min_params; // 'min_params' -- params to be performed minimizations
    min_params = gsl_vector_alloc (num_params);
    initial_step_size = gsl_vector_alloc (num_params);
    minParams_fval_curr.clear();
    std::vector<double> SS, SD;
    // parameter transformation: in order to satisfying the constraints that all params > 0
    //      ref.: http://cafim.sssup.it/~giulio/software/multimin/multimin.html
    for (int i = 0; i < log10_min_params_range.size(); ++i) {
        double SS_tmp = (log10_min_params_range[i].second + log10_min_params_range[i].first) / 2.;
        double SD_tmp = (log10_min_params_range[i].second - log10_min_params_range[i].first) / 2.;
        SS.push_back( SS_tmp );
        SD.push_back( SD_tmp );
        double x = SS_tmp;
        double y = atanh( (x - SS_tmp) / SD_tmp );
        gsl_vector_set (min_params, i, y);
        gsl_vector_set (initial_step_size, i, 0.01);
        minParams_fval_curr.push_back( minParams_init[i] );
    }
    minParams_fval_curr.push_back(0.); // initial final val
    
    // set up the minimizer
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; // method used for minimization: Nelder-Mead
    gsl_multimin_fminimizer *s = nullptr; // initial the minimizer
    gsl_multimin_function min_func; // 'min_func' -- functions to be minimized over
    min_func.n = num_params;
    min_func.f = gsl_min_fitted_params_SSR_func;
    std::pair<std::vector<double>, std::vector<double>> params_tmp = {SS, SD};
    min_func.params = &params_tmp;
    
    s = gsl_multimin_fminimizer_alloc (T, num_params);
    gsl_multimin_fminimizer_set (s, &min_func, min_params, initial_step_size);
    
    // set up stopping criteria
    std::size_t iter_max = 1000; // max number of iterations
    std::size_t iter = 0;  // current iter
    int status;
    double SSR_curr = 0.;
    double epsabs = 10.;
    
    // iterate:
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) // break when error occurs
            break;
        
        SSR_curr = s->fval;
        status = gsl_multimin_test_size (SSR_curr, epsabs);
        
        minParams_fval_curr.clear();
        std::cout << iter << " ";
        for (int i = 0; i < num_params; ++i) {
            double y = gsl_vector_get (s->x, i);
            double original_param = pow(10, SS[i] + SD[i] * tanh(y) );
            minParams_fval_curr.push_back (original_param);
            std::cout << original_param << " ";
        }
        minParams_fval_curr.push_back(SSR_curr);
        std::cout << SSR_curr << std::endl;
        
        
        if (status == GSL_SUCCESS){
            std::cout << "converged to a minimum!" << std::endl;
        }
        
    }
    while (status == GSL_CONTINUE && iter < iter_max);
    
    gsl_vector_free(min_params);
    gsl_vector_free(initial_step_size);
    gsl_multimin_fminimizer_free (s);
    
    
}

double E_DES::ComputeSSR(const std::vector<std::vector<double>> &param_est_data_set,
                  const std::vector<double> &glucoses,
                  const std::vector<double> &insulins){
    if (param_est_data_set.size() != glucoses.size()) {
        std::cout << "ERROR(E_DES::ComputeSSR): two vector sizes do NOT match!" << std::endl;
        return 0.;
    }
    double SSR = 0.;
    double glu_tmp = 0., ins_tmp = 0.;
    for (int i = 1; i < param_est_data_set.size(); ++i) { // skip the initial data point
        glu_tmp = pow( (glucoses[i] - param_est_data_set[i][1]) / param_est_data_set[i][2], 2.);
        ins_tmp = pow( (insulins[i] - param_est_data_set[i][3]) / param_est_data_set[i][4], 2.);
        SSR += glu_tmp + ins_tmp;
    }
    return SSR;
}


double E_DES::gsl_min_fitted_params_SSR_func (const gsl_vector *v, void *params){
    double SSR = 0.;
    std::pair<std::vector<double>, std::vector<double>> * paramsP = static_cast<std::pair<std::vector<double>, std::vector<double>> *>(params);
    std::vector<double> SS_tmp = paramsP->first;
    std::vector<double> SD_tmp = paramsP->second;
    
    // assign values to the varied params
    std::vector<double *> min_paramsP = {&k1, &k2, &k3, &k4, &k5, &k6, &k7, &k8, &k9, &k10, &k11, &k12, &sigma, &KM};
    for (int i = 0; i < min_paramsP.size(); ++i) {
        *min_paramsP[i] = pow( 10, SS_tmp[i] + SD_tmp[i] * tanh( gsl_vector_get(v, i)) );
    }
    
    for (int i = 0; i < input_param_sets.size(); ++i) {
        ClearPreRuns();
        // set up input_params for the current run
        SetInputParams(input_param_sets[i]);
        // prepare for evolution: set check_pts
        for (auto param_est_data_set_row: param_est_data_sets[i]) {
            check_pts.push_back(param_est_data_set_row[0]);
        }
        // set up initial conditions
        double Gpl_init_tmp = param_est_data_sets[i][0][1];
        double Ipl_init_tmp = param_est_data_sets[i][0][3];
        std::vector<double> init_conditions = {check_pts[0], 0., Gpl_init_tmp, Ipl_init_tmp, 0., 0.};
        SetInitConditions(init_conditions);
        // evolve
        Solver_gsl();
        // calculate SSR
        SSR += ComputeSSR(param_est_data_sets[i], glucoses, insulins);
    }
    return SSR;
}



// ++++++++++++++++++++++++++++++   Initialization of static memembers  ++++++++++++++++++++++++++++++

// subject type
int E_DES::type = 0;

// pre-setted fitted-params: healthy person
double E_DES::k1_H = 0.0183626;
double E_DES::k2_H = 0.0937043;
double E_DES::k3_H = 0.00428242;
double E_DES::k4_H = 0.000218343;
double E_DES::k5_H = 0.357796;
double E_DES::k6_H = 0.0577189;
double E_DES::k7_H = 0.248337;
double E_DES::k8_H = 6.70016;
double E_DES::k9_H = 0.;         // short-acting insulin
double E_DES::k10_H = 0.;        // short-acting insulin
double E_DES::k11_H = 0.013054;
double E_DES::k12_H = 0.669478;
double E_DES::sigma_H = 1.42996;
double E_DES::KM_H = 14.7465;

// pre-setted fitted-params: D1 (currently, same as D2)
double E_DES::k1_D1 = 0.0171219;
double E_DES::k2_D1 = 0.0682668;
double E_DES::k3_D1 = 0.00374951;
double E_DES::k4_D1 = 0.000510048;
double E_DES::k5_D1 = 0.211746;
double E_DES::k6_D1 = 0.0187956;
double E_DES::k7_D1 = 0.172134;
double E_DES::k8_D1 = 2.40894;
double E_DES::k9_D1 = 0.;         // short-acting insulin
double E_DES::k10_D1 = 0.;        // short-acting insulin
double E_DES::k11_D1 = 0.000422715;
double E_DES::k12_D1 = 0.0191808;
double E_DES::sigma_D1 = 1.52822;
double E_DES::KM_D1 = 18.937;

//// pre-setted fitted-params: D2
//double E_DES::k1_D2 = 0.0171219;
//double E_DES::k2_D2 = 0.0682668;
//double E_DES::k3_D2 = 0.00374951;
//double E_DES::k4_D2 = 0.000510048;
//double E_DES::k5_D2 = 0.211746;
//double E_DES::k6_D2 = 0.0187956;
//double E_DES::k7_D2 = 0.172134;
//double E_DES::k8_D2 = 2.40894;
//double E_DES::k9_D2 = 0.;         // short-acting insulin
//double E_DES::k10_D2 = 0.;       // short-acting insulin
//double E_DES::k11_D2 = 0.000422715;
//double E_DES::k12_D2 = 0.0191808;
//double E_DES::sigma_D2 = 1.52822;
//double E_DES::KM_D2 = 18.937;

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

// Manually tweaked fitted-params for D2, comparing with those in original E_DES paper (7/1/2018)
double E_DES::k1_D2 = 1.45E-2;
double E_DES::k2_D2 = 2.76E-1;
double E_DES::k3_D2 = 0.005; // tweaked params
double E_DES::k4_D2 = 2.35E-4;
double E_DES::k5_D2 = 0.01; // tweaked params
double E_DES::k6_D2 = 0.4; // tweaked params
double E_DES::k7_D2 = 0.3; // tweaked params
double E_DES::k8_D2 = 2.; // tweaked params
double E_DES::k9_D2 = 0.;         // short-acting insulin
double E_DES::k10_D2 = 0.;        // short-acting insulin
double E_DES::k11_D2 = 0.05; // tweaked params
double E_DES::k12_D2 = 2.84E-1;
double E_DES::sigma_D2 = 1.34;
double E_DES::KM_D2 = 16.65; // tweaked params

// input params
double E_DES::Dmeal = 75E3;
double E_DES::Mb = 75.;

// fitted-params: default healthy person
double E_DES::k1 = 0.0183626;
double E_DES::k2 = 0.0937043;
double E_DES::k3 = 0.00428242;
double E_DES::k4 = 0.000218343;
double E_DES::k5 = 0.357796;
double E_DES::k6 = 0.0577189;
double E_DES::k7 = 0.248337;
double E_DES::k8 = 6.70016;
double E_DES::k9 = 0.;         // short-acting insulin
double E_DES::k10 = 0.;        // short-acting insulin
double E_DES::k11 = 0.013054;
double E_DES::k12 = 0.669478;
double E_DES::sigma = 1.42996;
double E_DES::KM = 14.7465;

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

// time_offset
double E_DES::time_offset = 0.;

// initial_conditions
double E_DES::t_init = 0.;      // min
double E_DES::MGgut_init= 0.;  // mg
double E_DES::Gpl_init = 5.;   // mmol/L
double E_DES::Ipl_init = 8.;   // U/L
double E_DES::Jpl_init = 0.;   // mU/L/min
double E_DES::Iif_init = 0.;   // mU/L

// evolved parameters at the current time instant
double E_DES::t_curr = 0.;      // min
double E_DES::MGgut_curr = 0.;  // mg
double E_DES::Gpl_curr = 5.;   // mmol/L
double E_DES::Ipl_curr = 8.;   // U/L
double E_DES::Jpl_curr = 0.;   // mU/L/min
double E_DES::Iif_curr = 0.;   // mU/L

// default initialization for other parameters
std::vector<double> E_DES::check_pts;
std::vector<double> E_DES::glucoses;
std::vector<double> E_DES::E_DES::insulins;
std::vector<double> E_DES::time_instants;
std::vector<E_DES::data_set> E_DES::param_est_data_sets;
std::vector<std::vector<double>> E_DES::input_param_sets;
std::vector<double> E_DES::minParams_init;
std::vector<double> E_DES::minParams_fval_curr;


