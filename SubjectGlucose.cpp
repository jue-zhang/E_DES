//
//  SubjectGlucose.cpp
//  E_DES
//
//  Created by Jue on 7/11/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "SubjectGlucose.hpp"

void SubjectGlucose::SetSubjectSpecifics(const std::tuple<int, double, double, double> &specifics){
    type        = std::get<0>(specifics);
    bodyMass    = std::get<1>(specifics);
    fastGlucose = std::get<2>(specifics);
    fastInsulin = std::get<3>(specifics);
    SetFittedParams_EDES(type);
    SetFittedParams_EDES_Ex(type);
}

// ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++

void SubjectGlucose::ClearPreEvolutionOutput(){
    time_instants.clear();
    glucoses.clear();
    insulins.clear();
    // time_instants_interp
    if (time_instants_interp != nullptr) delete [] time_instants_interp;
    // glucose
    if (glucose_spline != nullptr && glucoses_interp != nullptr) {
        gsl_spline_free (glucose_spline);
        delete [] glucoses_interp;
    }
    if (glucose_acc != nullptr) {
        gsl_interp_accel_free (glucose_acc);
    }
    // insulin
    if (insulin_spline != nullptr && insulins_interp != nullptr) {
        gsl_spline_free (insulin_spline);
        delete [] insulins_interp;
    }
    if (insulin_acc != nullptr) {
        gsl_interp_accel_free (insulin_acc);
    }
}

std::pair<double, double> SubjectGlucose::FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseEvents &f_ex_ents){
    double latest_meal_tmp = 0., latest_eat_time_tmp = 0.;
    auto it = f_ex_ents.rbegin();
    while (it != f_ex_ents.rend()) {
        if (std::get<1>(*it) > 0.) {
            latest_meal_tmp = std::get<1>(*it);
            latest_eat_time_tmp = std::get<0>(*it);
            break;
        }
        ++it;
    }
    return {latest_meal_tmp, latest_eat_time_tmp};
}

void SubjectGlucose::CombineFoodIntakeExerciseEvents(){
    
    TypeDailyFoodIntakeExerciseEvents &ret = foodIntakeExerciseEvents;
    ret.clear();
    
    // split exerciseEvents into exercise info at particular time instants
    std::vector<std::tuple<double, double, int>> ex_time_instants;
    double t_F_prev = -100.; // take care of the case with two CONSECTIVE different exercises (careful treatment at the overlap)
    for (int i = 0; i < exerciseEvents.size(); ++i) {
        auto t_I = std::get<0>(exerciseEvents[i]);
        auto t_F = std::get<1>(exerciseEvents[i]);
        auto intensity = std::get<2>(exerciseEvents[i]);
        auto exType = std::get<3>(exerciseEvents[i]);
        if ( fabs(t_I - t_F_prev) < 1E-5 ) { // two consective different exericses, modify the previous element
            std::get<1>(ex_time_instants[ex_time_instants.size()-1]) = intensity;
            std::get<2>(ex_time_instants[ex_time_instants.size()-1]) = exType;
        }
        else
            ex_time_instants.push_back({t_I, intensity, exType});
        ex_time_instants.push_back({t_F, 0., 0});
        t_F_prev = t_F;
    }
    
    // merge foodIntake and ex_time_instants, according to the orders of their time instants
    int i_f = 0, i_e = 0;
    auto sz_f = foodIntakeEvents.size();
    auto sz_e = ex_time_instants.size();
    while (i_f < sz_f && i_e < sz_e) {
        auto t_f = foodIntakeEvents[i_f].first;
        auto t_e = std::get<0>(ex_time_instants[i_e]);
        auto intensity = std::get<1>(ex_time_instants[i_e]);
        auto exType = std::get<2>(ex_time_instants[i_e]);
        auto latest_food_intake = FindTheLatestFoodIntakeEvents(ret);// find out the latest food intake event
        double latest_meal_tmp = latest_food_intake.first;
        double latest_eat_time_tmp = latest_food_intake.second;
        if (t_f < t_e) {
            // check if the food intake is still in an exercise
            bool stillInEx = false;
            if (ret.size() > 0 && std::get<2>(ret[ret.size()-1]) > 0. ) stillInEx = true;
            if (stillInEx) {
                ret.push_back({t_f, foodIntakeEvents[i_f].second, latest_meal_tmp, t_f - latest_eat_time_tmp,
                    std::get<4>(ret[ret.size()-1]), std::get<5>(ret[ret.size()-1])});
            }
            else
                ret.push_back({t_f, foodIntakeEvents[i_f].second, latest_meal_tmp, t_f -latest_eat_time_tmp, 0., 0. });
            ++i_f;
        }
        if (t_f == t_e) {
            ret.push_back({t_f, foodIntakeEvents[i_f].second, latest_meal_tmp, t_f - latest_eat_time_tmp, intensity, exType });
            ++i_f;
            ++i_e;
        }
        if (t_f > t_e) {
            ret.push_back({t_e, 0, latest_meal_tmp, t_e - latest_eat_time_tmp, intensity, exType });
            ++i_e;
        }
    }
    while (i_f < sz_f) {
        auto latest_food_intake = FindTheLatestFoodIntakeEvents(ret);// find out the latest food intake event
        double latest_meal_tmp = latest_food_intake.first;
        double latest_eat_time_tmp = latest_food_intake.second;
        ret.push_back({foodIntakeEvents[i_f].first, foodIntakeEvents[i_f].second, latest_meal_tmp,
            foodIntakeEvents[i_f].first - latest_eat_time_tmp, 0., 0. });
        ++i_f;
    }
    while (i_e < sz_e) {
        auto latest_food_intake = FindTheLatestFoodIntakeEvents(ret);// find out the latest food intake event
        double latest_meal_tmp = latest_food_intake.first;
        double latest_eat_time_tmp = latest_food_intake.second;
        auto t_e = std::get<0>(ex_time_instants[i_e]);
        auto intensity = std::get<1>(ex_time_instants[i_e]);
        auto exType = std::get<2>(ex_time_instants[i_e]);
        ret.push_back({t_e, 0, latest_meal_tmp, t_e - latest_eat_time_tmp, intensity, exType });
        ++i_e;
    }
    
//    // check
//    for (int i = 0; i < ret.size(); ++i) {
//        std::cout << std::get<0>(ret[i])/60.  << " " << std::get<1>(ret[i]) << " "
//        << std::get<2>(ret[i]) << " " << std::get<3>(ret[i]) << " " << std::get<4>(ret[i])
//        << " " << std::get<5>(ret[i]) << std::endl;
//    }
    
}

void SubjectGlucose::EvolutionUnderFoodIntakeExerciseEvents(){
    
    CombineFoodIntakeExerciseEvents();
    
    // ----- EDES_EX ------
    // check if 'foodIntakeExerciseEvents' are set properly
    if (foodIntakeExerciseEvents.size() < 2) {
        std::cout << "ERROR(SubjectGlucose::EvolutionUnderFoodIntakeExerciseEvents): the size of 'foodIntakeExerciseEvents' must be larger than 1" << std::endl;
        return;
    }
    
    ClearPreEvolutionOutput();

    // Consideration of 'exerciseEvents'
    auto num_evols = foodIntakeExerciseEvents.size() - 1; // total number of evolutions
    std::vector<double> initialConditions_tmp;
    const int evol_steps = 80;
    for (int evol = 0; evol < num_evols ; ++evol) {
        // set up params for the current evolution
        // set evolution specifices
        auto this_meal  = std::get<1>(foodIntakeExerciseEvents[evol]);
        auto latest_meal    = std::get<2>(foodIntakeExerciseEvents[evol]);
        auto delta_meal         = std::get<3>(foodIntakeExerciseEvents[evol]);
        auto ex_intensity  = std::get<4>(foodIntakeExerciseEvents[evol]);
        auto ex_type    = std::get<5>(foodIntakeExerciseEvents[evol]);
        eDES_ex.SetEvolutionSpecifics({type, bodyMass, fastGlucose, fastInsulin, this_meal,
            latest_meal, delta_meal, ex_type});
        // initial conditions
        auto tI_tmp = std::get<0>(foodIntakeExerciseEvents[evol]);
        auto tF_tmp = std::get<0>(foodIntakeExerciseEvents[evol+1]);
        if (evol == 0) initialConditions_tmp = {tI_tmp, 0., fastGlucose, fastInsulin, 0., 0., 0., ex_intensity};
        else{
            initialConditions_tmp = eDES_ex.GetCurrentEvolvedParams();
            initialConditions_tmp[6] = initialConditions_tmp[7]; // in the new run, ExPre_init = (Ex from the previous run)
            initialConditions_tmp[7] += ex_intensity; // add new ex_intensity
        }
        eDES_ex.SetInitConditions(initialConditions_tmp);
        eDES_ex.SetCheckPts(tI_tmp, tF_tmp, evol_steps);
        // perform evolution
        eDES_ex.ClearPreRunOutput();
        eDES_ex.Solver_gsl();
        // store the evolved results (Note: avoid duplications on the boundary points between two consective evolutions)
        int overlap = 1;
        if (evol == num_evols - 1) overlap = 0;
        time_instants.insert(time_instants.end(), eDES_ex.time_instants.begin(), eDES_ex.time_instants.end()-overlap);
        glucoses.insert(glucoses.end(), eDES_ex.glucoses.begin(), eDES_ex.glucoses.end()-overlap);
        insulins.insert(insulins.end(), eDES_ex.insulins.begin(), eDES_ex.insulins.end()-overlap);
    }
    
//        // ----- EDES ------
//    // check if 'foodIntakeEvents' are set properly
//    if (foodIntakeEvents.size() < 2) {
//        std::cout << "ERROR(SubjectGlucose::EvolutionUnderFoodIntakeExerciseEvents): the input size of 'foodIntakeEvents' must be larger than 1" << std::endl;
//        return;
//    }
//    
//    ClearPreEvolutionOutput();
//    
//    // No consideration of 'exerciseEvents'
//    auto num_evols = foodIntakeEvents.size() - 1; // total number of evolutions
//    double foodIntake_tmp, exercise_tmp, tI_tmp, tF_tmp;
//    std::vector<double> initialConditions_tmp;
//    const int evol_steps = 40;
//    for (int evol = 0; evol < num_evols ; ++evol) {
//        // set up params for the current evolution
//        foodIntake_tmp = foodIntakeEvents[evol].second;
//        exercise_tmp = 0;
//        eDES.SetEvolutionSpecifics({type, bodyMass, foodIntake_tmp, exercise_tmp, fastGlucose, fastInsulin});
//        tI_tmp = foodIntakeEvents[evol].first;
//        tF_tmp = foodIntakeEvents[evol+1].first;
//        if (evol == 0) initialConditions_tmp = {tI_tmp, 0., fastGlucose, fastInsulin, 0., 0.};
//        else initialConditions_tmp = eDES.GetCurrentEvolvedParams();
//        eDES.SetInitConditions(initialConditions_tmp);
//        eDES.SetCheckPts(tI_tmp, tF_tmp, evol_steps);
//        // perform evolution
//        eDES.ClearPreRunOutput();
//        eDES.Solver_gsl();
//        // store the evolved results (Note: avoid duplications on the boundary points between two consective evolutions)
//        int overlap = 1;
//        if (evol == num_evols - 1) overlap = 0;
//        time_instants.insert(time_instants.end(), eDES.time_instants.begin(), eDES.time_instants.end()-overlap);
//        glucoses.insert(glucoses.end(), eDES.glucoses.begin(), eDES.glucoses.end()-overlap);
//        insulins.insert(insulins.end(), eDES.insulins.begin(), eDES.insulins.end()-overlap);
//    }
    
    // set up for interpolation: ref. https://www.gnu.org/software/gsl/doc/html/interp.html
    auto sz = time_instants.size();
    glucose_spline = gsl_spline_alloc (gsl_interp_cspline, sz);
    glucose_acc = gsl_interp_accel_alloc ();
    insulin_spline = gsl_spline_alloc (gsl_interp_cspline, sz);
    insulin_acc = gsl_interp_accel_alloc ();
    time_instants_interp = new double[sz];
    glucoses_interp = new double[sz];
    insulins_interp = new double[sz];
    for (int i = 0; i < sz; ++i) {
        time_instants_interp[i] = time_instants[i];
        glucoses_interp[i] = glucoses[i];
        insulins_interp[i] = insulins[i];
    }
    gsl_spline_init (glucose_spline, time_instants_interp, glucoses_interp, sz);
    gsl_spline_init (insulin_spline, time_instants_interp, insulins_interp, sz);
    
}

std::vector<std::pair<double, double>> SubjectGlucose::GlucoseUnderFoodIntakeExerciseEvents(double timeInterval){
    
    double t_init = std::get<0>(foodIntakeExerciseEvents[0]);
    double t_end = std::get<0>(foodIntakeExerciseEvents[foodIntakeExerciseEvents.size()-1]);
    double t_tmp = t_init, glucose_tmp = 0.;
    std::vector<std::pair<double, double>> ret;
    while (t_tmp < t_end) {
        glucose_tmp = gsl_spline_eval (glucose_spline, t_tmp, glucose_acc);
        ret.push_back({t_tmp, glucose_tmp});
        t_tmp += timeInterval;
    }
    // add the last time instant
    ret.push_back({t_end, gsl_spline_eval (glucose_spline, t_end, glucose_acc)});
    return ret;
}


std::vector<std::pair<double, double>> SubjectGlucose::InsulinUnderFoodIntakeExerciseEvents(double timeInterval){
    
    double t_init = std::get<0>(foodIntakeExerciseEvents[0]);
    double t_end = std::get<0>(foodIntakeExerciseEvents[foodIntakeExerciseEvents.size()-1]);
    double t_tmp = t_init, insulin_tmp = 0.;
    std::vector<std::pair<double, double>> ret;
    while (t_tmp < t_end) {
        insulin_tmp = gsl_spline_eval (insulin_spline, t_tmp, insulin_acc);
        ret.push_back({t_tmp, insulin_tmp});
        t_tmp += timeInterval;
    }
    // add the last time instant
    ret.push_back({t_end, gsl_spline_eval (insulin_spline, t_end, insulin_acc)});
    return ret;
}




std::vector<std::pair<double, double>> SubjectGlucose::GlucoseUnderFoodIntakeExerciseEvents(const std::vector<double> &check_pts){
    std::vector<std::pair<double, double>> ret;
    for (int i = 0; i < check_pts.size(); ++i) {
        auto glucose_tmp = gsl_spline_eval (glucose_spline, check_pts[i], glucose_acc);
        ret.push_back({check_pts[i], glucose_tmp});
    }
    return ret;
}


std::vector<std::pair<double, double>> SubjectGlucose::InsulinUnderFoodIntakeExerciseEvents(const std::vector<double> &check_pts){
    std::vector<std::pair<double, double>> ret;
    for (int i = 0; i < check_pts.size(); ++i) {
        auto insulin_tmp = gsl_spline_eval (insulin_spline, check_pts[i], insulin_acc);
        ret.push_back({check_pts[i], insulin_tmp});
    }
    return ret;
}



//std::vector<std::pair<double, double>> SubjectGlucose::GetGlucoses(){
//    std::vector<std::pair<double, double>> ret;
//    for (int i = 0; i < time_instants.size(); ++i) {
//        ret.push_back({time_instants[i], glucoses[i]});
//    }
//    return ret;
//}


// ++++++++++++++++++++++++++++++   Methods for estimating fitted-params  ++++++++++++++++++++++++++++++

void SubjectGlucose::SetDatasetsForParameterEstimation(const std::vector<std::string> &dailyDataSets_i,
                                              const std::vector<TypeSubjectSpecifics> &dailySubjectSpecifics_i,
                                              const std::vector<TypeDailyFoodIntakeEvents> &dailyFoodIntakeEvents_i,
                                              const std::vector<TypeDailyExerciseEvents> &dailyExerciseEvents_i){
    
    dailySubjectSpecifics = dailySubjectSpecifics_i;
    dailyFoodIntakeEvents = dailyFoodIntakeEvents_i;
    dailyExerciseEvents = dailyExerciseEvents_i;
    
    std::ifstream ifile_gi;
    TypeDataSet dailyDataSet_tmp;
    TypeDataSetRow dataSetRow_tmp;
    double ti, glu_i, glu_err_i, ins_i, ins_err_i;
    for (auto i = 0; i < dailyDataSets_i.size(); ++i) {
        dailyDataSet_tmp.clear();
        ifile_gi.open(dailyDataSets_i[i], std::ifstream::in);
        while (ifile_gi >> ti >> glu_i >> glu_err_i >> ins_i >> ins_err_i){
            dataSetRow_tmp = {ti, glu_i, glu_err_i, ins_i, ins_err_i};
            dailyDataSet_tmp.push_back(dataSetRow_tmp);
        }
        dailyDataSets.push_back(dailyDataSet_tmp);
        ifile_gi.close();
    }
    
}

void SubjectGlucose::ClearDatasetsForParameterEstimation(){
    dailyDataSets.clear();
    dailySubjectSpecifics.clear();
    dailyFoodIntakeEvents.clear();
    dailyExerciseEvents.clear();
}

std::map<std::string, double> SubjectGlucose::EstimateFittedParameters_EDES (const std::vector<std::string> &chosenFittedParams_str){
    
    // obtain the initial fitted-params before fitting
    int type = std::get<0>(dailySubjectSpecifics[0]);
    eDES.SetFittedParams(type);
    std::map<std::string, double> fittedParamsInit = eDES.GetFittedParamsDict();
    
    // set up the dictionary for the fitted params
    std::map<std::string, std::pair<double, double>> dict_lower_upper_bounds;
    // default setting for the allowed ranges of 'fp', [0.01 * currVal, 100 * currVal]
    double lower_tmp, upper_tmp;
    for (int i = 0; i < chosenFittedParams_str.size(); ++i) {
        lower_tmp = 0.01 * fittedParamsInit[chosenFittedParams_str[i]];
        upper_tmp = 100 * fittedParamsInit[chosenFittedParams_str[i]];
        dict_lower_upper_bounds.insert({ chosenFittedParams_str[i], {lower_tmp, upper_tmp} });
    }
    // special treatment on some fitted parameters
    // k1, k2, sigma, KM
    if (dict_lower_upper_bounds.find("k1") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k1"] = { 0.005, 0.035 };
    }
    if (dict_lower_upper_bounds.find("k2") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k2"] = { 0.05, 0.8 };
    }
    if (dict_lower_upper_bounds.find("sigma") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["sigma"] = { 1., 2. };
    }
    if (dict_lower_upper_bounds.find("KM") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["KM"] = { 5., 30. };
    }
    
    // set up the minimized params (gsl_vector) used in the gsl subroutine of minimization
    auto num_params = chosenFittedParams_str.size();
    gsl_vector *initial_step_size, *min_params; // 'min_params' -- params to be performed minimizations
    min_params = gsl_vector_alloc (num_params);
    initial_step_size = gsl_vector_alloc (num_params);
    // parameter transformation, as in GSL no constraints are imposed on 'min_params'
    //      1. transform fitted-params 'fp' into 'log10(fp)'
    //      2. for each 'log10(fp)', if it should be within [a, b], then 'min_params' is given
    //          log10(fp) = SS + SD * tanh(min_params), where SS = (a+b)/2, SD = (b-a)/2
    //          ref.: http://cafim.sssup.it/~giulio/software/multimin/multimin.html
    //    std::vector<double> SS, SD; // store the transformation info., used for later conversion from 'min_params' to 'fitted-paras'
    double SS_tmp, SD_tmp, log10_lower_tmp, log10_upper_tmp;
    int index = 0;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        gsl_vector_set (min_params, index, 0); // centering the initial values of 'min_params' at 0
        gsl_vector_set (initial_step_size, index, 0.01); // initial step size for 'min_params'
        ++index;
    }
    
    // set up the minimizer
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; // method used for minimization: Nelder-Mead
    gsl_multimin_fminimizer *s = nullptr; // initial the minimizer
    gsl_multimin_function min_func; // 'min_func' -- functions to be minimized over
    min_func.n = num_params;
    min_func.f = gsl_EstimateFittedParameters_EDES;
    // prepare the params to be passed to 'gsl_min_fitted_params_SSR_func'
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = {this, &dict_lower_upper_bounds};
    min_func.params = &params_tmp;
    
    s = gsl_multimin_fminimizer_alloc (T, num_params);
    gsl_multimin_fminimizer_set (s, &min_func, min_params, initial_step_size);
    
    
    // set up stopping criteria
    std::size_t iter_max = 1000; // max number of iterations
    std::size_t iter = 0;  // current iter
    int status;
    double SSR_curr = 0.;
    int num_data_pts = 0;
    for (int i = 0; i < dailyDataSets.size(); ++i) {
        num_data_pts += dailyDataSets[0].size() - 1;
    }
    std::cout << "num_data_pts = " << num_data_pts << std::endl;
    double epsabs = num_data_pts * 2 * 0.5; // average fitting within about 0.5 sigma for all data points, '2' for both glucoses and insulins
    
    // iterate:
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) // break when error occurs
            break;
        
        SSR_curr = s->fval;
        
        status = gsl_multimin_test_size (SSR_curr, epsabs);
        
        std::cout << iter << " ";
        int index = 0;
        for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
            auto fittedParamName = it->first;
            log10_lower_tmp = log10((it->second).first);
            log10_upper_tmp = log10((it->second).second);
            SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
            SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
            auto fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(s->x, index)) );
            std::cout << fittedParamName << " " << fitted_param_tmp << " ";
            ++index;
        }
        std::cout << SSR_curr << std::endl;
        
        if (status == GSL_SUCCESS){
            std::cout << "converged to a minimum!" << std::endl;
        }
        
        
    }
    while (status == GSL_CONTINUE && iter < iter_max);
    
    // output all fitted-parameters after fitting
    index = 0;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        auto fittedParamName = it->first;
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        auto fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(s->x, index)) );
        fittedParamsInit[fittedParamName] = fitted_param_tmp;
        ++index;
    }
    
    gsl_vector_free(min_params);
    gsl_vector_free(initial_step_size);
    gsl_multimin_fminimizer_free (s);
    
    return fittedParamsInit;
}

double SubjectGlucose::gsl_EstimateFittedParameters_EDES (const gsl_vector *v, void *paramsP){

    // passing parameters
    auto paramsP_tmp = *static_cast< std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>>*> *>(paramsP);
    auto this2 = std::get<0>(paramsP_tmp);
    auto dict_lower_upper_bounds = *std::get<1>(paramsP_tmp);
    
    // recover the fitted params from gsl_vector to their normal forms
    std::map<std::string, double> fittedParamsNew;
    int index = 0;
    std::string fittedParamName;
    double log10_lower_tmp, log10_upper_tmp, SS_tmp, SD_tmp, fitted_param_tmp;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        fittedParamName = it->first;
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(v, index)) );
        fittedParamsNew[fittedParamName] = fitted_param_tmp;
        ++index;
    }
    
    // compute SSR
    double SSR = 0.;
    SubjectGlucose sGlucose_tmp;
    // set up sGlucose_tmp: passing the daily datasets and other info
    sGlucose_tmp.ClearDatasetsForParameterEstimation();
    sGlucose_tmp.dailyDataSets = this2->dailyDataSets;
    sGlucose_tmp.dailySubjectSpecifics = this2->dailySubjectSpecifics;
    sGlucose_tmp.dailyFoodIntakeEvents = this2->dailyFoodIntakeEvents;
    sGlucose_tmp.dailyExerciseEvents = this2->dailyExerciseEvents;

    auto num_data_sets = sGlucose_tmp.dailyDataSets.size();
    for (int i = 0; i < num_data_sets; ++i) {
        sGlucose_tmp.eDES.ClearPreRunOutput();
        sGlucose_tmp.SetSubjectSpecifics(sGlucose_tmp.dailySubjectSpecifics[i]); // still default parameters according to 'type'
        sGlucose_tmp.SetFittedParams_EDES(fittedParamsNew); // update the minimized fitted-params
        sGlucose_tmp.SetFoodIntakeEvents(sGlucose_tmp.dailyFoodIntakeEvents[i]);
        sGlucose_tmp.SetExerciseEvents(sGlucose_tmp.dailyExerciseEvents[i]);
        sGlucose_tmp.EvolutionUnderFoodIntakeExerciseEvents();
        auto check_pts = sGlucose_tmp.ObtainCheckPtsFromDailyDataSets(sGlucose_tmp.dailyDataSets[i]);
        auto glucoses_tmp = sGlucose_tmp.GlucoseUnderFoodIntakeExerciseEvents(check_pts);
        auto insulins_tmp = sGlucose_tmp.InsulinUnderFoodIntakeExerciseEvents(check_pts);
        SSR += sGlucose_tmp.ComputeSSR(sGlucose_tmp.dailyDataSets[i], glucoses_tmp, insulins_tmp);
    }
    
    return SSR;
}


std::map<std::string, double> SubjectGlucose::EstimateFittedParameters_EDES_Ex (const std::vector<std::string> &chosenFittedParams_str){
    
    // obtain the initial fitted-params before fitting
    std::map<std::string, double> fittedParamsInit = eDES.GetFittedParamsDict();
    
    // set up the dictionary for the fitted params
    std::map<std::string, std::pair<double, double>> dict_lower_upper_bounds;
    // default setting for the allowed ranges of 'fp', [0.01 * currVal, 100 * currVal]
    double lower_tmp, upper_tmp;
    for (int i = 0; i < chosenFittedParams_str.size(); ++i) {
        lower_tmp = 0.01 * fittedParamsInit[chosenFittedParams_str[i]];
        upper_tmp = 100 * fittedParamsInit[chosenFittedParams_str[i]];
        dict_lower_upper_bounds.insert({ chosenFittedParams_str[i], {lower_tmp, upper_tmp} });
    }
    // special treatment on some fitted parameters
    // k1, k2, sigma, KM
    if (dict_lower_upper_bounds.find("k1") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k1"] = { 0.005, 0.035 };
    }
    if (dict_lower_upper_bounds.find("k2") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k2"] = { 0.05, 0.8 };
    }
    if (dict_lower_upper_bounds.find("sigma") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["sigma"] = { 1., 2. };
    }
    if (dict_lower_upper_bounds.find("KM") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["KM"] = { 5., 30. };
    }
    // k4e, k5e, k8e, lam:
    if (dict_lower_upper_bounds.find("k4e") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k4e"] = { 0.000049, 0.00005 };
    }
    if (dict_lower_upper_bounds.find("k5e") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k5e"] = { 0.000049, 0.00005 };
    }
    if (dict_lower_upper_bounds.find("k8e") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k8e"] = { 0.000049, 0.00005 };
    }
    if (dict_lower_upper_bounds.find("lam") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["lam"] = { 1/120., 1/119.9 };
    }
    
    // set up the minimized params (gsl_vector) used in the gsl subroutine of minimization
    auto num_params = chosenFittedParams_str.size();
    gsl_vector *initial_step_size, *min_params; // 'min_params' -- params to be performed minimizations
    min_params = gsl_vector_alloc (num_params);
    initial_step_size = gsl_vector_alloc (num_params);
    // parameter transformation, as in GSL no constraints are imposed on 'min_params'
    //      1. transform fitted-params 'fp' into 'log10(fp)'
    //      2. for each 'log10(fp)', if it should be within [a, b], then 'min_params' is given
    //          log10(fp) = SS + SD * tanh(min_params), where SS = (a+b)/2, SD = (b-a)/2
    //          ref.: http://cafim.sssup.it/~giulio/software/multimin/multimin.html
    //    std::vector<double> SS, SD; // store the transformation info., used for later conversion from 'min_params' to 'fitted-paras'
    double SS_tmp, SD_tmp, log10_lower_tmp, log10_upper_tmp;
    int index = 0;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        gsl_vector_set (min_params, index, 0); // centering the initial values of 'min_params' at 0
        gsl_vector_set (initial_step_size, index, 0.01); // initial step size for 'min_params'
        ++index;
    }
    
    // set up the minimizer
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; // method used for minimization: Nelder-Mead
    gsl_multimin_fminimizer *s = nullptr; // initial the minimizer
    gsl_multimin_function min_func; // 'min_func' -- functions to be minimized over
    min_func.n = num_params;
    min_func.f = gsl_EstimateFittedParameters_EDES_Ex;
    // prepare the params to be passed to 'gsl_min_fitted_params_SSR_func'
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = {this, &dict_lower_upper_bounds};
    min_func.params = &params_tmp;
    
    s = gsl_multimin_fminimizer_alloc (T, num_params);
    gsl_multimin_fminimizer_set (s, &min_func, min_params, initial_step_size);
    
    
    // set up stopping criteria
    std::size_t iter_max = 500; // max number of iterations
    std::size_t iter = 0;  // current iter
    int status;
    double SSR_curr = 0.;
    int num_data_pts = 0;
    for (int i = 0; i < dailyDataSets.size(); ++i) {
        num_data_pts += dailyDataSets[0].size() - 1;
    }
    double epsabs = num_data_pts * 2 * 0.5; // average fitting within about 0.5 sigma for all data points, '2' for both glucoses and insulins
    
    // iterate:
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) // break when error occurs
            break;
        
        SSR_curr = s->fval;
        
        status = gsl_multimin_test_size (SSR_curr, epsabs);
        
        std::cout << iter << " ";
        int index = 0;
        for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
            auto fittedParamName = it->first;
            log10_lower_tmp = log10((it->second).first);
            log10_upper_tmp = log10((it->second).second);
            SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
            SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
            auto fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(s->x, index)) );
            std::cout << fittedParamName << " " << fitted_param_tmp << " ";
            ++index;
        }
        std::cout << SSR_curr << std::endl;
        
        if (status == GSL_SUCCESS){
            std::cout << "converged to a minimum!" << std::endl;
        }
        
        
    }
    while (status == GSL_CONTINUE && iter < iter_max);
    
    // output all fitted-parameters after fitting
    index = 0;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        auto fittedParamName = it->first;
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        auto fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(s->x, index)) );
        fittedParamsInit[fittedParamName] = fitted_param_tmp;
        ++index;
    }
    
    gsl_vector_free(min_params);
    gsl_vector_free(initial_step_size);
    gsl_multimin_fminimizer_free (s);
    
    return fittedParamsInit;
}


double SubjectGlucose::gsl_EstimateFittedParameters_EDES_Ex (const gsl_vector *v, void *paramsP){
    
    // passing parameters
    auto paramsP_tmp = *static_cast< std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>>*> *>(paramsP);
    auto this2 = std::get<0>(paramsP_tmp);
    auto dict_lower_upper_bounds = *std::get<1>(paramsP_tmp);
    
    // recover the fitted params from gsl_vector to their normal forms
    std::map<std::string, double> fittedParamsNew;
    int index = 0;
    std::string fittedParamName;
    double log10_lower_tmp, log10_upper_tmp, SS_tmp, SD_tmp, fitted_param_tmp;
    for (auto it = dict_lower_upper_bounds.begin(); it != dict_lower_upper_bounds.end(); ++it) {
        fittedParamName = it->first;
        log10_lower_tmp = log10((it->second).first);
        log10_upper_tmp = log10((it->second).second);
        SS_tmp = ( log10_lower_tmp + log10_upper_tmp) / 2.;
        SD_tmp = ( - log10_lower_tmp + log10_upper_tmp) / 2.;
        fitted_param_tmp = pow( 10, SS_tmp + SD_tmp * tanh( gsl_vector_get(v, index)) );
        fittedParamsNew[fittedParamName] = fitted_param_tmp;
        ++index;
    }
    
    // compute SSR
    double SSR = 0.;
    SubjectGlucose sGlucose_tmp;
    // set up sGlucose_tmp: passing the daily datasets and other info
    sGlucose_tmp.ClearDatasetsForParameterEstimation();
    sGlucose_tmp.dailyDataSets = this2->dailyDataSets;
    sGlucose_tmp.dailySubjectSpecifics = this2->dailySubjectSpecifics;
    sGlucose_tmp.dailyFoodIntakeEvents = this2->dailyFoodIntakeEvents;
    sGlucose_tmp.dailyExerciseEvents = this2->dailyExerciseEvents;
    
    auto num_data_sets = sGlucose_tmp.dailyDataSets.size();
    for (int i = 0; i < num_data_sets; ++i) {
        sGlucose_tmp.eDES_ex.ClearPreRunOutput();
        sGlucose_tmp.SetSubjectSpecifics(sGlucose_tmp.dailySubjectSpecifics[i]); // still default parameters according to 'type'
        sGlucose_tmp.SetFittedParams_EDES_Ex(fittedParamsNew); // update the minimized fitted-params
        sGlucose_tmp.SetFoodIntakeEvents(sGlucose_tmp.dailyFoodIntakeEvents[i]);
        sGlucose_tmp.SetExerciseEvents(sGlucose_tmp.dailyExerciseEvents[i]);
        sGlucose_tmp.CombineFoodIntakeExerciseEvents();
        sGlucose_tmp.EvolutionUnderFoodIntakeExerciseEvents();
        auto check_pts = sGlucose_tmp.ObtainCheckPtsFromDailyDataSets(sGlucose_tmp.dailyDataSets[i]);
        auto glucoses_tmp = sGlucose_tmp.GlucoseUnderFoodIntakeExerciseEvents(check_pts);
        auto insulins_tmp = sGlucose_tmp.InsulinUnderFoodIntakeExerciseEvents(check_pts);
        SSR += sGlucose_tmp.ComputeSSR(sGlucose_tmp.dailyDataSets[i], glucoses_tmp, insulins_tmp);
    }
    
    return SSR;
}


double SubjectGlucose::ComputeSSR(const TypeDataSet &dailyDataSet,
                                  std::vector<std::pair<double, double>> &glucoses,
                                  std::vector<std::pair<double, double>> &insulins){
    if (dailyDataSet.size() != glucoses.size()) {
        std::cout << "ERROR(E_DES::ComputeSSR): two vector sizes do NOT match!" << std::endl;
        return 0.;
    }
    double SSR = 0.;
    double glu, glu_err, ins, ins_err, glu_evol, ins_evol;
    double glu_SSR = 0., ins_SSR = 0.;
    for (int i = 1; i < dailyDataSet.size(); ++i) { // skip the initial data point
        glu     = std::get<1>(dailyDataSet[i]);
        glu_err = std::get<2>(dailyDataSet[i]);
        ins     = std::get<3>(dailyDataSet[i]);
        ins_err = std::get<4>(dailyDataSet[i]);
        glu_evol = std::get<1>(glucoses[i]);
        ins_evol = std::get<1>(insulins[i]);
        if (fabs(glu) > 1E-5) // if there exists a measurement on glucose
            glu_SSR = pow( (glu_evol - glu) / glu_err, 2.);
        if (fabs(ins) > 1E-5) // if there exists a measurement on insulin
            ins_SSR = pow( (ins_evol - ins) / ins_err, 2.);
        SSR += glu_SSR + ins_SSR;
    }
    return SSR;
}


std::vector<double> SubjectGlucose::ObtainCheckPtsFromDailyDataSets( const TypeDataSet &dataSet ){
    std::vector<double> ret;
    for (int i = 0; i < dataSet.size(); ++i) {
        ret.push_back(std::get<0>(dataSet[i]));
    }
    return ret;
}

// ++++++++++++++++++++++++++++++   Methods for class itself  ++++++++++++++++++++++++++++++

SubjectGlucose::~SubjectGlucose(){
    // time_instants_interp
    if (time_instants_interp != nullptr) delete [] time_instants_interp;
    // glucose
    if (glucose_spline != nullptr && glucoses_interp != nullptr) {
        gsl_spline_free (glucose_spline);
        delete [] glucoses_interp;
    }
    if (glucose_acc != nullptr) {
        gsl_interp_accel_free (glucose_acc);
    }
    // insulin
    if (insulin_spline != nullptr && insulins_interp != nullptr) {
        gsl_spline_free (insulin_spline);
        delete [] insulins_interp;
    }
    if (insulin_acc != nullptr) {
        gsl_interp_accel_free (insulin_acc);
    }
    
}

