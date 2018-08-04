//
//  SubjectGlucose.cpp
//  E_DES
//
//  Created by Jue on 7/11/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include "SubjectGlucose.hpp"

//void SubjectGlucose::SetSubjectSpecifics(const std::tuple<int, double, double, double> &specifics){
//    type        = std::get<0>(specifics);
//    bodyMass    = std::get<1>(specifics);
//    fastGlucose = std::get<2>(specifics);
//    fastInsulin = std::get<3>(specifics);
//    SetFittedParams_EDES(type);
//    SetFittedParams_EDES_Ex(type);
//}

void SubjectGlucose::SetSubjectSpecifics(const TypeSubjectSpecifics &specifics){
    type        = std::get<0>(specifics);
    bodyMass    = std::get<1>(specifics);
    fastGlucose = std::get<2>(specifics);
    fastInsulin = std::get<3>(specifics);
    Gpl_basal   = std::get<4>(specifics);
    Ipl_basal   = std::get<5>(specifics);
    SetFittedParams_EDES(type);
    SetFittedParams_EDES_Ex(type);
    SetFittedParams_EDES_Ex_Med(type);
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
    return std::make_pair(latest_meal_tmp, latest_eat_time_tmp);
}

std::pair<double, double> SubjectGlucose::FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseSaIEvents &f_ex_saI_ents){
    double latest_meal_tmp = 0., latest_eat_time_tmp = 0.;
    auto it = f_ex_saI_ents.rbegin();
    while (it != f_ex_saI_ents.rend()) {
        if (std::get<1>(*it) > 0.) {
            latest_meal_tmp = std::get<1>(*it);
            latest_eat_time_tmp = std::get<0>(*it);
            break;
        }
        ++it;
    }
    return std::make_pair(latest_meal_tmp, latest_eat_time_tmp);
}

std::pair<double, double> SubjectGlucose::FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseSaILaIEvents &f_ex_saI_laI_ents){
    double latest_meal_tmp = 0., latest_eat_time_tmp = 0.;
    auto it = f_ex_saI_laI_ents.rbegin();
    while (it != f_ex_saI_laI_ents.rend()) {
        if (std::get<1>(*it) > 0.) {
            latest_meal_tmp = std::get<1>(*it);
            latest_eat_time_tmp = std::get<0>(*it);
            break;
        }
        ++it;
    }
    return std::make_pair(latest_meal_tmp, latest_eat_time_tmp);
}

std::tuple<int, double, double> SubjectGlucose::FindTheLatestLaIEvents(const TypeDailyFoodIntakeExerciseSaILaIEvents &f_ex_saI_laI_ents){
    int latest_laI_type_tmp = 0;
    double latest_laI_tmp = 0., latest_laI_time_tmp = 0.;
    auto it = f_ex_saI_laI_ents.rbegin();
    while (it != f_ex_saI_laI_ents.rend()) {
        if (std::get<9>(*it) > 0.) {
            latest_laI_type_tmp = std::get<8>(*it);
            latest_laI_tmp = std::get<9>(*it);
            latest_laI_time_tmp = std::get<0>(*it);
            break;
        }
        ++it;
    }
    return std::make_tuple(latest_laI_type_tmp, latest_laI_tmp, latest_laI_time_tmp);
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
            ex_time_instants.push_back(std::make_tuple(t_I, intensity, exType));
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
//    std::cout << "Combine food and Ex events:" << std::endl;
//    for (int i = 0; i < ret.size(); ++i) {
//        std::cout << std::get<0>(ret[i])  << " " << std::get<1>(ret[i]) << " "
//        << std::get<2>(ret[i]) << " " << std::get<3>(ret[i]) << " " << std::get<4>(ret[i])
//        << " " << std::get<5>(ret[i]) << std::endl;
//    }
    
}


void SubjectGlucose::CombineFoodIntakeExerciseSaIEvents(){
    
    CombineFoodIntakeExerciseEvents();
    
    TypeDailyFoodIntakeExerciseSaIEvents &ret = foodIntakeExerciseSaIEvents;
    ret.clear();
    
    // combine foodEx with saI
    auto sz_fe = foodIntakeExerciseEvents.size();
    auto sz_saI = SAInsulinEvents.size();
    std::size_t i_fe = 0, i_saI = 0;
    while (i_fe < sz_fe && i_saI < sz_saI) {
        auto t_fe = std::get<0>(foodIntakeExerciseEvents[i_fe]);
        auto t_saI = std::get<0>(SAInsulinEvents[i_saI]);
        if (t_fe < t_saI) { // foodEx event first, insert foodEx event
            auto saITuple = std::make_tuple(0, 0.);
            auto foodExMedTuple = std::tuple_cat(foodIntakeExerciseEvents[i_fe], saITuple);
            ret.push_back(foodExMedTuple);
            ++i_fe;
        }
        if ( fabs(t_fe - t_saI) < 1E-5 ) { // foodEx == saI, combine both
            auto sa_insulin_type = std::get<1>(SAInsulinEvents[i_saI]);
            auto sa_insulin = std::get<2>(SAInsulinEvents[i_saI]);
            auto saITuple = std::make_tuple(sa_insulin_type, sa_insulin);
            auto foodExMedTuple = std::tuple_cat(foodIntakeExerciseEvents[i_fe], saITuple);
            ret.push_back(foodExMedTuple);
            ++i_fe;
            ++i_saI;
        }
        if (t_fe > t_saI) { //foodEx < saI, insert saI
            auto sa_insulin_type = std::get<1>(SAInsulinEvents[i_saI]);
            auto sa_insulin = std::get<2>(SAInsulinEvents[i_saI]);
            auto latest_meal = FindTheLatestFoodIntakeEvents(ret);
            auto foodExMedTuple = std::make_tuple(t_saI, 0., latest_meal.first, t_saI - latest_meal.second, 0., 0, // foodEx
                                                  sa_insulin_type, sa_insulin);
            ret.push_back(foodExMedTuple);
            ++i_saI;
        }
    }
    while (i_fe < sz_fe) { // additional foodEx events
        auto medTuple = std::make_tuple(0, 0.);
        auto foodExMedTuple = std::tuple_cat(foodIntakeExerciseEvents[i_fe], medTuple);
        ret.push_back(foodExMedTuple);
        ++i_fe;
    }
    while (i_saI < sz_saI){ // additional saI
        auto t_saI = std::get<0>(SAInsulinEvents[i_saI]);
        auto sa_insulin_type = std::get<1>(SAInsulinEvents[i_saI]);
        auto sa_insulin = std::get<2>(SAInsulinEvents[i_saI]);
        auto latest_meal = FindTheLatestFoodIntakeEvents(ret);
        auto foodExMedTuple = std::make_tuple(t_saI, 0., latest_meal.first, t_saI - latest_meal.second, 0., 0, // foodEx
                                              sa_insulin_type, sa_insulin);
        ret.push_back(foodExMedTuple);
        ++i_saI;
    }
    
//    // check
//    std::cout << "Combine food, Ex, SaI events: " << std::endl;
//    for (int i = 0; i < ret.size(); ++i) {
//        std::cout << std::get<0>(ret[i])  << " " << std::get<1>(ret[i]) << " "
//        << std::get<2>(ret[i]) << " " << std::get<3>(ret[i]) << " " << std::get<4>(ret[i])
//        << " " << std::get<5>(ret[i]) << " " << std::get<6>(ret[i]) << " " << std::get<7>(ret[i]) << std::endl;
//    }
    
}

void SubjectGlucose::CombineFoodIntakeExerciseSaILaIEvents(){
    
    CombineFoodIntakeExerciseSaIEvents();
    
    TypeDailyFoodIntakeExerciseSaILaIEvents &ret = foodIntakeExerciseSaILaIEvents;
    ret.clear();
    
    auto sz_fes = foodIntakeExerciseSaIEvents.size();
    auto sz_l = LAInsulinEvents.size();
    std::size_t i_fes = 0, i_l = 0;
    while (i_fes < sz_fes && i_l < sz_l) {
        auto t_fes = std::get<0>(foodIntakeExerciseSaIEvents[i_fes]);
        auto t_l = std::get<0>(LAInsulinEvents[i_l]);
        auto latest_laI = FindTheLatestLaIEvents(ret);
        auto latest_meal = FindTheLatestFoodIntakeEvents(ret);
        if (t_fes < t_l) { // insert t_fes
            auto laI_tuple = std::make_tuple(0, 0., std::get<0>(latest_laI), std::get<1>(latest_laI), t_fes - std::get<2>(latest_laI));
            auto tmp = std::tuple_cat(foodIntakeExerciseSaIEvents[i_fes], laI_tuple);
            ret.push_back(tmp);
            ++i_fes;
        }
        if (t_fes == t_l) { // t_fes = t_l
            auto laI_curr_type = std::get<1>(LAInsulinEvents[i_l]);
            auto laI_curr = std::get<2>(LAInsulinEvents[i_l]);
            auto laI_tuple = std::make_tuple(laI_curr_type, laI_curr, std::get<0>(latest_laI), std::get<1>(latest_laI), t_fes - std::get<2>(latest_laI));
            auto tmp = std::tuple_cat(foodIntakeExerciseSaIEvents[i_fes], laI_tuple);
            ret.push_back(tmp);
            ++i_fes;
            ++i_l;
        }
        if (t_fes > t_l) { // insert t_l
            auto foodSaITuple = std::make_tuple(t_l, 0., latest_meal.first, t_l - latest_meal.second, 0., 0, 0, 0.);
            auto laI_curr_type = std::get<1>(LAInsulinEvents[i_l]);
            auto laI_curr = std::get<2>(LAInsulinEvents[i_l]);
            auto laI_tuple = std::make_tuple(laI_curr_type, laI_curr, std::get<0>(latest_laI), std::get<1>(latest_laI),
                                             t_l - std::get<2>(latest_laI));
            auto tmp = std::tuple_cat(foodSaITuple, laI_tuple);
            ret.push_back(tmp);
            ++i_l;
        }
    }
    while (i_fes < sz_fes) { // insert t_fes
        auto latest_laI = FindTheLatestLaIEvents(ret);
        auto t_fes = std::get<0>(foodIntakeExerciseSaIEvents[i_fes]);
        auto laI_tuple = std::make_tuple(0, 0., std::get<0>(latest_laI), std::get<1>(latest_laI), t_fes - std::get<2>(latest_laI));
        auto tmp = std::tuple_cat(foodIntakeExerciseSaIEvents[i_fes], laI_tuple);
        ret.push_back(tmp);
        ++i_fes;
    }
    while (i_l < sz_l) { // insert t_l
        auto t_l = std::get<0>(foodIntakeExerciseSaIEvents[i_l]);
        auto latest_meal = FindTheLatestFoodIntakeEvents(ret);
        auto latest_laI = FindTheLatestLaIEvents(ret);
        auto foodSaITuple = std::make_tuple(t_l, 0., latest_meal.first, t_l - latest_meal.second, 0., 0, 0, 0.);
        auto laI_curr_type = std::get<1>(LAInsulinEvents[i_l]);
        auto laI_curr = std::get<2>(LAInsulinEvents[i_l]);
        auto laI_tuple = std::make_tuple(laI_curr_type, laI_curr, std::get<0>(latest_laI), std::get<1>(latest_laI),
                                         t_l - std::get<2>(latest_laI));
        auto tmp = std::tuple_cat(foodSaITuple, laI_tuple);
        ret.push_back(tmp);
        ++i_l;
    }
    
//    // check
//    std::cout << "Combine food, Ex, SaI, LaI events: " << std::endl;
//    for (int i = 0; i < ret.size(); ++i) {
//        std::cout << std::get<0>(ret[i])  << " " << std::get<1>(ret[i]) << " "
//        << std::get<2>(ret[i]) << " " << std::get<3>(ret[i]) << " " << std::get<4>(ret[i])
//        << " " << std::get<5>(ret[i]) << " " << std::get<6>(ret[i]) << " " << std::get<7>(ret[i])
//        << " " << std::get<8>(ret[i]) << " " << std::get<9>(ret[i])
//        << " " << std::get<10>(ret[i]) << " " << std::get<11>(ret[i]) << " " << std::get<12>(ret[i])  << std::endl;
//    }
}

void SubjectGlucose::EvolutionUnderFoodIntakeExerciseEvents(){
    
    // ----- EDES_EX ------
    
    CombineFoodIntakeExerciseEvents();
    
    // check if 'foodIntakeExerciseEvents' are set properly
    if (foodIntakeExerciseEvents.size() < 2) {
        std::cout << "ERROR(SubjectGlucose::EvolutionUnderFoodIntakeExerciseEvents): the size of 'foodIntakeExerciseEvents' must be larger than 1" << std::endl;
        return;
    }
    
    ClearPreEvolutionOutput();

    // Consideration of 'exerciseEvents'
    auto num_evols = foodIntakeExerciseEvents.size() - 1; // total number of evolutions
    std::vector<double> initialConditions_tmp;
    const int evol_steps = 50;
    for (int evol = 0; evol < num_evols ; ++evol) {
        // set up params for the current evolution
        // set evolution specifices
        auto this_meal  = std::get<1>(foodIntakeExerciseEvents[evol]);
        auto latest_meal    = std::get<2>(foodIntakeExerciseEvents[evol]);
        auto delta_meal         = std::get<3>(foodIntakeExerciseEvents[evol]);
        auto ex_intensity  = std::get<4>(foodIntakeExerciseEvents[evol]);
        auto ex_type    = std::get<5>(foodIntakeExerciseEvents[evol]);
        eDES_ex.SetEvolutionSpecifics(std::make_tuple(type, bodyMass, fastGlucose, fastInsulin, this_meal,
                                                  latest_meal, delta_meal, ex_type));
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

void SubjectGlucose::EvolutionUnderFoodExInsulinEvents(){
    
    // ----- EDES_EX_MED ------
    
    CombineFoodIntakeExerciseSaILaIEvents();
    
    // check if 'foodIntakeExerciseEvents' are set properly
    if (foodIntakeExerciseSaILaIEvents.size() < 2) {
        std::cout << "ERROR(SubjectGlucose::EvolutionUnderFoodExInsulinEvents): the size of 'foodIntakeExerciseSaILaIEvents' must be larger than 1" << std::endl;
        return;
    }
    
    ClearPreEvolutionOutput();
    //
    //    <t, this_meal, latest_meal, delta_eat, exIntensity, exType, insulin_sa_type, insulin_sa, insulin_la_type_curr, insulin_la_curr, insulin_la_type_latest, insulin_la_latest, delta_insulin_la>
    
    // Consideration of 'foodIntakeExerciseEvents'
    auto num_evols = foodIntakeExerciseSaILaIEvents.size() - 1; // total number of evolutions
    std::vector<double> initialConditions_tmp;
    const int evol_steps = 40;
    for (int evol = 0; evol < num_evols ; ++evol) {
        // set up params for the current evolution
        // set evolution specifices
        int ex_type, insulin_sa_type, insulin_la_type_curr, insulin_la_type_latest, oralMed_type;
        double t_start, exIntensity, this_meal, latest_meal, delta_eat, insulin_sa, insulin_la_curr, insulin_la_latest, delta_insulin, oralMed;
        std::tie (t_start, this_meal, latest_meal, delta_eat, exIntensity, ex_type, insulin_sa_type, insulin_sa,
                  insulin_la_type_curr, insulin_la_curr, insulin_la_type_latest, insulin_la_latest, delta_insulin) = foodIntakeExerciseSaILaIEvents[evol];
        
        //        eDES_ex_med.SetEvolutionSpecifics(std::make_tuple(type, bodyMass, Gpl_basal, Ipl_basal,
        //                                                          this_meal, latest_meal, delta_eat, ex_type,
        //                                                          insulin_sa_type,
        //                                                          insulin_la_type_curr, insulin_la_type_latest,
        //                                                          insulin_la_curr, insulin_la_latest, delta_insulin,
        //                                                          oralMed_type, oralMed  ));
        eDES_ex_med.SetEvolutionSpecifics(std::make_tuple(type, bodyMass, fastGlucose, fastInsulin,
                                                          this_meal, latest_meal, delta_eat, ex_type,
                                                          insulin_sa_type,
                                                          insulin_la_type_curr, insulin_la_type_latest,
                                                          insulin_la_curr, insulin_la_latest, delta_insulin,
                                                          oralMed_type, oralMed  ));
        // initial conditions
        auto tI_tmp = std::get<0>(foodIntakeExerciseSaILaIEvents[evol]);
        auto tF_tmp = std::get<0>(foodIntakeExerciseSaILaIEvents[evol+1]);
        
        if (evol == 0) initialConditions_tmp = {tI_tmp, 0., fastGlucose, fastInsulin, fastGlucose-Gpl_basal, 0., exIntensity, insulin_sa, 0.};
        else{
            initialConditions_tmp = eDES_ex_med.GetCurrentEvolvedParams();
            initialConditions_tmp[5] = initialConditions_tmp[6]; // in the new run, ExPre_init = (Ex from the previous run)
            initialConditions_tmp[6] += exIntensity; // add new ex_intensity
            initialConditions_tmp[7] += insulin_sa; // add new insulin_sa
        }
        eDES_ex_med.SetInitConditions(initialConditions_tmp);
        eDES_ex_med.SetCheckPts(tI_tmp, tF_tmp, evol_steps);
        // perform evolution
        eDES_ex_med.ClearPreRunOutput();
        eDES_ex_med.Solver_gsl();
        // store the evolved results (Note: avoid duplications on the boundary points between two consective evolutions)
        int overlap = 1;
        if (evol == num_evols - 1) overlap = 0;
        time_instants.insert(time_instants.end(), eDES_ex_med.time_instants.begin(), eDES_ex_med.time_instants.end()-overlap);
        glucoses.insert(glucoses.end(), eDES_ex_med.glucoses.begin(), eDES_ex_med.glucoses.end()-overlap);
        insulins.insert(insulins.end(), eDES_ex_med.insulins.begin(), eDES_ex_med.insulins.end()-overlap);
    }
    
    
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


std::vector<std::pair<double, double>> SubjectGlucose::ObtainGlucose(double timeInterval){
    
    double t_init = std::get<0>(foodIntakeExerciseSaILaIEvents[0]);
    double t_end = std::get<0>(foodIntakeExerciseSaILaIEvents[foodIntakeExerciseSaILaIEvents.size()-1]);
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


std::vector<std::pair<double, double>> SubjectGlucose::ObtainInsulin(double timeInterval){
    
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




std::vector<std::pair<double, double>> SubjectGlucose::ObtainGlucose(const std::vector<double> &check_pts){
    std::vector<std::pair<double, double>> ret;
    for (int i = 0; i < check_pts.size(); ++i) {
        auto glucose_tmp = gsl_spline_eval (glucose_spline, check_pts[i], glucose_acc);
        ret.push_back({check_pts[i], glucose_tmp});
    }
    return ret;
}


std::vector<std::pair<double, double>> SubjectGlucose::ObtainInsulin(const std::vector<double> &check_pts){
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

void SubjectGlucose::SetDatasetsForParameterEstimation(const std::vector<std::string> &dailyDataSets_i,
                                                       const std::vector<TypeSubjectSpecifics> &dailySubjectSpecifics_i,
                                                       const std::vector<TypeDailyFoodIntakeEvents> &dailyFoodIntakeEvents_i,
                                                       const std::vector<TypeDailyExerciseEvents> &dailyExerciseEvents_i,
                                                       const std::vector<TypeDailySAInsulinEvents> &dailySAInsulinEvents_i,
                                                       const std::vector<TypeDailyLAInsulinEvents> &dailyLAInsulinEvents_i){
    
    dailySubjectSpecifics = dailySubjectSpecifics_i;
    dailyFoodIntakeEvents = dailyFoodIntakeEvents_i;
    dailyExerciseEvents = dailyExerciseEvents_i;
    dailySAInsulinEvents = dailySAInsulinEvents_i;
    dailyLAInsulinEvents = dailyLAInsulinEvents_i;
    
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
    dailySAInsulinEvents.clear();
    dailyLAInsulinEvents.clear();
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
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = std::make_tuple(this, &dict_lower_upper_bounds);
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
        auto glucoses_tmp = sGlucose_tmp.ObtainGlucose(check_pts);
        auto insulins_tmp = sGlucose_tmp.ObtainInsulin(check_pts);
        SSR += sGlucose_tmp.ComputeSSR(sGlucose_tmp.dailyDataSets[i], glucoses_tmp, insulins_tmp);
    }
    
    return SSR;
}


std::map<std::string, double> SubjectGlucose::EstimateFittedParameters_EDES_Ex (const std::vector<std::string> &chosenFittedParams_str){
    
    // obtain the initial fitted-params before fitting
    std::map<std::string, double> fittedParamsInit = eDES_ex.GetFittedParamsDict();
    
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
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = std::make_tuple(this, &dict_lower_upper_bounds);
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
        auto glucoses_tmp = sGlucose_tmp.ObtainGlucose(check_pts);
        auto insulins_tmp = sGlucose_tmp.ObtainInsulin(check_pts);
        SSR += sGlucose_tmp.ComputeSSR(sGlucose_tmp.dailyDataSets[i], glucoses_tmp, insulins_tmp);
    }
    
    return SSR;
}


std::map<std::string, double> SubjectGlucose::EstimateFittedParameters_EDES_Ex_Med (const std::vector<std::string> &chosenFittedParams_str){
    
    // obtain the initial fitted-params before fitting
    std::map<std::string, double> fittedParamsInit = eDES_ex_med.GetFittedParamsDict();
    
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
    // k5, k6, k7, k8, KM, from healthy to D2
//    if (dict_lower_upper_bounds.find("k5") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k5"] = { 0.1*fittedParamsInit["k5"], fittedParamsInit["k5"] };
//    }
//    if (dict_lower_upper_bounds.find("k6") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k6"] = { 0.01*fittedParamsInit["k6"], fittedParamsInit["k6"] };
//    }
//    if (dict_lower_upper_bounds.find("k7") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k7"] = { 0.01*fittedParamsInit["k7"], fittedParamsInit["k7"] };
//    }
//    if (dict_lower_upper_bounds.find("k8") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k8"] = { 0.01*fittedParamsInit["k8"], fittedParamsInit["k8"] };
//    }
//    if (dict_lower_upper_bounds.find("KM") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["KM"] = { fittedParamsInit["KM"], 30. };
//    }
    //    // k4e, k5e, k8e, lam:
    //    if (dict_lower_upper_bounds.find("k4e") != dict_lower_upper_bounds.end()) {
    //        dict_lower_upper_bounds["k4e"] = { 0.000049, 0.00005 };
    //    }
    //    if (dict_lower_upper_bounds.find("k5e") != dict_lower_upper_bounds.end()) {
    //        dict_lower_upper_bounds["k5e"] = { 0.000049, 0.00005 };
    //    }
    //    if (dict_lower_upper_bounds.find("k8e") != dict_lower_upper_bounds.end()) {
    //        dict_lower_upper_bounds["k8e"] = { 0.000049, 0.00005 };
    //    }
    //    if (dict_lower_upper_bounds.find("lam") != dict_lower_upper_bounds.end()) {
    //        dict_lower_upper_bounds["lam"] = { 1/120., 1/119.9 };
    //    }
    
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
    min_func.f = gsl_EstimateFittedParameters_EDES_Ex_Med;
    // prepare the params to be passed to 'gsl_min_fitted_params_SSR_func'
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = std::make_tuple(this, &dict_lower_upper_bounds);
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
    
    // check iteration: terminate if the results are converging
    const int num_check_iter = 10; // num of iteration results to be averaged over
    double pcnt_curr = 1. / num_check_iter;
    double pcnt_prev = 1- pcnt_curr;
    double avg_check_iter_prev = 100000.;
    double avg_check_iter_curr = 0.;
    double precision_check_iter = 0.5;
    
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
        
        // check iteration: terminate if the results are converging
        avg_check_iter_curr = pcnt_prev * avg_check_iter_prev + pcnt_curr * SSR_curr;
        if (fabs(avg_check_iter_curr - avg_check_iter_prev) < precision_check_iter) break;
        avg_check_iter_prev = avg_check_iter_curr;
        
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


double SubjectGlucose::gsl_EstimateFittedParameters_EDES_Ex_Med (const gsl_vector *v, void *paramsP){
    
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
    sGlucose_tmp.dailyLAInsulinEvents = this2->dailyLAInsulinEvents;
    sGlucose_tmp.dailySAInsulinEvents = this2->dailySAInsulinEvents;
    
    auto num_data_sets = sGlucose_tmp.dailyDataSets.size();
    for (int i = 0; i < num_data_sets; ++i) {
        sGlucose_tmp.eDES_ex_med.ClearPreRunOutput();
        sGlucose_tmp.SetSubjectSpecifics(sGlucose_tmp.dailySubjectSpecifics[i]); // still default parameters according to 'type'
        sGlucose_tmp.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
        sGlucose_tmp.SetFoodIntakeEvents(sGlucose_tmp.dailyFoodIntakeEvents[i]);
        sGlucose_tmp.SetExerciseEvents(sGlucose_tmp.dailyExerciseEvents[i]);
        sGlucose_tmp.SetLaIEvents(sGlucose_tmp.dailyLAInsulinEvents[i]);
        sGlucose_tmp.SetSaIEvents(sGlucose_tmp.dailySAInsulinEvents[i]);
        sGlucose_tmp.EvolutionUnderFoodExInsulinEvents();
        auto check_pts = sGlucose_tmp.ObtainCheckPtsFromDailyDataSets(sGlucose_tmp.dailyDataSets[i]);
        auto glucoses_tmp = sGlucose_tmp.ObtainGlucose(check_pts);
        auto insulins_tmp = sGlucose_tmp.ObtainInsulin(check_pts);
        SSR += sGlucose_tmp.ComputeSSR(sGlucose_tmp.dailyDataSets[i], glucoses_tmp, insulins_tmp);
//        // test
//        for (auto glucose : glucoses_tmp) {
//            std::cout << glucose.first << " " << glucose.second << std::endl;
//        }
    }
    
    return SSR;
}

std::tuple<double, std::map<std::string, double>, std::vector<std::pair<double, double>>> SubjectGlucose::ReconGlucoseCurve(int type, double bodyMass, double glu_pre, double glu_post, double glu_post_time, double food, double timeInterval){
    
    ClearDatasetsForParameterEstimation();
    
    dailyDataSets = {{ {0., glu_pre, 0.1, 0., 0.}, {glu_post_time, glu_post, 0.1, 0., 0.} }};
    dailySubjectSpecifics = {{type, bodyMass, glu_pre, 8., glu_pre, 8.}};
    dailyFoodIntakeEvents = {{ {0., food}, {240., 0.} }};
    dailyExerciseEvents = { {} };
    dailySAInsulinEvents = { {} };
    dailyLAInsulinEvents = { {} } ;
    
    SetFittedParams_EDES_Ex_Med(type);
    
//    std::vector<std::string> chosenFittedParams_str = {"k1", "k2", "k4", "k5", "k6", "k7", "k8", "k9", "KM"};
    
    std::vector<std::string> chosenFittedParams_str = {"k6", "k7", "KM"};
    
    // obtain the initial fitted-params before fitting
    std::map<std::string, double> fittedParamsInit = eDES_ex_med.GetFittedParamsDict();
    
    // set up the dictionary for the fitted params
    std::map<std::string, std::pair<double, double>> dict_lower_upper_bounds;
    // default setting for the allowed ranges of 'fp', [0.01 * currVal, 100 * currVal]
    double lower_tmp, upper_tmp;
    for (int i = 0; i < chosenFittedParams_str.size(); ++i) {
        lower_tmp = 0.01 * fittedParamsInit[chosenFittedParams_str[i]];
        upper_tmp = 100 * fittedParamsInit[chosenFittedParams_str[i]];
        dict_lower_upper_bounds.insert({ chosenFittedParams_str[i], {lower_tmp, upper_tmp} });
    }
    // "k1", "k2", "k4", "k5", "k6", "k7", "k8", "k9"
//    if (dict_lower_upper_bounds.find("k1") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k1"] = { 0.01, 0.1 };
//    }
//    if (dict_lower_upper_bounds.find("k2") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k2"] = { 0.01, 1. };
//    }
//    if (dict_lower_upper_bounds.find("k4") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k4"] = { 0.001, 0.004 };
//    }
//    if (dict_lower_upper_bounds.find("k5") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k5"] = { 0.001, 0.02 };
//    }
    if (dict_lower_upper_bounds.find("k6") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k6"] = { 0.005, 1.5 };
    }
    if (dict_lower_upper_bounds.find("k7") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["k7"] = { 5E-6, 1.0 };
    }
//    if (dict_lower_upper_bounds.find("k8") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k8"] = { 0.1, 20 };
//    }
//    if (dict_lower_upper_bounds.find("k9") != dict_lower_upper_bounds.end()) {
//        dict_lower_upper_bounds["k9"] = { 0.005, 0.1 };
//    }
    if (dict_lower_upper_bounds.find("KM") != dict_lower_upper_bounds.end()) {
        dict_lower_upper_bounds["KM"] = { 5., 100. };
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
    min_func.f = gsl_EstimateFittedParameters_EDES_Ex_Med;
    // prepare the params to be passed to 'gsl_min_fitted_params_SSR_func'
    std::tuple<SubjectGlucose *, std::map<std::string, std::pair<double, double>> *> params_tmp = std::make_tuple(this, &dict_lower_upper_bounds);
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
    
    // check iteration: terminate if the results are converging
    const int num_check_iter = 5; // num of iteration results to be averaged over
    double pcnt_curr = 1. / num_check_iter;
    double pcnt_prev = 1- pcnt_curr;
    double avg_check_iter_prev = 100000.;
    double avg_check_iter_curr = 0.;
    double precision_check_iter = 0.1;
    
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
        
        // check iteration: terminate if the results are converging
        avg_check_iter_curr = pcnt_prev * avg_check_iter_prev + pcnt_curr * SSR_curr;
        if (fabs(avg_check_iter_curr - avg_check_iter_prev) < precision_check_iter) break;
        avg_check_iter_prev = avg_check_iter_curr;
        
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
    
    
    SetSubjectSpecifics(dailySubjectSpecifics[0]);
    SetFittedParams_EDES_Ex_Med(fittedParamsInit);
    SetFoodIntakeEvents(dailyFoodIntakeEvents[0]);
    SetExerciseEvents(dailyExerciseEvents[0]);
    SetSaIEvents(dailySAInsulinEvents[0]);
    SetLaIEvents(dailyLAInsulinEvents[0]);
    
    EvolutionUnderFoodExInsulinEvents();
    
    std::vector<std::pair<double, double>> glucoses = ObtainGlucose(timeInterval);
    
    return std::make_tuple(SSR_curr, fittedParamsInit, glucoses);

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

// ++++++++++++++++++++++++++++++   Methods for analyzing the glucose curve  ++++++++++++++++++++++++++++++

std::tuple<double, double, double, std::vector<double>> SubjectGlucose::AnalyzeGlucoseCurve(const std::vector<std::pair<double, double>> &glucoses){
    
    if (glucoses.size() < 2) {
        std::cout << "ERROR(SubjectGlucose::AnalyzeFourHourGlucose)" << std::endl;
        throw std::runtime_error("ERROR(SubjectGlucose::AnalyzeFourHourGlucose)");
    }
    
    // output: mean, amplitude, glu_range_time_pcnt(<3.9, 3.9-7.8, 7.8-10, >10.)
    double AUC = 0.;
    double min = 100., max = 0.;
    double less_3P9 = 0., btn_3P9_7P8 = 0., btn_7p8_10 = 0., gtr_10 = 0.;
    double timeInterval = glucoses[1].first - glucoses[0].first;
    for (auto iter = glucoses.begin(); iter < glucoses.end()-1; ++iter) {
        double glu_tmp = (iter->second + (iter+1)->second) / 2.;
        AUC += glu_tmp * timeInterval;
        if (glu_tmp < min) min = glu_tmp;
        if (glu_tmp > max) max = glu_tmp;
        if (glu_tmp <= 3.9) less_3P9 += timeInterval;
        if (glu_tmp > 3.9 && glu_tmp <= 7.8) btn_3P9_7P8 += timeInterval;
        if (glu_tmp > 7.8 && glu_tmp <= 10.) btn_7p8_10 += timeInterval;
        if (glu_tmp > 10.) gtr_10 += timeInterval;
    }
    
    auto total_time = (glucoses.size() - 1) * timeInterval;
    double mean = AUC / total_time;
    double amplitude = max - min;
//    std::map<std::string, double> glu_range_time_pcnt = {
//        {"<= 3.9 mmol/L", less_3P9 / total_time},
//        {"3.9 - 7.8 mmol/L", btn_3P9_7P8 / total_time},
//        {"7.8 - 10 mmol/L", btn_7p8_10 / total_time},
//        {"> 10 mmol/L", gtr_10 / total_time}
//    };
    std::vector<double> glu_range_time_pcnt = {less_3P9 / total_time, btn_3P9_7P8 / total_time, btn_7p8_10 / total_time, gtr_10 / total_time};
    
    // HbA1c
    double FBG = glucoses[0].second;
    double mean_24h = (FBG * 12. +  mean * 12.) / 24.;
    double HbA1c = (mean_24h + 0.582) / 1.198; // ref. 糖尿病防治指南
    
//    // GL
//    double GL = (mean - FBG);
    
    return std::make_tuple(mean, HbA1c, amplitude, glu_range_time_pcnt);
    
}

std::tuple<int, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>> SubjectGlucose::EvaluateSeverityFittedGlucoses(double bodyMass, double food, double glu_pre, const std::vector<std::pair<double, double>> &fitted_glucoses){
    
    std::map<std::string, double> modelParams_H = {
        {"k1", 0.144307930814521826E-01},
        {"k2", 0.379549420156179096E+01},
        {"k3", 0.692175703520740747E-03},
        {"k4", 0.202396936120216369E-02},
        {"k5", 0.706987528575793098E-02},
        {"k6", 0.113520049403267409E+01},
        {"k7", 0.439296116630208511E+00},
        {"k8", 0.341210766161778345E+01},
        {"k9", 0.716506564790328065E-01},
        {"k10", 0.05},
        {"k11", 0},
        {"sigma", 0.135575143777190310E+01},
        {"KM", 0.131571530479106595E+02},
        {"k4e", 0},
        {"k5e", 10},
        {"k6e", 0.05},
        {"k8e", 0.5},
        {"lam", 0.00833333}
    };
    
    std::map<std::string, double> modelParams_GFT = {
        {"k1", 0.144307930814521826E-01},
        {"k2", 0.379549420156179096E+01},
        {"k3", 0.692175703520740747E-03},
        {"k4", 0.202396936120216369E-02},
        {"k5", 0.624663357697790612E-02},
        {"k6", 0.586344805285301973E+00},
        {"k7", 0.176202986299337089E-01},
        {"k8", 0.340429479929708823E+01},
        {"k9", 0.716506564790328065E-01},
        {"k10", 0.05},
        {"k11", 0},
        {"sigma", 0.135575143777190310E+01},
        {"KM", 0.325326568142694939E+02},
        {"k4e", 0},
        {"k5e", 10},
        {"k6e", 0.05},
        {"k8e", 0.5},
        {"lam", 0.00833333}
    };
    
    std::map<std::string, double> modelParams_T2Light = {
        {"k1", 0.144307930814521826E-01},
        {"k2", 0.379549420156179096E+01},
        {"k3", 0.692175703520740747E-03},
        {"k4", 0.202396936120216369E-02},
        {"k5", 0.616149969195030692E-02},
        {"k6", 0.586133719919008533E+00},
        {"k7", 0.131910649054100878E-01},
        {"k8", 0.335281340038196740E+01},
        {"k9", 0.716506564790328065E-01},
        {"k10", 0.05},
        {"k11", 0},
        {"sigma", 0.135575143777190310E+01},
        {"KM", 0.348813785020632423E+02},
        {"k4e", 0},
        {"k5e", 10},
        {"k6e", 0.05},
        {"k8e", 0.5},
        {"lam", 0.00833333}
    };
    
    std::map<std::string, double> modelParams_T2Mid = {
        {"k1", 0.144307930814521826E-01},
        {"k2", 0.379549420156179096E+01},
        {"k3", 0.692175703520740747E-03},
        {"k4", 0.202396936120216369E-02},
        {"k5", 0.558564373452744272E-02},
        {"k6", 0.355383834174836988E+00},
        {"k7", 0.143221305922955957E-03},
        {"k8", 0.335281340038196740E+01},
        {"k9", 0.716506564790328065E-01},
        {"k10", 0.05},
        {"k11", 0},
        {"sigma", 0.135575143777190310E+01},
        {"KM", 0.545867398329898208E+02},
        {"k4e", 0},
        {"k5e", 10},
        {"k6e", 0.05},
        {"k8e", 0.5},
        {"lam", 0.00833333}
    };
    
    std::map<std::string, double> modelParams_T2Heavy = {
        {"k1", 0.144307930814521826E-01},
        {"k2", 0.379549420156179096E+01},
        {"k3", 0.692175703520740747E-03},
        {"k4", 0.202396936120216369E-02},
        {"k5", 0.558440318801594650E-02},
        {"k6", 0.120312241302089853E+00},
        {"k7", 0.313625573849650083E-04},
        {"k8", 0.335281340038196740E+01},
        {"k9", 0.716506564790328065E-01},
        {"k10", 0.05},
        {"k11", 0},
        {"sigma", 0.135575143777190310E+01},
        {"KM", 0.546147125708596803E+02},
        {"k4e", 0},
        {"k5e", 10},
        {"k6e", 0.05},
        {"k8e", 0.5},
        {"lam", 0.00833333}
    };
    
    auto f_sz = fitted_glucoses.size();
    double start_time = fitted_glucoses[0].first;
    double final_time = fitted_glucoses[f_sz - 1].first;
    double timeInterval = fitted_glucoses[1].first - fitted_glucoses[0].first;
    
    std::tuple<int, double, double, double, double, double> subjectSpecifics = std::make_tuple(2, bodyMass, glu_pre, 8., glu_pre, 8.);
    
    std::vector<std::pair<double, double>> foodIntakeEvents = {
        {start_time, food},
        {final_time, 0}
    };
    
    std::vector<std::tuple<double, double, double, int>> exerciseEvents = { };
    std::vector<std::tuple<double, int, double>> SAInsulinEvents = {};
    std::vector<std::tuple<double, int, double>> LAInsulinEvents = {};
    
    SetSubjectSpecifics(subjectSpecifics);
    SetFoodIntakeEvents(foodIntakeEvents);
    SetExerciseEvents(exerciseEvents);
    SetSaIEvents(SAInsulinEvents);
    SetLaIEvents(LAInsulinEvents);
    
    // healthy
    SetFittedParams_EDES_Ex_Med(modelParams_H);
    EvolutionUnderFoodExInsulinEvents();
    auto h_glucoses = ObtainGlucose(timeInterval);
    
    // GFT
    SetFittedParams_EDES_Ex_Med(modelParams_GFT);
    EvolutionUnderFoodExInsulinEvents();
    auto GFT_glucoses = ObtainGlucose(timeInterval);
    
    // T2Light
    SetFittedParams_EDES_Ex_Med(modelParams_T2Light);
    EvolutionUnderFoodExInsulinEvents();
    auto T2Light_glucoses = ObtainGlucose(timeInterval);
    
    // T2Mid
    SetFittedParams_EDES_Ex_Med(modelParams_T2Mid);
    EvolutionUnderFoodExInsulinEvents();
    auto T2Mid_glucoses = ObtainGlucose(timeInterval);
    
    // T2Heavy
    SetFittedParams_EDES_Ex_Med(modelParams_T2Heavy);
    EvolutionUnderFoodExInsulinEvents();
    auto T2Heavy_glucoses = ObtainGlucose(timeInterval);
    
    std::vector<double> closeness_glucoses;
    closeness_glucoses.push_back(GlucosesCloseness(fitted_glucoses, h_glucoses));
    closeness_glucoses.push_back(GlucosesCloseness(fitted_glucoses, GFT_glucoses));
    closeness_glucoses.push_back(GlucosesCloseness(fitted_glucoses, T2Light_glucoses));
    closeness_glucoses.push_back(GlucosesCloseness(fitted_glucoses, T2Mid_glucoses));
    closeness_glucoses.push_back(GlucosesCloseness(fitted_glucoses, T2Heavy_glucoses));
    
    int min_type;
    double min = 1E99;
    for (int i = 0; i < closeness_glucoses.size(); ++i) {
        if (min > closeness_glucoses[i]) {
            min = closeness_glucoses[i];
            min_type = i;
        }
    }
    
    auto ret = std::make_tuple(min_type, h_glucoses, GFT_glucoses, T2Light_glucoses, T2Mid_glucoses, T2Heavy_glucoses);
    
    return ret;
}

double SubjectGlucose::GlucosesCloseness(const std::vector<std::pair<double, double>> &fitted_glucoses, const std::vector<std::pair<double, double>> &standard_glucoses){
    
    double ret = 0.;
    
    // compute the SSR normailzed by the number of data points
    auto sz = fitted_glucoses.size();
    for (int i = 0; i < sz; ++i) {
        ret += pow((fitted_glucoses[i].second - standard_glucoses[i].second), 2);
    }
    return ret / sz;
}

// ++++++++++++++++++++++++++++++   Methods for single meal evaluation  ++++++++++++++++++++++++++++++

std::tuple<double, double, double, std::vector<std::pair<std::string, double>>, std::vector<std::pair<std::string, double>>, std::vector<SubjectGlucose::TypeGlucosesOutput>> SubjectGlucose::EvaluateSingleMeal_Glucose(const std::vector<std::pair<std::string, TypeFoodAttrs>> &foodInputLists, double bodyMass, double FBG, const std::map<std::string, double> &modelParams){
    
    TypeFoodAttrs foodTot = 0.;
    for (auto iter = foodInputLists.begin(); iter != foodInputLists.end(); ++iter) foodTot += iter->second;
    
    std::tuple<int, double, double, double, double, double> subjectSpecifics = std::make_tuple(2, bodyMass, FBG, 8., FBG, 8.);
    std::vector<std::tuple<double, double, double, int>> exerciseEvents = { };
    std::vector<std::tuple<double, int, double>> SAInsulinEvents = {};
    std::vector<std::tuple<double, int, double>> LAInsulinEvents = {};
    
    SetSubjectSpecifics(subjectSpecifics);
    SetExerciseEvents(exerciseEvents);
    SetSaIEvents(SAInsulinEvents);
    SetLaIEvents(LAInsulinEvents);
    SetFittedParams_EDES_Ex_Med(modelParams);
    
    double start_time = 0., final_time = 4*60.;
    double timeInterval = 5.;
    
    // all food
    std::vector<std::pair<double, double>> foodIntakeEvents_tmp = { {start_time, foodTot}, {final_time, 0} };
    SetFoodIntakeEvents(foodIntakeEvents_tmp);
    EvolutionUnderFoodExInsulinEvents();
    auto glucoses_allFood = ObtainGlucose(timeInterval);
    auto analysisRes_glu_allFood = AnalyzeGlucoseCurve(glucoses_allFood);
    auto mean = std::get<0>(analysisRes_glu_allFood);
    auto amplitude = std::get<1>(analysisRes_glu_allFood);
    auto PBG2h = glucoses_allFood[glucoses_allFood.size()/2].second;
    
    // loop over each food
    std::vector<std::pair<std::string, double>> mean_foods, amplitude_foods;
    std::vector<TypeGlucosesOutput> glucoses_food_collection;
    for (int i = 0; i < foodInputLists.size(); ++i) {
        foodIntakeEvents_tmp = { {start_time, foodInputLists[i].second}, {final_time, 0} };
        SetFoodIntakeEvents(foodIntakeEvents_tmp);
        EvolutionUnderFoodExInsulinEvents();
        auto glucoses_food = ObtainGlucose(timeInterval);
        auto analysisRes_glu_food = AnalyzeGlucoseCurve(glucoses_food);
        mean_foods.push_back(std::make_pair(foodInputLists[i].first, std::get<0>(analysisRes_glu_food)));
        amplitude_foods.push_back(std::make_pair(foodInputLists[i].first, std::get<1>(analysisRes_glu_food)));
        glucoses_food_collection.push_back(glucoses_food);
    }
    
    // sort mean_foods and amplitude_foods
    std::sort(mean_foods.begin(), mean_foods.end(), [](const std::pair<std::string, double> &a,
                                                       const std::pair<std::string, double> &b) { return a.second > b.second;});
    std::sort(amplitude_foods.begin(), amplitude_foods.end(), [](const std::pair<std::string, double> &a,
                                                       const std::pair<std::string, double> &b) { return a.second > b.second;});
    
    auto ret = std::make_tuple(mean, amplitude, PBG2h, mean_foods, amplitude_foods, glucoses_food_collection);
    
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


