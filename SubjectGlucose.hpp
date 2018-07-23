//
//  SubjectGlucose.hpp
//  E_DES
//
//  Created by Jue on 7/11/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#ifndef SubjectGlucose_hpp
#define SubjectGlucose_hpp

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <map>

#include "E_DES.hpp"
#include "E_DES_Ex.hpp"
#include "E_DES_Ex_Med.hpp"

class SubjectGlucose{
    
    // ++++++++++++++++++++++++++++++   Short-hand type declarations  ++++++++++++++++++++++++++++++
public:
//    // type of subject specifics: <type, bodyMass, fastGlucose, fastInsulin>
//    using TypeSubjectSpecifics = std::tuple<int, double, double, double>;
    
    // type of subject specifics: <type, bodyMass, fastGlucose, fastInsulin, Gpl_basal, Ipl_basal>
    using TypeSubjectSpecifics = std::tuple<int, double, double, double, double, double>;
    
    // type of a food intake event: <t_intake, foodIntake>
    using TypeFoodIntakeEvent = std::pair<double, double>;
    // type of an exercise event: <t_start, t_finish, intensity, type>, where 't_start' ('t_finish') is the time
    //      instant when the exercise starts (finishes), 'intensity' is the level of the intensity of the exercise,
    //      and 'type' is the type of exercise.
    using TypeExerciseEvent = std::tuple<double, double, double, int>;
    // type of a combined food intake and exercise event: <t, this_meal, latest_meal, delta_eat, exIntensity, exType>
    using TypeFoodIntakeExerciseEvent = std::tuple<double, double, double, double, double, int>;
    // type of a short-acting insulin event: <t_intake, sa_insulin_type, sa_insulin>
    using TypeSAInsulinEvent = std::tuple<double, int, double>;
    // type of a long-acting insulin event: <t_intake, la_insulin_type, la_insulin>
    using TypeLAInsulinEvent = std::tuple<double, int, double>;
    // type of a combined food intake and exercise event: <t, this_meal, latest_meal, delta_eat, exIntensity, exType, insulin_sa_type, insulin_sa>
    using TypeFoodExSaIEvent = std::tuple<double, double, double, double, double, int, int, double>;
    // type of a combined food intake and exercise event: <t, this_meal, latest_meal, delta_eat, exIntensity, exType, insulin_sa_type, insulin_sa, insulin_la_type_curr, insulin_la_curr, insulin_la_type_latest, insulin_la_latest, delta_insulin_la>
    using TypeFoodExSaILaIEvent = std::tuple<double, double, double, double, double, int, int, double, int, double, int, double, double>;
    
    // type of a collection of daily food intake events
    using TypeDailyFoodIntakeEvents = std::vector<TypeFoodIntakeEvent>;
    // type of a collection of daily exercise events
    using TypeDailyExerciseEvents = std::vector<TypeExerciseEvent>;
    // type of a collection of short-acting insulin events
    using TypeDailySAInsulinEvents = std::vector<TypeSAInsulinEvent>;
    // type of a collection of long-acting insulin events
    using TypeDailyLAInsulinEvents = std::vector<TypeLAInsulinEvent>;
    // type of a collection of daily food intake and exercise events
    using TypeDailyFoodIntakeExerciseEvents = std::vector<TypeFoodIntakeExerciseEvent>;
    // type of a collection of daily food intake, exercise, saI events
    using TypeDailyFoodIntakeExerciseSaIEvents = std::vector<TypeFoodExSaIEvent>;
    // type of a collection of daily food intake, exercise, med events
    using TypeDailyFoodIntakeExerciseSaILaIEvents = std::vector<TypeFoodExSaILaIEvent>;
    
    // type of each row in personal data set: <t[i], glucose[i], glu_err[i], insulin[i], ins_err[i]>
    using TypeDataSetRow = std::tuple<double, double, double, double, double>;
    // type of a single personal data set
    using TypeDataSet = std::vector<TypeDataSetRow>;
    
    // ++++++++++++++++++++++++++++++   Methods for evolution  ++++++++++++++++++++++++++++++
public:
    
    // Method:
    //      Perform the evolution under the given foodIntakeEvents and exerciseEvents
    // Outcome:
    //      Store the evolution results in 'spline' for later interpolation
    void EvolutionUnderFoodIntakeExerciseEvents();

    // Method:
    //      Obtain the glucose levels in the user-specified intervals under a set of consective food intake and excercise events
    // Input:
    //      timeInterval -- the user specified timeInterval
    // Output:
    //      the glucose levels, in the form of <time, glucose>, in the user-specified intervals,
    //      i.e., {t_init, t_init + timeInterval, t_init + 2 * timeInterval, ..., t_end }.
    std::vector<std::pair<double, double>> GlucoseUnderFoodIntakeExerciseEvents(double timeInterval);
    std::vector<std::pair<double, double>> InsulinUnderFoodIntakeExerciseEvents(double timeInterval);
    
    // Method:
    //      Obtain the glucose levels according to 'check_pts'
    std::vector<std::pair<double, double>> GlucoseUnderFoodIntakeExerciseEvents(const std::vector<double> &check_pts);
    std::vector<std::pair<double, double>> InsulinUnderFoodIntakeExerciseEvents(const std::vector<double> &check_pts);
    
    
//    std::vector<std::pair<double, double>> GetGlucoses(); // return the obtained glucoses at the specied time instants
    
private:
    void ClearPreEvolutionOutput();
    void CombineFoodIntakeExerciseEvents();
    void CombineFoodIntakeExerciseSaIEvents();
    void CombineFoodIntakeExerciseSaILaIEvents();
    std::pair<double, double> FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseEvents &f_ex_ents);
    std::pair<double, double> FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseSaIEvents &f_ex_saI_ents);
    std::pair<double, double> FindTheLatestFoodIntakeEvents(const TypeDailyFoodIntakeExerciseSaILaIEvents &f_ex_saI_laI_ents);
    std::tuple<int, double, double> FindTheLatestLaIEvents(const TypeDailyFoodIntakeExerciseSaILaIEvents &f_ex_saI_laI_ents);
    
    
    // ++++++++++++++++++++++++++++++   Methods for modifying/extracting model parameters  ++++++++++++++++++++++++++++++
    
public:
    // set subject specifics <type, bodyMass, fastGlucose, fastInsulin>, 'foodIntakeEvents', 'exerciseEvents'
    void SetSubjectSpecifics(const TypeSubjectSpecifics &specifics);
    void SetFoodIntakeEvents(const std::vector<TypeFoodIntakeEvent> &foodIntakeEvents_i) { foodIntakeEvents =foodIntakeEvents_i; }
    void SetExerciseEvents (const std::vector<TypeExerciseEvent> &exerciseEvents_i) { exerciseEvents = exerciseEvents_i;}
    void SetSaIEvents(const TypeDailySAInsulinEvents &SAInsulinEvents_i) {SAInsulinEvents = SAInsulinEvents_i;}
    void SetLaIEvents(const TypeDailyLAInsulinEvents &LAInsulinEvents_i) {LAInsulinEvents = LAInsulinEvents_i;}

public:
    // set the fitted params in EDES
    void SetFittedParams_EDES(const int &type) { eDES.SetFittedParams(type); }
    void SetFittedParams_EDES(std::ifstream &param_file) {eDES.SetFittedParams(param_file);}
    void SetFittedParams_EDES(const std::map<std::string, double> &dict_fitted_params) { eDES.SetFittedParams(dict_fitted_params);}
    
    // set the fitted params in EDES_Ex
    void SetFittedParams_EDES_Ex(const int &type) { eDES_ex.SetFittedParams(type); }
    //    void SetFittedParams_EDES_Ex(std::ifstream &param_file) {eDES_ex.SetFittedParams(param_file);}
    void SetFittedParams_EDES_Ex(const std::string &param_file_path) {eDES_ex.SetFittedParams(param_file_path);}
    void SetFittedParams_EDES_Ex(const std::map<std::string, double> &dict_fitted_params) {
        eDES_ex.SetFittedParams(dict_fitted_params);
    }
    
    // set the fitted params in EDES_Ex_Med
    void SetFittedParams_EDES_Ex_Med(const int &type) { eDES_ex_med.SetFittedParams(type); }
    void SetFittedParams_EDES_Ex_Med(const std::string &param_file_path) {eDES_ex_med.SetFittedParams(param_file_path);}
    void SetFittedParams_EDES_Ex_Med(const std::map<std::string, double> &dict_fitted_params) {
        eDES_ex_med.SetFittedParams(dict_fitted_params);
    }
    
    // ++++++++++++++++++++++++++++++   Methods for estimating fitted-params  ++++++++++++++++++++++++++++++
public:
    
    // Method:
    //      Set up the personal data sets used for estimating fitted-params
    // Input:
    //      dataSets -- paths of data files
    //          file format: every row "i", "t[i] glucose[i] glu_err[i] insulin[i] ins_err[i]"
    //      subjectSpecificsCol -- the subject specific params, e.g., bodyMass, for all data sets
    //      foodIntakeEventsCol -- foodIntakeEvents for all data sets
    //      exerciseEventsCol   -- exerciseEvents for all data sets
    // Output:
    //      'minParams_fval_curr' stores all the final estimated fitted-params
    void SetDatasetsForParameterEstimation(const std::vector<std::string> &dailyDataSets_i,
                                           const std::vector<TypeSubjectSpecifics> &dailySubjectSpecifics_i,
                                           const std::vector<TypeDailyFoodIntakeEvents> &dailyFoodIntakeEvents_i,
                                           const std::vector<TypeDailyExerciseEvents> &dailyExerciseEvents_i);
    void SetDatasetsForParameterEstimation(const std::vector<std::string> &dailyDataSets_i,
                                           const std::vector<TypeSubjectSpecifics> &dailySubjectSpecifics_i,
                                           const std::vector<TypeDailyFoodIntakeEvents> &dailyFoodIntakeEvents_i,
                                           const std::vector<TypeDailyExerciseEvents> &dailyExerciseEvents_i,
                                           const std::vector<TypeDailySAInsulinEvents> &dailySAInsulinEvents_i,
                                           const std::vector<TypeDailyLAInsulinEvents> &dailyLAInsulinEvents_i);
    
    void ClearDatasetsForParameterEstimation();
    
    
//    // Method:
//    //      Estimate the chosen set of fitted parameters
//    // Input:
//    //      chosenFittedParams -- the chosen set of fitted params, a vector of boolean values corresponding to each fitted-params
//    //      Order of fitted-params:
//    //      0  1  2  3  4  5  6  7  8  9   10  11  12    13
//    //      |  |  |  |  |  |  |  |  |  |   |   |   |     |
//    //      k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 sigma KM
//    void EstimateFittedParameters_EDES (const std::vector<bool> &chosenFittedParams);
    
    // Directly input the names of the fitted parameters that are to be minimized
    std::map<std::string, double> EstimateFittedParameters_EDES (const std::vector<std::string> &chosenFittedParams_str);
    std::map<std::string, double> EstimateFittedParameters_EDES_Ex (const std::vector<std::string> &chosenFittedParams_str);
    std::map<std::string, double> EstimateFittedParameters_EDES_Ex_Med (const std::vector<std::string> &chosenFittedParams_str);
    
private:
    // internal usage for minimizing the fitted parameters
    static double gsl_EstimateFittedParameters_EDES (const gsl_vector *v, void *paramsP);
    static double gsl_EstimateFittedParameters_EDES_Ex (const gsl_vector *v, void *paramsP);
    static double gsl_EstimateFittedParameters_EDES_Ex_Med (const gsl_vector *v, void *paramsP);
    double ComputeSSR(const TypeDataSet &param_est_data_set,
                      std::vector<std::pair<double, double>> &glucoses,
                      std::vector<std::pair<double, double>> &insulins);
    std::vector<double> ObtainCheckPtsFromDailyDataSets( const TypeDataSet &dataSet);

    // ++++++++++++++++++++++++++++++   Methods for class itself  ++++++++++++++++++++++++++++++
public:
    SubjectGlucose() = default;
    ~SubjectGlucose();
    
    // ++++++++++++++++++++++++++++++   SubjectGlucose Parameters  ++++++++++++++++++++++++++++++
private:
    
    // ---------------- specifics of the subject
    
    int type = 0; // type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    double bodyMass = 75.;
    double fastGlucose = 5.; // the fast glucose level, mmol/L
    double fastInsulin = 8.; // the fast insulin level, U/L
    double Gpl_basal = 5.; // targeted basal glucose level
    double Ipl_basal = 8.; // targeted basal insulin level
    
    // --------------- input info about daily foodIntakeEvents and exerciseEvents
    //
    // foodIntakeEvents -- a list of food intake events in the form of <time, foodIntake>, where "time" is the food intake
    //              time instant, and 'foodIntake' is the amount of intaked food, in unit of 'mg'. The first element of
    //              'foodIntakeEvents' is the initial time, i.e., foodIntakeEvents[0] = {t_init, initial_food_intake},
    //              while the last element is the final time instant. Note: If we want results at 8h after the final food intake,
    //              we should set the last element in 'foodIntakeEvents' as {tF + 8*60, 0.}, where 'tF' is the time instant
    //              of the final NON-ZERO food intake event.
    // exerciseEvents -- a list of exercise events in the form of <t_start, t_finish, intensity, type>,
    //              where 't_start' ('t_finish') is the time instant when the exercise starts (finishes),
    //              'intensity' is the level of the intensity of the exercise, and 'type' is the type of exercise.
    // foodIntakeExerciseEvents -- a list of cominbed foodIntake and ex. events: <t, foodIntake, isExercise, exIntensity, exType>
    TypeDailyFoodIntakeEvents foodIntakeEvents;
    TypeDailyExerciseEvents exerciseEvents;
    TypeDailySAInsulinEvents SAInsulinEvents;
    TypeDailyLAInsulinEvents LAInsulinEvents;
    TypeDailyFoodIntakeExerciseEvents foodIntakeExerciseEvents;
    TypeDailyFoodIntakeExerciseSaIEvents foodIntakeExerciseSaIEvents;
    TypeDailyFoodIntakeExerciseSaILaIEvents foodIntakeExerciseSaILaIEvents;
    
    // --------------- Model used to perform evolution
    E_DES eDES;
    E_DES_Ex eDES_ex;
    E_DES_Ex_Med eDES_ex_med;
    
    //---------------- evolved glucose and insulin levels under the foodIntakeEvents and exerciseEvents
    
    // output the obtained glucoses and insulins at the specified time instants in 'time_instants'
    std::vector<double> time_instants;
    std::vector<double> glucoses;
    std::vector<double> insulins;
    // for interpolation
    double *time_instants_interp = nullptr;
    double *glucoses_interp = nullptr;
    gsl_spline *glucose_spline = nullptr; // store the spline for the later interpolation
    gsl_interp_accel *glucose_acc = nullptr; // for fast interpolation
    double *insulins_interp = nullptr;
    gsl_spline *insulin_spline = nullptr; // store the spline for the later interpolation
    gsl_interp_accel *insulin_acc = nullptr; // for fast interpolation
    
    
    // ---------------- Vars used for estimating fitted-params
    
    // store a collection of personal data sets and the corresponding subject specifics, foodIntakeEvents, exerciseEvents
    std::vector<TypeDataSet> dailyDataSets;
    std::vector<TypeSubjectSpecifics> dailySubjectSpecifics;
    std::vector<TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents;
    std::vector<TypeDailyExerciseEvents> dailyExerciseEvents;
    std::vector<TypeDailySAInsulinEvents> dailySAInsulinEvents;
    std::vector<TypeDailyLAInsulinEvents> dailyLAInsulinEvents;
    
    
    
};


#endif /* SubjectGlucose_hpp */
