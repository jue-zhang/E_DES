//
//  main.cpp
//  E_DES
//
//  Created by Jue on 6/24/18.
//  Copyright © 2018 Jue. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "SubjectGlucose.hpp"

using namespace std;


int main(int argc, const char * argv[]) {
    
    // +++++++++++++++++   Using the fitted curve in Figure 14 of Dr. Mass's thesis
    
    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/fig14_d2_glucose_insulin_with_error.dat";
    
    SubjectGlucose sGlucose;
    vector<string> dataset_files;
    
    dataset_files = {dataset_file_path};
    
    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
        {2 , 83.3, 5.46, 8.63, 4.55, 16.}
    };
    
    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
        { {0., 75E3}, {360., 0.} }
    };
    
    
    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
    
    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
    
    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;

    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);
    
//    sGlucose.SetFittedParams_EDES_Ex_Med(0);
    
    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig17_h.dat";
    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
    
    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
    vector<string> chosenFittedParams_str = {"k5", "k6", "k7", "k8", "KM"};
    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
    
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig14_d2.dat", ofstream::out);
    
    // +++++++++++++++++   Using the fitted curve in Figure 14 of Dr. Mass's thesis
    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/fig14_h_glucose_insulin_with_error.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {0 , 75., 4.55, 16., 4.55, 16.}
//    };
//    
//    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
//        { {0., 75E3}, {360., 0.} }
//    };
//    
//    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;
//
//    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);
//
//    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM"};
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "sigma", "KM"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
    
//    // +++++++++++++++++   Using the fitted curve in Figure 17 of Dr. Mass's thesis
//    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/fig17_fitted_glucose_insulin_with_error.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {0 , 75., 5., 7.4, 5., 7.4}
//    };
//    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
//        { {0., 75E3}, {360., 0.} }
//    };
//    
//    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;
//    
//    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);
//    
//    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "sigma", "KM"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig17_h.dat", ofstream::out);
    
    cout << "optimized_params: " << endl;
    for (int i = 0; i < fittedParamsAll_str.size(); ++i) {
        cout << fittedParamsAll_str[i] << " " << fitted_params[fittedParamsAll_str[i]] << endl;
        output << fitted_params[fittedParamsAll_str[i]] << endl;
    }
    
    output.close();

    
    
    return 0;
}







