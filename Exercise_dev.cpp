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
    
    // +++++++++++++++++   Using the data sets in Fig.2 of the following paper
    // ref. https://academic.oup.com/jcem/article/99/9/3334/2538522
    
    // ------- set up the file path
    string root_dir = "/Users/Jue/Desktop/precision_health/";
    string dpe_file_dir = root_dir + "models/E_DES/E_DES/datasets_for_param_estimation/exercise/";
    vector<string> dpe_file_names_con = {"Karstoft_con_glucose_insulin_mean_error.dat"};
    vector<string> dpe_file_names_IW = {"Karstoft_IW_glucose_insulin_mean_error_M_60.dat"};
    vector<string> dpe_file_names_CW = {"Karstoft_CW_glucose_insulin_mean_error.dat"};
    vector<string> dpe_glucose_insulin_files_con, dpe_glucose_insulin_files_IW, dpe_glucose_insulin_files_CW;
    for (auto dpe_file_name: dpe_file_names_con){
        dpe_glucose_insulin_files_con.push_back(dpe_file_dir + dpe_file_name);
    }
    for (auto dpe_file_name: dpe_file_names_IW){
        dpe_glucose_insulin_files_IW.push_back(dpe_file_dir + dpe_file_name);
    }
    for (auto dpe_file_name: dpe_file_names_CW){
        dpe_glucose_insulin_files_CW.push_back(dpe_file_dir + dpe_file_name);
    }
    
    SubjectGlucose sGlucose;
    vector<string> dataset_files;
    
//    // ------- control
//    
//    // set the vars. specific to the suject for all data sets
//    
//    dataset_files = dpe_glucose_insulin_files_con;
//    int type = 2;
//    double bodyMass = 85.9;
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {type, bodyMass, 7.96, 12.84}
//    };
//    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
//        { {0., 0.}, {105., 75E3},  {345., 0.} }
//    };
//    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = {
//        { {60., 90., 1., 0} }
//    };
//
//    
//    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "k12", "sigma", "KM"};
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k11", "k12", "sigma", "KM"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2_Karstoft_con.dat", ofstream::out);
//    
//    cout << "optimized_params: " << endl;
//    for (int i = 0; i < fittedParamsAll_str.size(); ++i) {
//        cout << fittedParamsAll_str[i] << " " << fitted_params[fittedParamsAll_str[i]] << endl;
//        output << fitted_params[fittedParamsAll_str[i]] << endl;
//    }
//    
//    output.close();
    
    // ------- IW
    
    // set the vars. specific to the suject for all data sets
    
    dataset_files = dpe_glucose_insulin_files_IW;
    int type = 2;
    double bodyMass = 85.9;
    
    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
        {type, bodyMass, 6.7, 7.65}
    };
    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
        { {0., 0.}, {105., 75E3},  {345., 0.} }
    };
    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = {
        { {0., 60., 1., 0} }
    };
    
    
    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents);
    
    // starting from parameters obtained from fitting "con"
    string param_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2_Karstoft_con.dat";
    sGlucose.SetFittedParams_EDES_Ex(param_file_path);
    
    
    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "k11", "k12", "sigma", "KM", "k4e", "k5e", "k8e", "lam"};
    vector<string> chosenFittedParams_str = {"k4e", "k5e", "k8e", "lam"};
    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex(chosenFittedParams_str);
    
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2_Karstoft_IW.dat", ofstream::out);
    
    cout << "optimized_params: " << endl;
    for (int i = 0; i < fittedParamsAll_str.size(); ++i) {
        cout << fittedParamsAll_str[i] << " " << fitted_params[fittedParamsAll_str[i]] << endl;
        output << fitted_params[fittedParamsAll_str[i]] << endl;
    }
    
    output.close();

    
    
    return 0;
}







