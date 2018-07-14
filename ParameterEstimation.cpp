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
    
    // +++++++++++++++++   Using the data sets in Figure 14 of Dr. Mass's thesis
    
    // ------- set up the file path
    string root_dir = "/Users/Jue/Desktop/precision_health/";
    string dataset_dir = root_dir + "/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/";
    vector<string> dataset_names_h = {"fig14_h_glucose_insulin_with_error.dat"};
    vector<string> dataset_names_d2 = {"fig14_d2_glucose_insulin_with_error.dat"};
    vector<string> dataset_files_h, dataset_files_d2;
    for (auto dataset_name_h: dataset_names_h){
        dataset_files_h.push_back(dataset_dir + dataset_name_h );
    }
    for (auto dataset_name_d2: dataset_names_d2){
        dataset_files_d2.push_back(dataset_dir + dataset_name_d2 );
    }
    
    SubjectGlucose sGlucose;
    vector<string> dataset_files;
    
    // ------- healthy person: obtain fitted-params-H from manually tweaked set of fitted-params
    
    dataset_files = dataset_files_h;
    int type = 0;
    double bodyMass = 75.;
    
    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
        {type, bodyMass, 4.55, 16.}
    };
    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
        { {0., 75E3}, {360., 0.} }
    };
    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = {
        { {60., 90., 1., 0} }
    };
    
    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents);
    
    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "k11", "k12", "sigma", "KM"};
    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k11", "k12", "sigma", "KM"};
    auto fitted_params = sGlucose.EstimateFittedParameters_EDES(chosenFittedParams_str);
    
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H_test.dat", ofstream::out);
    
    cout << "optimized_params: " << endl;
    for (int i = 0; i < fittedParamsAll_str.size(); ++i) {
        cout << fittedParamsAll_str[i] << " " << fitted_params[fittedParamsAll_str[i]] << endl;
        output << fitted_params[fittedParamsAll_str[i]] << endl;
    }
    
    output.close();

    
    
    // ------- healthy person: obtain fitted-params-H from manually tweaked set of fitted-params
//
//    vector<vector<double>> input_parameter_sets = {
//        {75E3, 75.} // <foodIntake, bodyMass>
//    };
//    // starting from manually tweaked fitted-params
//    eDES.SetSubjectTypeFittedParams(0);
//    
//    vector<double> init_fitted_params = eDES.GetFittedParams();
//    
//    eDES.SetDataForParameterEstimation(dpe_glucose_insulin_files_h, input_parameter_sets);
//    
//    eDES.EstimateFittedParameters_gsl();
//    
//    vector<double> optimized_fitted_params = eDES.GetMinFittedParamsSSR();
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H.dat", ofstream::out);
//
//    cout << "init_fitted_params\t\t\toptimized_params: " << endl;
//    for (int i = 0; i < optimized_fitted_params.size() - 1; ++i) {
//        cout << init_fitted_params[i] << "\t\t\t" << optimized_fitted_params[i] << endl;
//        output << optimized_fitted_params[i] << endl;
//    }
//    
//    output.close();
    
//    // ------- D2: obtain fitted-params-D2 from optimized fitted-params-H
//
//    vector<vector<double>> input_parameter_sets = {
//        {75E3, 83.3} // <foodIntake, bodyMass>
//    };
//    // starting from fitted-params-H
//    ifstream fitted_params_H;
//    fitted_params_H.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H.dat", ofstream::in);
//    eDES.LoadFittedParams(fitted_params_H);
//    
//    vector<double> init_fitted_params = eDES.GetFittedParams();
//    
//    eDES.SetDataForParameterEstimation(dpe_glucose_insulin_files_d2, input_parameter_sets);
//    
//    eDES.EstimateFittedParameters_gsl();
//    
//    vector<double> optimized_fitted_params = eDES.GetMinFittedParamsSSR();
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2_test.dat",
//                ofstream::out);
//    
//    cout << "init_fitted_params\t\t\toptimized_params: " << endl;
//    for (int i = 0; i < optimized_fitted_params.size() - 1; ++i) {
//        cout << init_fitted_params[i] << "\t\t\t" << optimized_fitted_params[i] << endl;
//        output << optimized_fitted_params[i] << endl;
//    }
//    
//    output.close();
    
    // +++++++++++++++++   Using multiple data sets in Table 2 of the EDES paper
    
    //    string root_dir = "/Users/Jue/Desktop/precision_health/";
    
    //    ifstream data_para_est_file;
    //    string dpe_file_dir = root_dir + "/data/E_DES/data_parameter_estimation/clean_data/";
    //    vector<string> dpe_file_names = {"34_female", "34_male", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44"};
    //    string dpe_glucose_file_suffix = "_plasma_glucose_cleaned.dat";
    //    string dpe_insulin_file_suffix = "_plasma_insulin_cleaned.dat";
    //    vector<string> dpe_glucose_files, dpe_insulin_files;
    //    for (auto dpe_file_name: dpe_file_names){
    //        dpe_glucose_files.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_file_suffix);
    //        dpe_insulin_files.push_back(dpe_file_dir + dpe_file_name + dpe_insulin_file_suffix);
    //    }
    

    
    
    return 0;
}







