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

#include "E_DES.hpp"

using namespace std;


int main(int argc, const char * argv[]) {
    
    
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
    
    
    // +++++++++++++++++   Using the data sets in Figure 14 of Dr. Mass's thesis
    
    // set up the file path
    string root_dir = "/Users/Jue/Desktop/precision_health/";
    string dpe_file_dir = root_dir + "/data/E_DES/data_parameter_estimation/raw_data/";
    vector<string> dpe_file_names_h = {"model_h_glucose_insulin_with_error"};
    vector<string> dpe_file_names_d2 = {"model_d2_glucose_insulin_with_error"};
    string dpe_glucose_insulin_file_suffix = ".dat";
    vector<string> dpe_glucose_insulin_files_h, dpe_glucose_insulin_files_d2;
    for (auto dpe_file_name: dpe_file_names_h){
        dpe_glucose_insulin_files_h.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_insulin_file_suffix);
    }
    for (auto dpe_file_name: dpe_file_names_d2){
        dpe_glucose_insulin_files_d2.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_insulin_file_suffix);
    }
    
    E_DES eDES;
    
    // fit the data set for the healthy person
    vector<vector<double>> input_parameter_sets = {
        {75E3, 75.} // <foodIntake, bodyMass>
    };
    eDES.SetSubjectTypeFittedParams(0);
    
    
    vector<double> init_fitted_params = eDES.GetFittedParams();
    
    eDES.SetDataForParameterEstimation(dpe_glucose_insulin_files_h, input_parameter_sets);
    
    eDES.EstimateFittedParameters_gsl();
    
    vector<double> optimized_fitted_params = eDES.GetMinFittedParamsSSR();
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H_test.dat", ofstream::out);
    
    cout << "init_fitted_params\t\t\toptimized_params: " << endl;
    for (int i = 0; i < optimized_fitted_params.size() - 1; ++i) {
        cout << init_fitted_params[i] << "\t\t\t" << optimized_fitted_params[i] << endl;
        output << optimized_fitted_params[i] << endl;
    }
    
    output.close();

//    // type-II
//    vector<vector<double>> input_parameter_sets = {
//        {75E3, 83.3} // <foodIntake, bodyMass>
//    };
//    eDES.SetSubjectTypeFittedParams(2);
//    
//    ifstream fitted_params_f;
//    fitted_params_f.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H.dat", ofstream::in);
//    eDES.LoadFittedParams(fitted_params_f);
//    
//    vector<double> init_fitted_params = eDES.GetFittedParams();
//    
//    eDES.SetDataForParameterEstimation(dpe_glucose_insulin_files, input_parameter_sets);
//    
//    eDES.EstimateFittedParameters_gsl();
//    
//    vector<double> optimized_fitted_params = eDES.GetMinFittedParamsSSR();
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/optimized_fitted_params_D2.dat", ofstream::out);
//    
//    cout << "init_fitted_params\t\t\toptimized_params: " << endl;
//    for (int i = 0; i < optimized_fitted_params.size() - 1; ++i) {
//        cout << init_fitted_params[i] << "\t\t\t" << optimized_fitted_params[i] << endl;
//        output << optimized_fitted_params[i] << endl;
//    }
//    
//    output.close();

    
    
    return 0;
}







