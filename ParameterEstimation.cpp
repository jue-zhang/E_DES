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
    
    // +++++++++++++++++   Using the data sets in Figure 14 of Dr. Mass's thesis
    
    // ------- set up the file path
    string root_dir = "/Users/Jue/Desktop/precision_health/";
    string dpe_file_dir = root_dir + "/data/E_DES/data_parameter_estimation/raw_data/";
    vector<string> dpe_file_names_h = {"model_h_glucose_insulin_with_error_M"};
    vector<string> dpe_file_names_d2 = {"model_d2_glucose_insulin_with_error_M"};
    string dpe_glucose_insulin_file_suffix = ".dat";
    vector<string> dpe_glucose_insulin_files_h, dpe_glucose_insulin_files_d2;
    for (auto dpe_file_name: dpe_file_names_h){
        dpe_glucose_insulin_files_h.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_insulin_file_suffix);
    }
    for (auto dpe_file_name: dpe_file_names_d2){
        dpe_glucose_insulin_files_d2.push_back(dpe_file_dir + dpe_file_name + dpe_glucose_insulin_file_suffix);
    }
    
    E_DES eDES;
    
//    // ------- healthy person: obtain fitted-params-H from manually tweaked set of fitted-params
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
    
    // ------- D2: obtain fitted-params-D2 from optimized fitted-params-H

    vector<vector<double>> input_parameter_sets = {
        {75E3, 83.3} // <foodIntake, bodyMass>
    };
    // starting from fitted-params-H
    ifstream fitted_params_H;
    fitted_params_H.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_H.dat", ofstream::in);
    eDES.LoadFittedParams(fitted_params_H);
    
    vector<double> init_fitted_params = eDES.GetFittedParams();
    
    eDES.SetDataForParameterEstimation(dpe_glucose_insulin_files_d2, input_parameter_sets);
    
    eDES.EstimateFittedParameters_gsl();
    
    vector<double> optimized_fitted_params = eDES.GetMinFittedParamsSSR();
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2.dat",
                ofstream::out);
    
    cout << "init_fitted_params\t\t\toptimized_params: " << endl;
    for (int i = 0; i < optimized_fitted_params.size() - 1; ++i) {
        cout << init_fitted_params[i] << "\t\t\t" << optimized_fitted_params[i] << endl;
        output << optimized_fitted_params[i] << endl;
    }
    
    output.close();
    
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







