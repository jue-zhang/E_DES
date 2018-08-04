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
#include <time.h>

#include "SubjectGlucose.hpp"

using namespace std;


int main(int argc, const char * argv[]) {
    
    clock_t start, finish;
    
    start = clock();
    
    cout << "start...." << endl;
    
    // +++++++++++++++++   report_data: T2Mid
    
    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/T2Mid_data.dat";
    
    SubjectGlucose sGlucose;
    vector<string> dataset_files;
    
    dataset_files = {dataset_file_path};
    
    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
        {2 , 83.3, 8., 8., 8., 8.}
    };
    
    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
        { {0., 75E3}, {360., 0.} }
    };
    
    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
    
    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
    
    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;
    
    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);
    
    //    sGlucose.SetFittedParams_EDES_Ex_Med(0);
    
    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_h.dat";
    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
    
    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
    //    vector<string> chosenFittedParams_str = {"k5", "k6", "k7", "k8", "KM"};//k5, k6, k7, k8, KM
    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "sigma"};
    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
    
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_T2Mid.dat", ofstream::out);
    
//    // +++++++++++++++++   report_data: T2Light
//    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/T2Light_data.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 6.84, 8., 6.84, 8.}
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
//    //    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_h.dat";
//    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
//    //    vector<string> chosenFittedParams_str = {"k5", "k6", "k7", "k8", "KM"};//k5, k6, k7, k8, KM
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "sigma"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_T2Light.dat", ofstream::out);
    
//    // +++++++++++++++++   report_data: GFT
//    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/GFT_data.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 6.4, 8., 6.4, 8.}
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
//    //    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_h.dat";
//    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
////    vector<string> chosenFittedParams_str = {"k5", "k6", "k7", "k8", "KM"};//k5, k6, k7, k8, KM
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "sigma"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_GFT.dat", ofstream::out);

    
//    // +++++++++++++++++   report_data: health
//    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/health_data.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 5., 8., 5., 8.}
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
//    //    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig17_h.dat";
//    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "sigma"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_report_h.dat", ofstream::out);

    
//    // +++++++++++++++++   Using the fitted curve in Figure 14 of Dr. Mass's thesis
//    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/fig14_d2_glucose_insulin_with_error.dat";
//    
//    SubjectGlucose sGlucose;
//    vector<string> dataset_files;
//    
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 5.46, 8.63, 4.55, 16.}
//    };
//    
//    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
//        { {0., 75E3}, {360., 0.} }
//    };
//    
//    
//    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
//    
//    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;
//
//    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);
//    
////    sGlucose.SetFittedParams_EDES_Ex_Med(0);
//    
//    string fitted_params_h = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig17_h.dat";
//    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params_h);
//    
//    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
//        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
////    vector<string> chosenFittedParams_str = {"k5", "k6", "k7", "k8", "KM"};
//    vector<string> chosenFittedParams_str = {"k1", "k2", "k6", "k7", "k8", "k10"};
//    auto fitted_params = sGlucose.EstimateFittedParameters_EDES_Ex_Med(chosenFittedParams_str);
//    
//    ofstream output;
//    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_thesis_fig14_d2_test.dat", ofstream::out);
    
//    // +++++++++++++++++   Using the fitted curve in Figure 14 of Dr. Mass's thesis
//    
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

    sGlucose.SetSubjectSpecifics(dailySubjectSpecifics[0]);
    sGlucose.SetFittedParams_EDES_Ex_Med(fitted_params);
    sGlucose.SetFoodIntakeEvents(dailyFoodIntakeEvents[0]);
    sGlucose.SetExerciseEvents(dailyExerciseEvents[0]);
    sGlucose.SetSaIEvents(dailySAInsulinEvents[0]);
    sGlucose.SetLaIEvents(dailyLAInsulinEvents[0]);
    
    sGlucose.EvolutionUnderFoodExInsulinEvents();
    
    // extract the obtained glucose levels at the user specified time intervals
    double timeInterval = 5.;
    std::vector<std::pair<double, double>> glucoses = sGlucose.ObtainGlucose(timeInterval);
    
    // export results:
    string output_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/";
    output.open(output_file_path + "output.dat", ofstream::out);
    for (auto glucose: glucoses) {
        cout << glucose.first / 60.   << " " << glucose.second << endl;
        output << glucose.first  << " " << glucose.second << endl;
    }
    output.close();
    
    finish = clock();
    cout << "Time elapsed (sec): " << (double) (finish-start)/CLOCKS_PER_SEC << endl;
    
    return 0;
}







