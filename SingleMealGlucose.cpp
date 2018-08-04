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
    
    SubjectGlucose sGlucose;
    int type = 2;
    double bodyMass = 83.;
//    double glu_pre = 6.5, glu_post = 9.; // GFT
//    double glu_pre = 7., glu_post = 9.5; // T2Light
//    double glu_pre = 8., glu_post = 12.; // T2Mid
//    double glu_pre = 9., glu_post = 15.; // T2Heavy
    double glu_pre = 8.5, glu_post = 13.; // arbitrary
    double glu_post_time = 120.;
    double food = 75E3;
    double timeInterval = 5.;
    auto res = sGlucose.ReconGlucoseCurve(type, bodyMass, glu_pre, glu_post, glu_post_time, food, timeInterval);
    
    auto SSR = get<0>(res);
    auto fitted_params = get<1>(res);
    auto glucoses = get<2>(res);
    
    cout << "SSR: " << SSR << endl;
    
    
    // evaluating the severity of glucose levels
    auto evaluation = sGlucose.EvaluateSeverityFittedGlucoses(bodyMass, food, glu_pre, glucoses);
    cout << "The paitent is likely to be : " << get<0>(evaluation) << endl;
    
    ofstream std_glu_output;
    // healthy
    auto std_glucoses = get<1>(evaluation);
    std_glu_output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/output_h.dat", ofstream::out);
    for (auto glucose: std_glucoses) {
        std_glu_output << glucose.first  << " " << glucose.second << endl;
    }
    std_glu_output.close();
    // GFT
    std_glucoses = get<2>(evaluation);
    std_glu_output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/output_GFT.dat", ofstream::out);
    for (auto glucose: std_glucoses) {
        std_glu_output << glucose.first  << " " << glucose.second << endl;
    }
    std_glu_output.close();
    // T2Light
    std_glucoses = get<3>(evaluation);
    std_glu_output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/output_T2Light.dat", ofstream::out);
    for (auto glucose: std_glucoses) {
        std_glu_output << glucose.first  << " " << glucose.second << endl;
    }
    std_glu_output.close();
    // T2Mid
    std_glucoses = get<4>(evaluation);
    std_glu_output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/output_T2Mid.dat", ofstream::out);
    for (auto glucose: std_glucoses) {
        std_glu_output << glucose.first  << " " << glucose.second << endl;
    }
    std_glu_output.close();
    // T2Heavy
    std_glucoses = get<5>(evaluation);
    std_glu_output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/output_T2Heavy.dat", ofstream::out);
    for (auto glucose: std_glucoses) {
        std_glu_output << glucose.first  << " " << glucose.second << endl;
    }
    std_glu_output.close();
    
    
    ofstream output;
    output.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/single_meal_test.dat", ofstream::out);
    
    vector<string> fittedParamsAll_str = {"k1", "k2", "k3", "k4", "k5", "k6", "k7",
        "k8", "k9", "k10", "k11", "sigma", "KM", "k4e", "k5e", "k6e", "k8e", "lam"};
    
    cout << "optimized_params: " << endl;
    for (int i = 0; i < fittedParamsAll_str.size(); ++i) {
        cout << fittedParamsAll_str[i] << " " << fitted_params[fittedParamsAll_str[i]] << endl;
        output << fitted_params[fittedParamsAll_str[i]] << endl;
    }
    output.close();
    
    string output_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/";
    output.open(output_file_path + "output.dat", ofstream::out);
    for (auto glucose: glucoses) {
        cout << glucose.first / 60.   << " " << glucose.second << endl;
        output << glucose.first  << " " << glucose.second << endl;
    }
    output.close();
    
    
    cout << "single meal evaluation...." << endl;
    std::vector<std::pair<std::string, double>> foodInputLists = {
        {"foodA", 65E3},
        {"foodB", 10E3}
    };
    auto singleMealEvaluation = sGlucose.EvaluateSingleMeal_Glucose(foodInputLists, bodyMass, glu_pre, fitted_params);
    cout << "mean = " << get<0>(singleMealEvaluation) << endl;
    cout << "amplitude = " << get<1>(singleMealEvaluation) << endl;
    cout << "PBG2h = " << get<2>(singleMealEvaluation) << endl;
    cout << "food ordered by mean: " << endl;
    auto orderedFoodList_mean = get<3>(singleMealEvaluation);
    for (auto food_: orderedFoodList_mean) {
        cout << food_.first << " " << food_.second << endl;
    }
    cout << "food ordered by amp: " << endl;
    auto orderedFoodList_amp = get<4>(singleMealEvaluation);
    for (auto food_: orderedFoodList_amp) {
        cout << food_.first << " " << food_.second << endl;
    }
    
    
    finish = clock();
    cout << "Time elapsed (sec): " << (double) (finish-start)/CLOCKS_PER_SEC << endl;
    
    return 0;
}







