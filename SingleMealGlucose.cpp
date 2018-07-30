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
    double glu_pre = 8.;
    double glu_post = 12.;
    double glu_post_time = 120.;
    double food = 75E3;
    double timeInterval = 5.;
    auto res = sGlucose.ReconGlucoseCurve(type, bodyMass, glu_pre, glu_post, glu_post_time, food, timeInterval);
    
    auto SSR = get<0>(res);
    auto fitted_params = get<1>(res);
    auto glucoses = get<2>(res);
    
    
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
    
    cout << "SSR: " << SSR << endl;
    
    
    finish = clock();
    cout << "Time elapsed (sec): " << (double) (finish-start)/CLOCKS_PER_SEC << endl;
    
    return 0;
}







