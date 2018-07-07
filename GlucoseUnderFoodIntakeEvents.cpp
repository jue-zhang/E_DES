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
    
    E_DES eDES;
    double bodyMass = 83.3; // unit: kg
    double Gpl_init = 6.; // unit: mmol/L
    double Ipl_init = 8.; // unit: mU/L
    
//    std::vector<std::pair<double, double>> foodIntakeEvents = {
//        {8*60., 40E3},     // breakfast, eg., 8 am
//        {12*60., 40E3},   // lunch, eg., 12 am
//        {17*60., 40E3},   // dinner, eg., 5 pm
//        {24*60., 20E3},      // midnight snack, eg. 12 pm
//        {32*60., 0}       // the next morning at 8 am
//    };
    
    std::vector<std::pair<double, double>> foodIntakeEvents = {
        {7*60., 75E3},     // breakfast, eg., 7 am
        {13*60., 75E3},   // lunch, eg., 1 pm
        {19*60., 75E3},   // dinner, eg., 7 pm
        {23*60., 20E3},      // midnight snack, eg. 11 pm
        {31*60., 0}       // the next morning at 7 am
    };
    
//    // Set the subject type and fitted parameters
//    // 1. Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
//    // 2. Note: Type-I is NOT supported yet, and Type-II parameters are used instead.
//    eDES.SetSubjectTypeFittedParams(2);
    
    // Load fitted-params-D2
    ifstream fitted_params_D2;
    fitted_params_D2.open("/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2.dat", ofstream::in);
    eDES.LoadFittedParams(fitted_params_D2);
    
    std::vector<std::pair<double, double>> glucose_at_10_mins;
    ofstream output;
    string output_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/";
    output.open(output_file_path + "output.dat", ofstream::out);
    glucose_at_10_mins = eDES.GlucoseUnderFoodIntakeEvents(bodyMass, Gpl_init, Ipl_init, foodIntakeEvents);
    
    for (auto glucose_at_10_min: glucose_at_10_mins) {
        cout << glucose_at_10_min.first / 60 << " " << glucose_at_10_min.second << endl;
        output << glucose_at_10_min.first / 60 << " " << glucose_at_10_min.second << endl;
    }

    return 0;
}







