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
    
    SubjectGlucose subject;
    
    int type = 2;
    double bodyMass = 85.9; // unit: kg
//    double fastGlucose = 7.96; // unit: mmol/L
//    double fastInsulin = 12.84; // unit: mU/L
    double fastGlucose = 6.7; // unit: mmol/L
    double fastInsulin = 7.65; // unit: mU/L
    
    std::vector<std::pair<double, double>> foodIntakeEvents = {
        {0., 0.}, {105., 75E3},  {345., 0.}
    };
    
    std::vector<std::tuple<double, double, double, int>> exerciseEvents = {
        {0., 60., 1., 0}
//        {60., 90., 2., 1},
//        {200, 500, 3, 2}
    };
    
    // set up subject specifics and its daily activities
    std::tuple<int, double, double, double> subjectSpecifics = {type, bodyMass, fastGlucose, fastInsulin};
    subject.SetSubjectSpecifics(subjectSpecifics);
    subject.SetFoodIntakeEvents(foodIntakeEvents);
    subject.SetExerciseEvents(exerciseEvents);
    subject.CombineFoodIntakeExerciseEvents();
    
    
    string param_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/fitted_params/optimized_fitted_params_D2_Karstoft_con_for_IW.dat";
    
//    subject.SetFittedParams_EDES(param_file);
    subject.SetFittedParams_EDES_Ex(param_file_path);
    
    
    // perform evolution
    // Note: for each eaily activities, one can have to perform the 'EvolutionUnderFoodIntakeExerciseEvents' ONCE,
    //  while extracting the results multiple times with different timeIntervals. This is because
    //  after performing the evolution, an internal interpolation function on the resultant glucose levels is stored,
    //  so that to extract glucose levels with different timeIntervals, one can directly start from the interpolation function.
    subject.EvolutionUnderFoodIntakeExerciseEvents();
    
    // extract the obtained glucose levels at the user specified time intervals
    double timeInterval = 5;
    std::vector<std::pair<double, double>> glucoses = subject.GlucoseUnderFoodIntakeExerciseEvents(timeInterval);
    std::vector<std::pair<double, double>> insulins = subject.InsulinUnderFoodIntakeExerciseEvents(timeInterval);
    
    // export results:
    ofstream output;
    string output_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/";
    output.open(output_file_path + "output.dat", ofstream::out);
    for (int i = 0; i < glucoses.size(); ++i){
        cout << glucoses[i].first   << " " << glucoses[i].second << endl;
        output << glucoses[i].first  << " " << glucoses[i].second << " " << insulins[i].second << endl;
    }
    output.close();

    return 0;
}







