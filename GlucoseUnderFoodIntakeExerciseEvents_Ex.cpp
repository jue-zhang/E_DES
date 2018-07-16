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
    double bodyMass = 67.0; // unit: kg
    double fastGlucose = 6.8; // unit: mmol/L
    double fastInsulin = 8.; // unit: mU/L
    
    // Food intake events:
    // format: <t_eat, foodIntake>
//    std::vector<std::pair<double, double>> foodIntakeEvents = {
//        {7*60., 75E3},     // breakfast, eg., 7 am
//        {13*60., 75E3},   // lunch, eg., 1 pm
//        {19*60., 75E3},   // dinner, eg., 7 pm
//        {23*60., 20E3},      // midnight snack, eg. 11 pm
//        {31*60., 0}       // the next morning at 7 am
//    };
    std::vector<std::pair<double, double>> foodIntakeEvents = {
        {0., 0},     // breakfast, eg., 7 am
        {60., 27E3},   // lunch, eg., 1 pm
        {300., 45E3},   // dinner, eg., 7 pm
        {660., 42E3},      // midnight snack, eg. 11 pm
        {1440., 0}       // the next morning at 7 am
    };
    
    // Exercise events:
    // format: <t_start, t_finish, ex_intensity, ex_type>
    // NOTE: 1. EXERCISE TYPE IS NOT SUPPORTED YET!
    //       2. Currently, to obtain a reasonable result, please take
    //          only one exercise a day, 'ex_intensity' < 1., exercise duration < 1 hour
    std::vector<std::tuple<double, double, double, int>> exerciseEvents = {
//        {9.*60, 10.*60., 1., 0}
        {150, 165, 1., 0}    // exercise b/n 1pm and 2pm with 'intensity' = 1., 'type' = 0
//        {20.*60, 21.*60, 1., 0}
    };
    
    // set up subject specifics and its daily activities
    std::tuple<int, double, double, double> subjectSpecifics = {type, bodyMass, fastGlucose, fastInsulin};
    subject.SetSubjectSpecifics(subjectSpecifics);
    subject.SetFoodIntakeEvents(foodIntakeEvents);
    subject.SetExerciseEvents(exerciseEvents);
    
    // perform evolution
    // Note: for each eaily activities, one can have to perform the 'EvolutionUnderFoodIntakeExerciseEvents' ONCE,
    //  while extracting the results multiple times with different timeIntervals. This is because
    //  after performing the evolution, an internal interpolation function on the resultant glucose levels is stored,
    //  so that to extract glucose levels with different timeIntervals, one can directly start from the interpolation function.
    subject.EvolutionUnderFoodIntakeExerciseEvents();
    
    // extract the obtained glucose levels at the user specified time intervals
    double timeInterval = 10.;
    std::vector<std::pair<double, double>> glucoses = subject.GlucoseUnderFoodIntakeExerciseEvents(timeInterval);
    
    // export results:
    ofstream output;
    string output_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/output/";
    output.open(output_file_path + "output.dat", ofstream::out);
    for (auto glucose: glucoses) {
        cout << glucose.first   << " " << glucose.second << endl;
        output << glucose.first /60. << " " << glucose.second << endl;
    }
    output.close();

    return 0;
}







