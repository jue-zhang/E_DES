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
    double foodIntake = 30E3; // unit: mg
    double bodyMass = 83.3; // unit: kg
    double Gpl_init = 5.; // unit: mmol/L
    double Ipl_init = 8.; // unit: mU/L
    
    // Set the subject type and fitted parameters
    // 1. Type of subject: 0 - healthy person, 1 - Type-I diabetes, 2 - Type-II diabetes
    // 2. Note: Type-I is NOT supported yet, and Type-II parameters are used instead.
    eDES.SetSubjectTypeFittedParams(2);
    
    // Note:
    // 1. Special treatment for the healty person: Currently, the model (without manually modifications) predicts glucose values
    //      below the initial values after about three hours of food intake. This problem is more severe when the amount of
    //      foodIntake is large. For example, when foodIntake = 1000E3 mg, the final glucose can be around 2~3, which is not
    //      acceptable. To overcome this, I manually demand that as long as the glucose level is below the initial value, the
    //      initial value will be used for the following time intervals. Therefore, it looks like that the glucose level
    //      returns to the initial value eventually.
    // 2. Diabetes: In "EightHourGlucose", the above special treatment is removed.
    std::vector<double> glucoses = eDES.EightHourGlucose(foodIntake, bodyMass, Gpl_init, Ipl_init);
    cout << "eight hour glucose (10-min intervals): " << endl;
    int i = 1;
    for (auto glucose: glucoses) {
        cout << i*10 << " " << glucose << endl;
        ++i;
    }
    
    return 0;
}







