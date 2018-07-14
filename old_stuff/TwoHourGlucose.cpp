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
    double foodIntake = 75E3; // unit: mg
    double bodyMass = 75.; // unit: kg
    double Gpl_init = 5.; // unit: mmol/L
    double Ipl_init = 8.; // unit: mU/L
    cout << "two hour glucose: " << eDES.TwoHourGlucose(foodIntake, bodyMass, Gpl_init, Ipl_init) << endl;
    
    
    return 0;
}







