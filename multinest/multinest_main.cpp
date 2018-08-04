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

#include "multinest.h"

#include "E_DES_Ex_Med.hpp"
#include "SubjectGlucose.hpp"

using namespace std;

//double E_DES_Ex_Med::k1_H = 0.0163619;
//double E_DES_Ex_Med::k2_H = 0.18896;
//double E_DES_Ex_Med::k3_H = 0.0566132;
//double E_DES_Ex_Med::k4_H = 0.00143478;
//double E_DES_Ex_Med::k5_H = 0.000120797;
//double E_DES_Ex_Med::k6_H = 0.550253;
//double E_DES_Ex_Med::k7_H = 0.00749822;
//double E_DES_Ex_Med::k8_H = 2.78914;
//double E_DES_Ex_Med::k9_H = 0.0261024;
//double E_DES_Ex_Med::k10_H = 0.05;
//double E_DES_Ex_Med::k11_H = 0.;
//double E_DES_Ex_Med::sigma_H = 1.22447;
//double E_DES_Ex_Med::KM_H = 15.532;
//double E_DES_Ex_Med::k4e_H = 0.;        // exercise
//double E_DES_Ex_Med::k5e_H = 10.;        // exercise
//double E_DES_Ex_Med::k6e_H = 0.05;        // exercise
//double E_DES_Ex_Med::k8e_H = 0.5;        // exercise
//double E_DES_Ex_Med::lam_H = 1/120.;        // exercise

void MultiNest_LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context){
    
    // T2Mid -> T2Heavy
    auto sGlucose_PE = *static_cast<SubjectGlucose *>(context);
    
    double k5 = 0.558564373452744272E-02 * pow(10, -2 + 2 *Cube[0]);
    double k6 = 0.355383834174836988E+00 * pow(10, -2 + 2 *Cube[1]);
    double k7 = 0.143221305922955957E-03 * pow(10, -2 + 2 *Cube[2]);
    double k8 = 0.335281284603546714E+01 * pow(10, -2 + 2 *Cube[3]);
    double KM = 0.545867398329898208E+02 + 20 * Cube[4];
    
    map<string, double> fittedParamsNew = {
        {"k1", 0.01443 },
        {"k2", 3.795 },
        {"k3", 0.00069218 },
        {"k4", 0.0020239 },
        {"k5", k5 },
        {"k6", k6 },
        {"k7", k7 },
        {"k8", k8 },
        {"k9", 0.07165 },
        {"sigma", 1.3557 },
        {"KM", KM }
    };
    
    // compute SSR
    double SSR = 0.;
    auto num_data_sets = sGlucose_PE.dailyDataSets.size();
    for (int i = 0; i < num_data_sets; ++i) {
        sGlucose_PE.eDES_ex_med.ClearPreRunOutput();
        sGlucose_PE.SetSubjectSpecifics(sGlucose_PE.dailySubjectSpecifics[i]); // still default parameters according to 'type'
        sGlucose_PE.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
        sGlucose_PE.SetFoodIntakeEvents(sGlucose_PE.dailyFoodIntakeEvents[i]);
        sGlucose_PE.SetExerciseEvents(sGlucose_PE.dailyExerciseEvents[i]);
        sGlucose_PE.SetSaIEvents(sGlucose_PE.dailySAInsulinEvents[i]);
        sGlucose_PE.SetLaIEvents(sGlucose_PE.dailyLAInsulinEvents[i]);
        sGlucose_PE.EvolutionUnderFoodExInsulinEvents();
        auto check_pts = sGlucose_PE.ObtainCheckPtsFromDailyDataSets(sGlucose_PE.dailyDataSets[i]);
        auto glucoses_tmp = sGlucose_PE.ObtainGlucose(check_pts);
        auto insulins_tmp = sGlucose_PE.ObtainInsulin(check_pts);
        SSR += sGlucose_PE.ComputeSSR(sGlucose_PE.dailyDataSets[i], glucoses_tmp, insulins_tmp);
        //        SSR = 100.;
    }
    lnew = - SSR / 2.;
    
    Cube[0] = k5;
    Cube[1] = k6;
    Cube[2] = k7;
    Cube[3] = k8;
    Cube[4] = KM;
    Cube[5] = SSR;
    
//    // T2Light -> T2Mid
//    auto sGlucose_PE = *static_cast<SubjectGlucose *>(context);
//    
//    double k5 = 0.616149969195030692E-02 * pow(10, -2 + 2 *Cube[0]);
//    double k6 = 0.586133719919008533E+00 * pow(10, -2 + 2 *Cube[1]);
//    double k7 = 0.131910649054100878E-01 * pow(10, -2 + 2 *Cube[2]);
//    double k8 = 0.335281340038196740E+01 * pow(10, -2 + 2 *Cube[3]);
//    double KM = 0.348813785020632423E+02 + 20 * Cube[4];
//    
//    map<string, double> fittedParamsNew = {
//        {"k1", 0.01443 },
//        {"k2", 3.795 },
//        {"k3", 0.00069218 },
//        {"k4", 0.0020239 },
//        {"k5", k5 },
//        {"k6", k6 },
//        {"k7", k7 },
//        {"k8", k8 },
//        {"k9", 0.07165 },
//        {"sigma", 1.3557 },
//        {"KM", KM }
//    };
//    
//    // compute SSR
//    double SSR = 0.;
//    auto num_data_sets = sGlucose_PE.dailyDataSets.size();
//    for (int i = 0; i < num_data_sets; ++i) {
//        sGlucose_PE.eDES_ex_med.ClearPreRunOutput();
//        sGlucose_PE.SetSubjectSpecifics(sGlucose_PE.dailySubjectSpecifics[i]); // still default parameters according to 'type'
//        sGlucose_PE.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
//        sGlucose_PE.SetFoodIntakeEvents(sGlucose_PE.dailyFoodIntakeEvents[i]);
//        sGlucose_PE.SetExerciseEvents(sGlucose_PE.dailyExerciseEvents[i]);
//        sGlucose_PE.SetSaIEvents(sGlucose_PE.dailySAInsulinEvents[i]);
//        sGlucose_PE.SetLaIEvents(sGlucose_PE.dailyLAInsulinEvents[i]);
//        sGlucose_PE.EvolutionUnderFoodExInsulinEvents();
//        auto check_pts = sGlucose_PE.ObtainCheckPtsFromDailyDataSets(sGlucose_PE.dailyDataSets[i]);
//        auto glucoses_tmp = sGlucose_PE.ObtainGlucose(check_pts);
//        auto insulins_tmp = sGlucose_PE.ObtainInsulin(check_pts);
//        SSR += sGlucose_PE.ComputeSSR(sGlucose_PE.dailyDataSets[i], glucoses_tmp, insulins_tmp);
//        //        SSR = 100.;
//    }
//    lnew = - SSR / 2.;
//    
//    Cube[0] = k5;
//    Cube[1] = k6;
//    Cube[2] = k7;
//    Cube[3] = k8;
//    Cube[4] = KM;
//    Cube[5] = SSR;
    
    
//    // GFT -> T2Light
//    auto sGlucose_PE = *static_cast<SubjectGlucose *>(context);
//    
//    double k5 = 0.0062466 * pow(10, -2 + 2 *Cube[0]);
//    double k6 = 0.58634 * pow(10, -2 + 2 *Cube[1]);
//    double k7 = 0.01762 * pow(10, -2 + 2 *Cube[2]);
//    double k8 = 3.4042 * pow(10, -2 + 2 *Cube[3]);
//    double KM = 32.532 + 20 * Cube[4];
//    
//    map<string, double> fittedParamsNew = {
//        {"k1", 0.01443 },
//        {"k2", 3.795 },
//        {"k3", 0.00069218 },
//        {"k4", 0.0020239 },
//        {"k5", k5 },
//        {"k6", k6 },
//        {"k7", k7 },
//        {"k8", k8 },
//        {"k9", 0.07165 },
//        {"sigma", 1.3557 },
//        {"KM", KM }
//    };
//    
//    // compute SSR
//    double SSR = 0.;
//    auto num_data_sets = sGlucose_PE.dailyDataSets.size();
//    for (int i = 0; i < num_data_sets; ++i) {
//        sGlucose_PE.eDES_ex_med.ClearPreRunOutput();
//        sGlucose_PE.SetSubjectSpecifics(sGlucose_PE.dailySubjectSpecifics[i]); // still default parameters according to 'type'
//        sGlucose_PE.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
//        sGlucose_PE.SetFoodIntakeEvents(sGlucose_PE.dailyFoodIntakeEvents[i]);
//        sGlucose_PE.SetExerciseEvents(sGlucose_PE.dailyExerciseEvents[i]);
//        sGlucose_PE.SetSaIEvents(sGlucose_PE.dailySAInsulinEvents[i]);
//        sGlucose_PE.SetLaIEvents(sGlucose_PE.dailyLAInsulinEvents[i]);
//        sGlucose_PE.EvolutionUnderFoodExInsulinEvents();
//        auto check_pts = sGlucose_PE.ObtainCheckPtsFromDailyDataSets(sGlucose_PE.dailyDataSets[i]);
//        auto glucoses_tmp = sGlucose_PE.ObtainGlucose(check_pts);
//        auto insulins_tmp = sGlucose_PE.ObtainInsulin(check_pts);
//        SSR += sGlucose_PE.ComputeSSR(sGlucose_PE.dailyDataSets[i], glucoses_tmp, insulins_tmp);
//        //        SSR = 100.;
//    }
//    lnew = - SSR / 2.;
//    
//    Cube[0] = k5;
//    Cube[1] = k6;
//    Cube[2] = k7;
//    Cube[3] = k8;
//    Cube[4] = KM;
//    Cube[5] = SSR;
    
//    // health -> GFT
//    auto sGlucose_PE = *static_cast<SubjectGlucose *>(context);
//    
//    double k5 = 0.00706987 * pow(10, -2 + 2 *Cube[0]);
//    double k6 = 1.1352 * pow(10, -2 + 2 *Cube[1]);
//    double k7 = 0.4393 * pow(10, -2 + 2 *Cube[2]);
//    double k8 = 3.4121 * pow(10, -2 + 2 *Cube[3]);
//    double KM = 13.157153 + 20 * Cube[4];
//    
//    map<string, double> fittedParamsNew = {
//        {"k1", 0.01443 },
//        {"k2", 3.795 },
//        {"k3", 0.00069218 },
//        {"k4", 0.0020239 },
//        {"k5", k5 },
//        {"k6", k6 },
//        {"k7", k7 },
//        {"k8", k8 },
//        {"k9", 0.07165 },
//        {"sigma", 1.3557 },
//        {"KM", KM }
//    };
//    
//    // compute SSR
//    double SSR = 0.;
//    auto num_data_sets = sGlucose_PE.dailyDataSets.size();
//    for (int i = 0; i < num_data_sets; ++i) {
//        sGlucose_PE.eDES_ex_med.ClearPreRunOutput();
//        sGlucose_PE.SetSubjectSpecifics(sGlucose_PE.dailySubjectSpecifics[i]); // still default parameters according to 'type'
//        sGlucose_PE.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
//        sGlucose_PE.SetFoodIntakeEvents(sGlucose_PE.dailyFoodIntakeEvents[i]);
//        sGlucose_PE.SetExerciseEvents(sGlucose_PE.dailyExerciseEvents[i]);
//        sGlucose_PE.SetSaIEvents(sGlucose_PE.dailySAInsulinEvents[i]);
//        sGlucose_PE.SetLaIEvents(sGlucose_PE.dailyLAInsulinEvents[i]);
//        sGlucose_PE.EvolutionUnderFoodExInsulinEvents();
//        auto check_pts = sGlucose_PE.ObtainCheckPtsFromDailyDataSets(sGlucose_PE.dailyDataSets[i]);
//        auto glucoses_tmp = sGlucose_PE.ObtainGlucose(check_pts);
//        auto insulins_tmp = sGlucose_PE.ObtainInsulin(check_pts);
//        SSR += sGlucose_PE.ComputeSSR(sGlucose_PE.dailyDataSets[i], glucoses_tmp, insulins_tmp);
//        //        SSR = 100.;
//    }
//    lnew = - SSR / 2.;
//    
//    Cube[0] = k5;
//    Cube[1] = k6;
//    Cube[2] = k7;
//    Cube[3] = k8;
//    Cube[4] = KM;
//    Cube[5] = SSR;

    
//    auto sGlucose_PE = *static_cast<SubjectGlucose *>(context);
//    
//    double k1 = 0.0163619 * pow(10, -2 + 4 *Cube[0]);
//    double k2 = 0.18896 * pow(10, -2 + 4 *Cube[1]);
//    double k3 = 0.0566132 * pow(10, -2 + 4 *Cube[2]);
//    double k4 = 0.00143478 * pow(10, -2 + 4 *Cube[3]);
//    double k5 = 0.000120797 * pow(10, -2 + 4 *Cube[4]);
//    double k6 = 0.550253 * pow(10, -2 + 4 *Cube[5]);
//    double k7 = 0.00749822 * pow(10, -2 + 4 *Cube[6]);
//    double k8 = 2.78914 * pow(10, -2 + 4 *Cube[7]);
//    double k9 = 0.0261024 * pow(10, -2 + 4 *Cube[8]);
////    double sigma = 1.22447 * pow(10, -2 + 4 *Cube[9]);
////    double KM = 15.532 * pow(10, -2 + 4 *Cube[10]);
//    double sigma = 1. + Cube[9];
//    double KM = 5. + 40 * Cube[10];
//    
//    map<string, double> fittedParamsNew = {
//        {"k1", k1 },
//        {"k2", k2 },
//        {"k3", k3 },
//        {"k4", k4 },
//        {"k5", k5 },
//        {"k6", k6 },
//        {"k7", k7 },
//        {"k8", k8 },
//        {"k9", k9 },
////        {"k10", 0.05 * pow(10, -2 + 4 *Cube[9]) },
////        {"sigma", 1. + Cube[10]},
////        {"KM", 5. + 25. * Cube[11] }
//        {"sigma", sigma },
//        {"KM", KM }
//    };
//    
//    // compute SSR
//    double SSR = 0.;
//    auto num_data_sets = sGlucose_PE.dailyDataSets.size();
//    for (int i = 0; i < num_data_sets; ++i) {
//        sGlucose_PE.eDES_ex_med.ClearPreRunOutput();
//        sGlucose_PE.SetSubjectSpecifics(sGlucose_PE.dailySubjectSpecifics[i]); // still default parameters according to 'type'
//        sGlucose_PE.SetFittedParams_EDES_Ex_Med(fittedParamsNew); // update the minimized fitted-params
//        sGlucose_PE.SetFoodIntakeEvents(sGlucose_PE.dailyFoodIntakeEvents[i]);
//        sGlucose_PE.SetExerciseEvents(sGlucose_PE.dailyExerciseEvents[i]);
//        sGlucose_PE.SetSaIEvents(sGlucose_PE.dailySAInsulinEvents[i]);
//        sGlucose_PE.SetLaIEvents(sGlucose_PE.dailyLAInsulinEvents[i]);
//        sGlucose_PE.EvolutionUnderFoodExInsulinEvents();
//        auto check_pts = sGlucose_PE.ObtainCheckPtsFromDailyDataSets(sGlucose_PE.dailyDataSets[i]);
//        auto glucoses_tmp = sGlucose_PE.ObtainGlucose(check_pts);
//        auto insulins_tmp = sGlucose_PE.ObtainInsulin(check_pts);
//        SSR += sGlucose_PE.ComputeSSR(sGlucose_PE.dailyDataSets[i], glucoses_tmp, insulins_tmp);
//        //        SSR = 100.;
//    }
//    lnew = - SSR / 2.;
//    
//    Cube[0] = k1;
//    Cube[1] = k2;
//    Cube[2] = k3;
//    Cube[3] = k4;
//    Cube[4] = k5;
//    Cube[5] = k6;
//    Cube[6] = k7;
//    Cube[7] = k8;
//    Cube[8] = k9;
//    Cube[9] = sigma;
//    Cube[10] = KM;
//    Cube[11] = SSR;

}


void MultiNest_Dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context) {

    // convert the 2D Fortran arrays to C++ arrays
    
    
    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    
    int i, j;
    
    double postdist[nSamples][nPar + 2];
    for( i = 0; i < nPar + 2; i++ )
        for( j = 0; j < nSamples; j++ )
            postdist[j][i] = posterior[0][i * nSamples + j];
    
    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    
    double pLivePts[nlive][nPar + 1];
    for( i = 0; i < nPar + 1; i++ )
        for( j = 0; j < nlive; j++ )
            pLivePts[j][i] = physLive[0][i * nlive + j];

}


int main(int argc, const char * argv[]) {
    
    SubjectGlucose sGlucose;
    vector<string> dataset_files;
    
    // +++++++++++++++++   report_data: T2Heavy
    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/T2Heavy_data.dat";
    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
        {2 , 83.3, 9.15, 8., 9.15, 8.}
    };
    
//    // +++++++++++++++++   report_data: T2Mid
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/T2Mid_data.dat";
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 8., 8., 8., 8.}
//    };

    
//    // +++++++++++++++++   report_data: T2Light
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/T2Light_data.dat";
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 6.84, 8., 6.84, 8.}
//    };
    
//    // +++++++++++++++++   report_data: GFT
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/GFT_data.dat";
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {2 , 83.3, 6.4, 8., 6.4, 8.}
//    };
    
//    // +++++++++++++++++   report_data: health
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/report_data/health_data.dat";
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {0 , 75., 5., 8., 5., 8.}
//    };
    
//    string dataset_file_path = "/Users/Jue/Desktop/precision_health/models/E_DES/E_DES/datasets_for_param_estimation/Mass_thesis/fig14_h_glucose_insulin_with_error.dat"; // Using the fitted curve in Figure 14 of Dr. Mass's thesis
//    dataset_files = {dataset_file_path};
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {0 , 75., 4.55, 16., 4.55, 16.}
//    };
//    
    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
        { {0., 75E3}, {360., 0.} }
    };
    
//    string fileDir = "/Users/Jue/Desktop/precision_health/data/E_DES/data_parameter_estimation/clean_data/";
//    vector<string> fileNames= {"34_female", "34_male", "35", "36", "37", "38", "39", "40", "41",
//        "42", "43", "44"};
//    for (int i = 0; i < fileNames.size(); ++i) {
//        string tmp = fileDir + fileNames[i] + "_combine_with_error.dat";
//        dataset_files.push_back(tmp);
//    }
//    
//    vector<SubjectGlucose::TypeSubjectSpecifics> dailySubjectSpecifics = {
//        {{0, 75., 4.91, 7.72, 4.91, 7.72}, {0, 75., 5.03, 8.36, 5.03,
//            8.36}, {0, 75., 4.77, 11.86, 4.77, 11.86}, {0, 75., 4.67, 6.32,
//                4.67, 6.32}, {0, 75., 4.7, 12.48, 4.7, 12.48}, {0, 75., 5.38, 5.69,
//                    5.38, 5.69}, {0, 75., 4.81, 4.54, 4.81, 4.54}, {0, 75., 5.03, 4.76,
//                        5.03, 4.76}, {0, 75., 4.85, 3.84, 4.85, 3.84}, {0, 75., 4.98, 4.34,
//                            4.98, 4.34}, {0, 75., 4.69, 3.77, 4.69, 3.77}, {0, 75., 4.88, 5.54,
//                                4.88, 5.54}}
//    };
//    
//    vector<SubjectGlucose::TypeDailyFoodIntakeEvents> dailyFoodIntakeEvents = {
//        {{{0., 75000}, {180., 0.}}, {{0., 75000}, {180., 0.}}, {{0.,
//            75000}, {120., 0.}}, {{0., 75000}, {240., 0.}}, {{0.,
//                75000}, {120., 0.}}, {{0., 75000}, {120., 0.}}, {{0.,
//                    75000}, {120., 0.}}, {{0., 75000}, {120., 0.}}, {{0.,
//                        75000}, {120., 0.}}, {{0., 75000}, {120., 0.}}, {{0.,
//                            75000}, {120., 0.}}, {{0., 75000}, {120., 0.}}}
//    };
    
//    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {},{},{},{},{},{},{},{},{},{},{},{} };
//    
//    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {},{},{},{},{},{},{},{},{},{},{},{} };
//    
//    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {},{},{},{},{},{},{},{},{},{},{},{} } ;
    
    dataset_files = {dataset_file_path};
    
    vector<SubjectGlucose::TypeDailyExerciseEvents> dailyExerciseEvents = { {} };
    
    vector<SubjectGlucose::TypeDailySAInsulinEvents> dailySAInsulinEvents = { {} };
    
    vector<SubjectGlucose::TypeDailyLAInsulinEvents> dailyLAInsulinEvents = { {} } ;
    
    sGlucose.SetDatasetsForParameterEstimation(dataset_files, dailySubjectSpecifics, dailyFoodIntakeEvents, dailyExerciseEvents, dailySAInsulinEvents, dailyLAInsulinEvents);

    
    // set the MultiNest sampling parameters
    
    int IS = 0;					// do Nested Importance Sampling?
    
    int mmodal = 1;					// do mode separation?
    
    int ceff = 0;					// run in constant efficiency mode?
    
    int nlive = 5000;				// number of live points
    
    double efr = 0.8;				// set the required efficiency
    
    double tol = 0.5;				// tol, defines the stopping criteria
    
    int ndims = 5;					// dimensionality (no. of free parameters)
    
    int nPar = 6;					// total no. of parameters including free & derived parameters
    
    int nClsPar = 5;				// no. of parameters to do mode separation on
    
    int updInt = 5000;				// after how many iterations feedback is required & the output files should be updated
    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
    
    double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
    
    int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
    
    int pWrap[ndims];				// which parameters to have periodic boundary conditions?
    //for(int i = 0; i < ndims; i++) pWrap[i] = 0;
    
    char root[100] = "chains/E_DES_T2Heavy_";			// root for output files
    
    int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
    
    int fb = 1;					// need feedback on standard output?
    
    int resume = 0;					// resume from a previous job?
    
    int outfile = 1;				// write output files?
    
    int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
    // set it to F if you want your main program to handle MPI initialization
    
    double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
    
    int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
    // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
    
    void *context = 0;				// not required by MultiNest, any additional information user wants to pass
    
    context = &sGlucose;
    
    // calling MultiNest
    nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, MultiNest_LogLike, MultiNest_Dumper, context);

    
    
    return 0;
}







