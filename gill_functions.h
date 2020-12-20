#ifndef CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H
#define CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <fstream>

using namespace std;

/*
1-Helper functions
2-Calculation of Ecosystem Service provision
3-Calculation of events' propensities
4-Initialization functions
5-ODEs and solver
*/


///////////////////////////////////////////////////////////////////////////////
// 1- Helper functions:
//       - getNeighbourMatrix
//       - getNeighbours
//       - getNeighboursState
///////////////////////////////////////////////////////////////////////////////

void getNeighbourMatrix(vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, unsigned int d);
void getNeighbours(vector<unsigned int> &neighboursList, const vector<vector<unsigned int>> &neighbourMatrix, unsigned int i);
void getNeighboursState(vector<unsigned int> &neighboursState, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, unsigned int i, unsigned int state);

////////////////////////////////////////////////////////////////////////////////
// 2- Calculation of Ecosystem Service provision:
//       - getNaturalConnectedComponents
//       - getEcosystemServiceProvision
//       - getAgriculturalProduction
////////////////////////////////////////////////////////////////////////////////

void getNaturalConnectedComponents(vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape);
void updateNCCadding(vector<vector<int>> &naturalComponents, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, unsigned int i);
void updateNCCremoving(vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape, unsigned int i);
void getEcosystemServiceProvision(vector<double> &ecosystemServices, const vector<vector<int>> &naturalComponents, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, double sar);
void getAgriculturalProduction(vector<double> &agriculturalProduction, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double ksi);

////////////////////////////////////////////////////////////////////////////////
// 3- Calculation of events' propensities:
//       - getConsumptionDeficit
//       - getSpontaneousPropensity
//       - getAgroPropensity
//       - getAbandonmentPropensity
//       - getPropensityVector
////////////////////////////////////////////////////////////////////////////////

double getConsumptionDeficit(const vector<double> &agriculturalProduction, const vector<double> &population);
void getSpontaneousPropensity(vector<double> &recoveryPropensity, vector<double> &degradationPropensity, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double Tr, double Td);
void getAgroPropensity(vector<double> &expansionPropensity, vector<double> &intensePropensity, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, const vector<double> &agriculturalProduction, const vector<double> &population, double w, double a, double Tag);
void getAbandonmentPropensity(vector<double> &naturalAbandonPropensity, vector<double> &degradedAbandonPropensity, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double Tab);
void getPropensityVector(vector<double> &propensityVector, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, const vector<double> &agriculturalProduction, const vector<double> &population, double Tr, double Td, double w, double a, double Tag, double Tab);
////////////////////////////////////////////////////////////////////////////////
// 4- Initialization functions:
//       - initializeLandscape
//       - initializePopulation
//       - initializeSES
////////////////////////////////////////////////////////////////////////////////

void initializeLandscape( vector<unsigned int> &landscape, vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, double ao0, double ai0, double w, gsl_rng  *r);
void initializePopulation( vector<double> &population, const vector<double> &agriculturalProduction);
void initializeSES( vector<unsigned int> &landscape, vector<double> &population, vector<vector<int>> &naturalComponents, vector<double> &agriculturalProduction, vector<double> &ecosystemServices, vector<vector<unsigned int>> &neighbourMatrix, vector<vector<unsigned int>> &neighbourMatrixES, unsigned int n, double ao0, double ai0, double ksi, double sar, unsigned int d, double w,gsl_rng  *r);

////////////////////////////////////////////////////////////////////////////////
// 5- ODEs and solver:
//       - populationEquation
//       - rungeKutta4
////////////////////////////////////////////////////////////////////////////////

double populationEquation(double population, double agriculturalProduction);
void rungeKutta4(vector<double> &population, vector<double> &agriculturalProduction, double dt);

////////////////////////////////////////////////////////////////////////////////
// 6- Outputs:
//       - saveAggregated
//       - saveLandscape
//       - saveComponents
//       - saveLandMetrics
//       - saveRipley
////////////////////////////////////////////////////////////////////////////////

void saveAggregated(ofstream &file, double t, vector<double> &population, vector<unsigned int> &landscape);
void saveLandscape(ofstream &file, double t, vector<unsigned int> &landscape);
void saveComponents(ofstream &file, double t, vector<unsigned int> &landscape, vector<vector<int>> &naturalComponents);
void saveLandMetrics(ofstream &file, double t, vector<vector<int>> &naturalComponents, vector<double> &ecosystemServices);
void saveRipley(ofstream &file, double t, unsigned int n, vector<unsigned int> &landscape, double ripleyDistance);

#endif
