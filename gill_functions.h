#ifndef CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H
#define CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H

#include <vector>
#include <gsl/gsl_rng.h>
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
//       - getNeighbours
//       - getNeighboursState
///////////////////////////////////////////////////////////////////////////////

void getNeighbours(vector<unsigned int> &neighboursList, unsigned int i, unsigned int n);
void getNeighboursState(vector<unsigned int> &neighboursState, const vector<unsigned int> &landscape, unsigned int i, unsigned int state, unsigned int n);

////////////////////////////////////////////////////////////////////////////////
// 2- Calculation of Ecosystem Service provision:
//       - getNaturalConnectedComponents
//       - getEcosystemServiceProvision
//       - getAgriculturalProduction
////////////////////////////////////////////////////////////////////////////////

void getNaturalConnectedComponents(vector<vector<int>> &naturalComponents, unsigned int n, const vector<unsigned int> &landscape);
double getEcosystemServiceProvision(const vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape, unsigned int i, unsigned int n, double ess);
void getAgriculturalProduction(vector<double> &agriculturalProduction, const vector<unsigned int> &landscape, const vector<vector<int>> &naturalComponents, unsigned int n, double ys0, double yn0, double pSD, double ori, double ini, double ess, gsl_rng  *r);
////////////////////////////////////////////////////////////////////////////////
// 3- Calculation of events' propensities:
//       - getConsumptionDeficit
//       - getSpontaneousPropensity
//       - getAgroPropensity
//       - getAbandonmentPropensity
//       - getPropensityVector
////////////////////////////////////////////////////////////////////////////////

double getConsumptionDeficit(const vector<double> &agriculturalProduction, const vector<double> &population, double consumption);
void getSpontaneousPropensity(vector<double> &recoveryPropensity, vector<double> &degradationPropensity, const vector<unsigned int> &landscape, const vector<vector<int>> &naturalComponents, unsigned int n, double Tr, double Td, double ess);
void getAgroPropensity(vector<double> &expansionPropensity, vector<double> &intensePropensity, const vector<unsigned int> &landscape, const vector<double> &agriculturalProduction, const vector<double> &population, double consumption, unsigned int n, double w, double a, double g);
void getAbandonmentPropensity(vector<double> &naturalAbandonPropensity, vector<double> &degradedAbandonPropensity, const vector<unsigned int> &landscape, double Ta, double ori, double ini, double d, double b);
void getPropensityVector(vector<double> &propensityVector, const vector<unsigned int> &landscape, const vector<vector<int>> &naturalComponents, const vector<double> &agriculturalProduction, const vector<double> &population, double consumption, unsigned int n, double Tr, double Td, double w, double a, double g, double Ta, double ori, double ini, double ess, double d, double b);

////////////////////////////////////////////////////////////////////////////////
// 4- Initialization functions:
//       - initializeLandscape
//       - initializePopulation
//       - initializeSES
////////////////////////////////////////////////////////////////////////////////

void initializeLandscape( vector<unsigned int> &landscape, unsigned int n, double ao0, double ai0, gsl_rng  *r);
void initializePopulation( vector<double> &population, double consumption, const vector<double> &agriculturalProduction);
void initializeSES( vector<unsigned int> &landscape, vector<double> &population, vector<vector<int>> &naturalComponents, vector<double> &agriculturalProduction, double consumption, unsigned int n, double ao0, double ai0, gsl_rng  *r, double ys0, double yn0, double pSD, double ori, double ini, double ess);

////////////////////////////////////////////////////////////////////////////////
// 5- ODEs and solver:
//       - populationEquation
//       - rungeKutta4
////////////////////////////////////////////////////////////////////////////////

double populationEquation(double population, double consumption, double agriculturalProduction, double r0);
void rungeKutta4(vector<double> &population, double consumption, vector<double> &agriculturalProduction, double dt, double r0);
#endif
