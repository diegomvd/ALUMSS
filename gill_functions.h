#ifndef CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H
#define CHAPTER2_2ACTIONMODELCPP_GILL_FUNCTIONS_H

#include <vector>
#include <gsl/gsl_rng.h>
using namespace std;

// gets the indexes of all the neighbours of i and puts them in list
void getNeighbours(unsigned int i, unsigned int n, vector<unsigned int> &list);
// gets the indexes of the neighbours in a given state and puts them in state_indexes
void getStateNeighbours(unsigned int i, unsigned int state, unsigned int n, const vector<unsigned int> &landscape, vector<unsigned int> &state_indexes);
// for every natural patch it gets its connected component membership
void getNaturalCluster(unsigned int n, const vector<unsigned int> &landscape, vector<vector<int>> &natural_components);
// calculates the exposure to nature of patch i
double getExposure2Nature(unsigned int i, unsigned int n, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components);
// calculates the recovery probability pezr unit time of each patch
void getRecoveryPropensity(unsigned int n, double Tr, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &recovery_propensity);
// calculates the degradation probability per unit time of each patch
void getDegradationPropensity(unsigned int n, double Td, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &degradation_propensity);
// calculates the action probability per unit time for each patch
void getActionPropensity(unsigned int n, double w, double a, double g, double consumption_deficit,const vector<unsigned int> &landscape, vector<double> &cropping_propensity, vector<double> &restoring_propensity);
// calculates the abandonment probabilty per unit time of each patch
void getAbandonmentPropensity(unsigned int n, double m, double y0, double Ta, const vector<double> &agricultural_production, vector<double> &abandonment_propensity);
// merges all the probabilites per unit time in a single vector with a cummulative sum
void getPropensityVector(unsigned int n, double Tr, double Td, double w, double a, double g, double m, double y0, double Ta, double consumption_deficit, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, const vector<double> &agricultural_production, vector<double> &propensity_vector);
// initializes the landscape vector
void initializeLandscape(unsigned int n, double a0, gsl_rng  *r, vector<unsigned int> &landscape);
// initializes the agricultural production vector
void initializeProduction(unsigned int n, double y0, double phi, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &agricultural_production);
// initializes the population density
double initializePopulation(double cg0, const vector<double> &agricultural_production);
// calculates the agricultural production
void getAgriculturalProduction(unsigned int n, double y0, double phi, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &agricultural_production);
// calculates the population change
double populationEquation(double r0, double cg0, double population, double agricultural_production);
// solver for the population ODE
double RungeKutta4(double r0, double cg0, double dt, double population, double agricultural_production);

#endif
