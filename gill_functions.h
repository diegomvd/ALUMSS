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
void getNaturalCluster(vector<vector<int>> &natural_components, unsigned int n, const vector<unsigned int> &landscape);
// calculates the exposure to nature of patch i
double getExposure2Nature(unsigned int i, unsigned int n, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components);
// calculates the recovery probability pezr unit time of each patch
void getRecoveryPropensity(unsigned int n, double Tr, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &recovery_propensity);
// calculates the degradation probability per unit time of each patch
void getDegradationPropensity(unsigned int n, double Td, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &degradation_propensity);
// calculates the action probability per unit time for each patch
void getActionPropensity(unsigned int n, double w, double a, double g, double consumption_deficit,const vector<unsigned int> &landscape, vector<double> &cropping_propensity, vector<double> &intense_propensity);
// calculates the abandonment probabilty per unit time of each patch
void getAbandonmentPropensity(double Ta, const vector<unsigned int> &landscape, const vector<double> &agricultural_production, const vector<double> &maintenance_costs, vector<double> &abandonmentO_propensity, vector<double> &abandonmentI_propensity);
// merges all the probabilites per unit time in a single vector with a cummulative sum
void getPropensityVector(unsigned int n, double Tr, double Td, double w, double a, double g, double Ta, double consumption_deficit, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, const vector<double> &agricultural_production, const vector<double> &maintenance_costs, vector<double> &propensity_vector);
// initializes the landscape vector
void initializeLandscape(vector<unsigned int> &landscape, unsigned int n, double ao0, double ai0, gsl_rng  *r);
// initializes the agricultural production vector
void initializeProduction(unsigned int n, double y0, double phi, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &agricultural_production);
// initializes the population density
double initializePopulation(double cg0, const vector<double> &agricultural_production);
// initializes vector storing maintenance costs
void initializeMaintenanceCosts(vector<double> &maintenance_costs, double y0, double m, const vector<unsigned int> &landscape);
// calculates the agricultural production
void getAgriculturalProduction(unsigned int n, double y0, double phi, double k, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &agricultural_production);
// increments maintenance cost of intense patches
void updateMaintenanceIntense(vector<double> &maintenance_costs, double dtp, double y0, double m, double k, double Ti, const vector<unsigned int> &landscape);
// calculates the population change
double populationEquation(double r0, double agricultural_production, double population, double consumption);
// calculates the per capita consumption change
double consumptionEquation(double kg, double kd, double minimum_consumption, double agricultural_production, double population, double consumption);
// solver for the population-consumption ODE system
void RungeKutta4(double r0, double kg, double kd, double minimum_consumption, double dt, double agricultural_production, vector<double> &population, vector<double> &consumption);
#endif
