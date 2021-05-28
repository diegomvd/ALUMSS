// to compile: c++ gill_main.cpp gill_functions.cpp cokus3.c -o gillespie-ses -lgsl -lgslcblas -lm -Wall -Weffc++ --std=c++17 -lstdc++fs

#include "gill_functions.h"

// #include <sodium.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iostream> //Allows cin/cout
#include <sstream>
#include <stdlib.h> //Allows DOS Commands (I only use "CLS" to clear the screen)
#include <iterator>
#include <filesystem>
#include <ctime>
#include <math.h>
#include <chrono>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

namespace fs = std::filesystem;

extern double ranMT(void);
extern void seedMT(unsigned long int);
extern void seedMT2(void);

#define PI 3.14159265358979323846
#define LOAD_CONF 0
#define SEED_MT 1

///////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM
///////////////////////////////////////////////////////////////////////////////

// time at beginning
auto start = chrono::high_resolution_clock::now();
auto start2 = chrono::high_resolution_clock::now();

int main(int argc, const char * argv[]){

  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETER DECLARATION
  ///////////////////////////////////////////////////////////////////////////////

  int n;  // lenght of the sides of the square landscape: number of cells=n*n
  double SimTime; // total simulation time
  double dtp; // timestep for population dynamics

  double a0; //number of agricultural patches at beggining
  double d0; // number of degraded patches at beggining

  double ksi; // productivity of intense agriculture per es productivity
  double y0; // baseline productivity of low intense agri per es productivity
  double sar; // ecosystem service saturation exponent
  double a; // intensification probability
  double w; // agricultural clustering parameter
  double Tag; // action probability per unit time per unit of consumption deficit
  double Tab; // mean fertility loss time
  double Tr,Td; // mean recovery and degradation time for max and min exposure to nature
  double d; // decay distance for ecosystem service delivery

  double dtsave; // timestep for saving data

  int seed; // this is expid
  unsigned long int seed2;

  ///////////////////////////////////////////////////////////////////////////////
  // IMPORT PARAMETER VALUES
  ///////////////////////////////////////////////////////////////////////////////

  if (argc>1) {
        char * pEnd;

        // time and space specifications for the simulation
        SimTime = strtod(argv[1], &pEnd);
        dtp = strtod(argv[2], &pEnd);
        n = atoi(argv[3]);

        // initial values of agricultural land use and consumption
        a0 = strtod(argv[4],&pEnd);
        d0 = strtod(argv[5],&pEnd);

        // agricultural production parameters
        ksi = strtod(argv[6], &pEnd);
        y0 = strtod(argv[7], &pEnd);
        sar = strtod(argv[8], &pEnd);

        // human action parameters
        a = strtod(argv[9], &pEnd);
        w = strtod(argv[10], &pEnd);
        Tag = strtod(argv[11], &pEnd);

        // abandonment parameters
        Tab = strtod(argv[12], &pEnd);

        // spontaneous evolution parameters
        Tr = strtod(argv[13], &pEnd);
        Td = strtod(argv[14], &pEnd);

        // distance for es provision
        d = strtod(argv[15], &pEnd);

        // save timespace just in case
        dtsave = strtod(argv[16], &pEnd);

        // save seed
        seed = atoi(argv[17]);
  }

  /////////////////////////////////////////////////////////////////////////////
  // CREATION OF DATA FILES
  /////////////////////////////////////////////////////////////////////////////

  //creating data directory with today's date
  auto tt = time(nullptr);
  auto tm = *localtime(&tt);
  ostringstream oss;
  oss << put_time(&tm, "%d-%m-%Y");
  string str_date = oss.str();
  string dirname = "DATA_"+str_date;
  // not using this: adding a count if several simulation runs on the same day
  // unsigned int count=0;
  // string temp_dirname = dirname+"_"+to_string(count);
  // while (fs::exists(temp_dirname)){
  //   count+=1;
  //   temp_dirname=dirname+to_string(count);
  // }
  // dirname=temp_dirname;
  if (fs::exists(dirname)){
    // don't need to create it
  }
  else{
    fs::create_directory(dirname); // commenting to see if avoids problem with sensitivity OM
  }

  //creating vector of strings to store all the input arguments
  vector<string> allArgs(argv,argv+argc);
  string filename;
  if(argc>1){
    filename = "_T_"+allArgs[1];
    filename += "_dtp_"+allArgs[2];
    filename += "_n_"+allArgs[3];
    filename += "_a0_"+allArgs[4];
    filename += "_d0_"+allArgs[5];
    filename += "_ksi_"+allArgs[6];
    filename += "_y0_"+allArgs[7];
    filename += "_sar_"+allArgs[8];
    filename += "_a_"+allArgs[9];
    filename += "_w_"+allArgs[10];
    filename += "_Tag_"+allArgs[11];
    filename += "_Tab_"+allArgs[12];
    filename += "_Tr_"+allArgs[13];
    filename += "_Td_"+allArgs[14];
    filename += "_d_"+allArgs[15];
    filename += "_dtsave_"+allArgs[16];
    filename += "_expid_"+allArgs[17];
    filename+=".dat";
  }

  // string filename_AGRE=dirname+"/"+"DATA_AGRE"+filename;
  // string filename_AGRE="DATA_AGRE"+filename;
  string filename_AGRE="DATA_AGRE";
  // string filename_AGRE="DATA_AGRE";
  ofstream tofile_agre(filename_AGRE);
  tofile_agre.precision(5);
  tofile_agre.setf(ios::scientific,ios::floatfield);

  string filename_LAND=dirname+"/"+"DATA_LAND"+filename;
  ofstream tofile_land(filename_LAND);
  tofile_land.precision(5);
  tofile_land.setf(ios::scientific,ios::floatfield);

  string filename_CLUS=dirname+"/"+"DATA_CLUS"+filename;
  ofstream tofile_clus(filename_CLUS);
  tofile_clus.precision(5);
  tofile_clus.setf(ios::scientific,ios::floatfield);

  string filename_CONF=dirname+"/"+"DATA_CONF"+filename;
  ofstream tofile_conf(filename_CONF);
  tofile_conf.precision(5);
  tofile_conf.setf(ios::scientific,ios::floatfield);

  // string filename_SENS="sensitivityOut.dat";
  // ofstream tofile_sens(filename_SENS);
  // tofile_sens.precision(5);
  // tofile_sens.setf(ios::scientific,ios::floatfield);
  //
  // string filename_SPEX="SPEXOut.dat";
  // ofstream tofile_spex(filename_SPEX);
  // tofile_spex.precision(5);
  // tofile_spex.setf(ios::scientific,ios::floatfield);

  /////////////////////////////////////////////////////////////////////////////
  // seeding the random number generator

  seed2=abs(seed);
  // seeding the random double generator: used for gillespie
  if (SEED_MT==1){
    seedMT(seed2);
  }
  seedMT2();

  // creating the random integer generator and seeding it
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, seed2);

  /////////////////////////////////////////////////////////////////////////////
  // VARIABLE DECLARATION AND INITIALISATION
  /////////////////////////////////////////////////////////////////////////////

  double t=0;
  double t_save=0;
  double dtg;
  double dt=dtp;
  double x_rand;
  unsigned int reaction,patch;
  unsigned int i;

  // this vector has only one member and it is the population
  vector<double> population;
  // vector containing the landscape state
  vector<unsigned int> landscape;
  // vector containing neighbours
  vector<vector<unsigned int>> neighbourMatrixES;
  vector<vector<unsigned int>> neighbourMatrix;
  // vector containing the production of each patch
  vector<double> agriculturalProduction;
  // vector containing the production of each patch
  vector<double> ecosystemServices;
  // vector containing the natural connected components information
  vector<vector<int>> naturalComponents;
  // vector containing the event's propensities
  vector<double> propensityVector;
  // vector to store the number of events of each kind
  vector<unsigned int> count_events={0,0,0,0,0,0};

  ////////////////////////////////////////////////////////////////////////////
  // STATE INITIALISATION
  ////////////////////////////////////////////////////////////////////////////

  // BY CONF FILE
  if (LOAD_CONF==1){
    cout << "Starting from conf file \n";
    ifstream conf_file("DATA_CONF_T_2000_dtp_0.1_n_40_p0_5_cg0_1_r0_1_y0_5_a_0.45_w_0_m_1_g_1_l_0.1_Ta_1_phi_3_Tr_1.5_Td_50_dtsave_1.dat");
    if(conf_file.is_open()) {

      // first extract the time and population
      double pop;
      if (!(conf_file >> t >> pop)){
        cout << "Error: gill_main.cpp: time and population could not be loaded from CONF file. \n";
      }
      population.push_back(pop);
      SimTime+=t;
      // extracting by token moves forward the pointer in the file
      // now extract the landscape
      unsigned int state;
      i=0;
      while(conf_file >> state){
        landscape.push_back(state);
        i+=1;
      }
    }
  }
  else{ // WITH ARGV PARAMETERS
    getNeighbourMatrix(neighbourMatrix,n,1);
    getNeighbourMatrix(neighbourMatrixES,n,d);
    initializeSES(landscape,population,naturalComponents,agriculturalProduction,ecosystemServices,neighbourMatrix,neighbourMatrixES,n,a0,d0,a,ksi,y0,sar,w,r);
  }

  /////////////////////////////////////////////////////////////////////////////
  // BEGIN OF SIMULATION
  /////////////////////////////////////////////////////////////////////////////

  unsigned int nat_cells;
  unsigned int deg_cells;

  // calculate the number of natural cells for the nopop experiment
  nat_cells = 0;
  deg_cells = 0;
  for(i=0;i<landscape.size();i++){
    if(landscape[i]==0){
      nat_cells+=1;
    }
    else if(landscape[i]==1){
      deg_cells+=1;
    }
  }

  unsigned int nMin=0;
  unsigned int nMax=0;
  double pMin=0;
  double pMax=0;
  unsigned int first_time=0;

  // entering the time loop
  while(t<SimTime){

    // updating agricultural production
    getAgriculturalProduction(agriculturalProduction, landscape, ecosystemServices, ksi, y0);

    // STOPPING EXECUTION AS SOON AS LANDSCAPE IS FULLY NATURAL OR DEGRADED
    nat_cells = 0;
    deg_cells = 0;
    for(i=0;i<landscape.size();i++){
      if(landscape[i]==0){
        nat_cells+=1;
      }
      else if(landscape[i]==1){
        deg_cells+=1;
      }
    }
    if(nat_cells==landscape.size() || deg_cells==landscape.size()){
      break;
    }

    ///////////////////////////////////////////////////////////////////////////
    // CALCULATING THE MINIMUM AND MAXIMUM VARAIBLE VALUES TO GET CYCLES
    // FOR THIS EXPERIMENT I SET UP THE INITIAL CONDITIONS AT THE TRANSITION POINT
    //////////////////////////////////////////////////////////////////////////

    if(t>SimTime/6){ // let some time for a transient before the cycles

      if (first_time==0){ // to initialize the value after the transient
        nMax = nat_cells;
        nMin = nat_cells;
        pMax = population[0];
        pMin = population[0];
        first_time=1;
      }

      // i only save the natural area and population for instance
      if (nat_cells>nMax){
        nMax=nat_cells; // reset the maximum value
      }
      if (nat_cells<nMin){
        nMin=nat_cells; // reset the minimum value
      }
      if (population[0]>pMax){
        pMax=population[0]; // reset the maximum value
      }
      if (population[0]<pMin){
        pMin=population[0]; // reset the minimum value
      }

    }


    ///////////////////////////////////////////////////////////////////////////
    // SAVING DATA
    ///////////////////////////////////////////////////////////////////////////
    if(t>=t_save)
    {
      saveAggregated(tofile_agre,t,population,landscape,agriculturalProduction,naturalComponents,ecosystemServices,n,2,(double)nMax/landscape.size(),(double)nMin/landscape.size(),pMax,pMin);
      saveLandscape(tofile_land,t,landscape);
      saveComponents(tofile_clus,t,landscape,naturalComponents);

      t_save+=dtsave;
    }


    ///////////////////////////////////////////////////////////////////////////
    // CALCULATING PROPENSITY VECTOR
    ///////////////////////////////////////////////////////////////////////////
    getPropensityVector(propensityVector,neighbourMatrix,landscape,ecosystemServices,agriculturalProduction,population,Tr,Td,w,a,Tag,Tab);
    //cout << "size of pvector is " << propensityVector.size() << "\n";
    ///////////////////////////////////////////////////////////////////////////
    // TIME UNTIL NEXT EVENT
    ///////////////////////////////////////////////////////////////////////////
    dtg=-1/propensityVector.back()*log(ranMT());
    ///////////////////////////////////////////////////////////////////////////
    // LOOKING IF NEXT THING TO DO IS TO UPDATE POPULATION AND CONSUMPTION OR
    // THE REALIZATION OF A STOCHASTIC EVENT
    ///////////////////////////////////////////////////////////////////////////

    if (dtg>dt){ // if the time until next event is larger than the ODE timestep
      // update population and consumption
      if (population[0]>0){
        rungeKutta4(population,agriculturalProduction,dt);
      }
      else{
        population[0]=0;
        // break;
      }

      // update the time as well as the timestep for ODE solving
      t+=dt;
      dt=dtp;
    }
    else{ // if the time until next event is shorter than the ODE timestep
      // compute random number to select next reaction
      x_rand = ranMT()*propensityVector.back();
      // traverse the propensity vector and stop once reaching the selceted cell
      i=0;
      while(x_rand>propensityVector[i]){
        i++;
      }
      // calculate the corresponding reaction and patch from the selected index i
      reaction=(int)i/(n*n); //result from euclidian division
      patch=i%(n*n); // remainder from euclidian division

      // transform the landscape according to reaction and patch
      if (reaction==0){landscape[patch]=0;count_events[0]+=1; updateNCCadding(naturalComponents,neighbourMatrix,landscape,patch); getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrixES,landscape,sar); } //recovery
      else if(reaction==1) {landscape[patch]=1;count_events[1]+=1; updateNCCremoving(naturalComponents,landscape,patch); getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrixES,landscape,sar); } //degradation
      else if(reaction==2) {landscape[patch]=2;count_events[2]+=1;updateNCCremoving(naturalComponents,landscape,patch);getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrixES,landscape,sar);} //expansion
      else if(reaction==3) {landscape[patch]=3;count_events[3]+=1;} //intensification
      else if(reaction==4) {landscape[patch]=0;count_events[4]+=1; updateNCCadding(naturalComponents,neighbourMatrix,landscape,patch); getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrixES,landscape,sar); } //abandonment to natural
      else if(reaction==5) {landscape[patch]=1;count_events[5]+=1;} //abandonment to degraded
      else {cout << "Error: gill_main.cpp reaction " << reaction << " does not exist.\n";}


      // update the time and timestep for ODE solving
      t+=dtg;
      dt-=dtg;
    }
  }

  // saving CONF file to re start other simulations from this point
  // tofile_conf << t << " " << population[0];
  // for(i=0 ; i<landscape.size() ; i++){
  //   tofile_conf << " " << landscape[i];
  // }
  // tofile_conf << "\n";

  // saving files so ifdtsave was largest than execution time one gets the final
  // values for every output we look at
  // ofstream tofile_sens("DATA_SENSITIVITY");
  // tofile_sens.precision(5);
  // tofile_sens.setf(ios::scientific,ios::floatfield);
  // saveAggregated(tofile_sens,t,population,landscape,agriculturalProduction,naturalComponents,ecosystemServices,n,2,(double)nMax/landscape.size(),(double)nMin/landscape.size(),pMax,pMin);

  saveAggregated(tofile_agre,t,population,landscape,agriculturalProduction,naturalComponents,ecosystemServices,n,2,(double)nMax/landscape.size(),(double)nMin/landscape.size(),pMax,pMin);

  // saveLandscape(tofile_land,t,landscape);
  // saveComponents(tofile_clus,t,landscape,naturalComponents);

  // saving output for sensitivity analysis
  // saveSensitivityOutput(tofile_sens,n,1,population,naturalComponents,landscape,ecosystemServices);
  // saving output for pattern exploration space and origin exploration space
  // saveAggregated(tofile_spex,t,population,landscape,agriculturalProduction);

  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::minutes>(stop - start);
  // cout << "simulation time " << t << "\n";
  cout << "total execution time " << duration.count() << endl;

  return 0;
}
