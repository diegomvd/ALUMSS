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
  int SimTime; // total simulation time
  double dtp; // timestep for population dynamics

  double a0; //number of agricultural patches at beggining
  double d0; // number of degraded patches at beggining

  double ksi; // productivity per unit input or unit ecosystem service
  double sar; // ecosystem service saturation exponent
  double a; // intensification probability
  double w; // agricultural clustering parameter
  double Tag; // action probability per unit time per unit of consumption deficit
  double Tab; // mean fertility loss time
  double Tr,Td; // mean recovery and degradation time for max and min exposure to nature
  double e12,c12; // half saturation values for es provision and consumption deficit
  double d; // distance at which eecosystem services are delivered

  double dtsave; // timestep for saving data

  unsigned long int seed; // this is expid

  ///////////////////////////////////////////////////////////////////////////////
  // IMPORT PARAMETER VALUES
  ///////////////////////////////////////////////////////////////////////////////

  if (argc>1) {
        char * pEnd;

        // time and space specifications for the simulation
        SimTime = atoi(argv[1]);
        dtp = strtod(argv[2], &pEnd);
        n = atoi(argv[3]);

        // initial values of agricultural land use and consumption
        a0 = strtod(argv[4],&pEnd);
        d0 = strtod(argv[5],&pEnd);

        // agricultural production parameters
        ksi = strtod(argv[6], &pEnd);
        sar = strtod(argv[7], &pEnd);

        // human action parameters
        a = strtod(argv[8], &pEnd);
        w = strtod(argv[9], &pEnd);
        Tag = strtod(argv[10], &pEnd);

        // abandonment parameters
        Tab = strtod(argv[11], &pEnd);

        // spontaneous evolution parameters
        Tr = strtod(argv[12], &pEnd);
        Td = strtod(argv[13], &pEnd);

        // half values for saturation functions
        e12 = strtod(argv[14], &pEnd);
        c12 = strtod(argv[15], &pEnd);

        // distance for es provision
        d = strtod(argv[16], &pEnd);

        // save timespace just in case
        dtsave = strtod(argv[17], &pEnd);

        // save seed
        seed = abs(atoi(argv[18]));
  }

  // seeding the random double generator: used for gillespie
  if (SEED_MT==1){
    seedMT(seed);
  }
  seedMT2();

  // creating the random integer generator and seeding it
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, seed);

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
    fs::create_directory(dirname);
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
    filename += "_sar_"+allArgs[7];
    filename += "_a_"+allArgs[8];
    filename += "_w_"+allArgs[9];
    filename += "_Tag_"+allArgs[10];
    filename += "_Tab_"+allArgs[11];
    filename += "_Tr_"+allArgs[12];
    filename += "_Td_"+allArgs[13];
    filename += "_d_"+allArgs[14];
    filename += "_e12_"+allArgs[15];
    filename += "_c12_"+allArgs[16];
    filename += "_dtsave_"+allArgs[17];
    filename += "_expid_"+allArgs[18];
    filename+=".dat";
  }

  string filename_AGRE=dirname+"/"+"DATA_AGRE"+filename;
  ofstream tofile_agre(filename_AGRE);
  tofile_agre.precision(5);
  tofile_agre.setf(ios::scientific,ios::floatfield);

  string filename_LAND=dirname+"/"+"DATA_LAND"+filename;
  ofstream tofile_land(filename_LAND);
  tofile_land.precision(5);
  tofile_land.setf(ios::scientific,ios::floatfield);

  string filename_METR=dirname+"/"+"DATA_METR"+filename;
  ofstream tofile_metr(filename_METR);
  tofile_metr.precision(5);
  tofile_metr.setf(ios::scientific,ios::floatfield);

  string filename_CLUS=dirname+"/"+"DATA_CLUS"+filename;
  ofstream tofile_clus(filename_CLUS);
  tofile_clus.precision(5);
  tofile_clus.setf(ios::scientific,ios::floatfield);

  string filename_RIPL=dirname+"/"+"DATA_RIPL"+filename;
  ofstream tofile_ripl(filename_RIPL);
  tofile_ripl.precision(5);
  tofile_ripl.setf(ios::scientific,ios::floatfield);

  string filename_CONF=dirname+"/"+"DATA_CONF"+filename;
  ofstream tofile_conf(filename_CONF);
  tofile_conf.precision(5);
  tofile_conf.setf(ios::scientific,ios::floatfield);

  string filename_SENS="sensitivityOut.dat";
  ofstream tofile_sens(filename_SENS);
  tofile_sens.precision(5);
  tofile_sens.setf(ios::scientific,ios::floatfield);

  string filename_SPEX="SPEXOut.dat";
  ofstream tofile_spex(filename_SPEX);
  tofile_spex.precision(5);
  tofile_spex.setf(ios::scientific,ios::floatfield);


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
    initializeSES(landscape,population,naturalComponents,agriculturalProduction,ecosystemServices,neighbourMatrix,neighbourMatrixES,n,a0,d0,a,ksi,e12,sar,d,w,r);
  }

  /////////////////////////////////////////////////////////////////////////////
  // BEGIN OF SIMULATION
  /////////////////////////////////////////////////////////////////////////////

  // for the nopop experiment
  unsigned int nat_cells;

  // entering the time loop
  while(t<SimTime){

    // updating agricultural production
    getAgriculturalProduction(agriculturalProduction, landscape, ecosystemServices, ksi, e12);

    ///////////////////////////////////////////////////////////////////////////
    // SAVING DATA
    ///////////////////////////////////////////////////////////////////////////
    if(t>=t_save)
    {
      saveAggregated(tofile_agre,t,population,landscape,agriculturalProduction);
      saveLandscape(tofile_land,t,landscape);
      saveComponents(tofile_clus,t,landscape,naturalComponents);
      saveLandMetrics(tofile_metr,t,naturalComponents,ecosystemServices);
      saveRipley(tofile_ripl,t,n,landscape,1);

      t_save+=dtsave;
    }


    ///////////////////////////////////////////////////////////////////////////
    // CALCULATING PROPENSITY VECTOR
    ///////////////////////////////////////////////////////////////////////////
    getPropensityVector(propensityVector,neighbourMatrix,landscape,ecosystemServices,agriculturalProduction,population,Tr,Td,w,a,Tag,Tab,e12,c12);
    //cout << "size of pvector is " << propensityVector.size() << "\n";
    ///////////////////////////////////////////////////////////////////////////
    // TIME UNTIL NEXT EVENT
    ///////////////////////////////////////////////////////////////////////////
    dtg=-1/propensityVector.back()*log(ranMT());
    ///////////////////////////////////////////////////////////////////////////
    // LOOKING IF NEXT THING TO DO IS TO UPDATE POPULATION AND CONSUMPTION OR
    // THE REALIZATION OF A STOCHASTIC EVENT
    ///////////////////////////////////////////////////////////////////////////

    // calculate the number of natural cells for the nopop experiment
    for(i=0;i<landscape.size();i++){
      if(landscape[i]==0){
        nat_cells+=1;
      }
    }
    if(nat_cells==landscape.size()){
      break;
    }

    if (dtg>dt){ // if the time until next event is larger than the ODE timestep
      // update population and consumption
      if (population[0]>0){
        rungeKutta4(population,agriculturalProduction,dt);
      }
      else{
        population[0]=0;
        break;
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
  tofile_conf << t << " " << population[0];
  for(i=0 ; i<landscape.size() ; i++){
    tofile_conf << " " << landscape[i];
  }
  tofile_conf << "\n";

  // saving output for sensitivity analysis
  saveSensitivityOutput(tofile_sens,n,1,population,naturalComponents,landscape,ecosystemServices);
  // saving output for pattern exploration space and origin exploration space
  // saveAggregated(tofile_spex,t,population,landscape,agriculturalProduction);

  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::minutes>(stop - start);
  cout << "total execution time " << duration.count() << endl;

  return 0;
}
