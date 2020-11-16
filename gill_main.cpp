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

///////////////////////////////////////////////////////////////////////////////
// MAIN PROGRAM
///////////////////////////////////////////////////////////////////////////////

// time at beginning
auto start = chrono::high_resolution_clock::now();
auto start2 = chrono::high_resolution_clock::now();

int main(int argc, const char * argv[]){

  // seeding the random double generator: used for gillespie
  seedMT2();
  // creating the random integer generator and seeding it
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, 1425638);

  ///////////////////////////////////////////////////////////////////////////////
  // PARAMETER DECLARATION
  ///////////////////////////////////////////////////////////////////////////////

  int n;  // lenght of the sides of the square landscape: number of cells=n*n
  int SimTime; // total simulation time
  double dtp; // timestep for population dynamics

  double ao0, ai0; //number of agricultural patches at beggining
  double c0; // initial per capita consumption

  double r0; // population growth rate
  double ys0,yn0; // productivity per unit input or unit ecosystem service
  double pSD; // perturbation sd for production
  double ori,ini; // organic and intense inputs
  double ess; // ecosystem service saturation exponent
  double a; // expansion probability
  double w; // agricultural clustering parameter
  double m,m0; // maintenance cost per unit of input
  double g; // action probability per unit time per unit of consumption deficit
  double Ta; // mean abandonment time per maintenance deficit
  double Tr,Td; // mean recovery and degradation time for max and min exposure to nature
  double kg,kd,cmin; // growth and decrease rates of per capita consumption
  double dtsave; // timestep for saving data

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
        ao0 = strtod(argv[4],&pEnd);
        ai0 = strtod(argv[5],&pEnd);
        c0 = strtod(argv[6], &pEnd);

        // population growth rate
        r0 = strtod(argv[7], &pEnd);

        // consumption growth and decrease rates
        kg = strtod(argv[8], &pEnd);
        kd = strtod(argv[9], &pEnd);
        cmin = strtod(argv[10], &pEnd);

        // agricultural production parameters
        ys0 = strtod(argv[11], &pEnd);
        yn0 = strtod(argv[12], &pEnd);
        pSD = strtod(argv[13], &pEnd);
        ori = strtod(argv[14], &pEnd);
        ini = strtod(argv[15], &pEnd);
        ess = strtod(argv[16], &pEnd);

        // human action parameters
        a = strtod(argv[17], &pEnd);
        w = strtod(argv[18], &pEnd);
        g = strtod(argv[19], &pEnd);

        // abandonment parameters
        m0= strtod(argv[20], &pEnd);
        m= strtod(argv[21], &pEnd);
        Ta = strtod(argv[22], &pEnd);

        // spontaneous evolution parameters
        Tr = strtod(argv[23], &pEnd);
        Td = strtod(argv[24], &pEnd);

        // save timespace just in case
        dtsave = strtod(argv[25], &pEnd);

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
    fs::create_directory(dirname);
  }

  //creating vector of strings to store all the input arguments
  vector<string> allArgs(argv,argv+argc);
  string filename;
  if(argc>1){
    filename = "_T_"+allArgs[1];
    filename += "_dtp_"+allArgs[2];
    filename += "_n_"+allArgs[3];
    filename += "_ao0_"+allArgs[4];
    filename += "_ai0_"+allArgs[5];
    filename += "_c0_"+allArgs[6];
    filename += "_r0_"+allArgs[7];
    filename += "_kg_"+allArgs[8];
    filename += "_kd_"+allArgs[9];
    filename += "_cmin_"+allArgs[10];
    filename += "_ys0_"+allArgs[11];
    filename += "_yn0_"+allArgs[12];
    filename += "_pSD_"+allArgs[13];
    filename += "_ori_"+allArgs[14];
    filename += "_ini_"+allArgs[15];
    filename += "_ess_"+allArgs[16];
    filename += "_a_"+allArgs[17];
    filename += "_w_"+allArgs[18];
    filename += "_g_"+allArgs[19];
    filename += "_m0_"+allArgs[20];
    filename += "_m_"+allArgs[21];
    filename += "_Ta_"+allArgs[22];
    filename += "_Tr_"+allArgs[23];
    filename += "_Td_"+allArgs[24];
    filename += "_dtsave_"+allArgs[25];
    filename += "_expid_"+allArgs[26];
    filename+=".dat";
  }

  string filename_POPU=dirname+"/"+"DATA_POPU"+filename;
  ofstream tofile_popu(filename_POPU);
  tofile_popu.precision(5);
  tofile_popu.setf(ios::scientific,ios::floatfield);

  string filename_LAND=dirname+"/"+"DATA_LAND"+filename;
  ofstream tofile_land(filename_LAND);
  tofile_land.precision(5);
  tofile_land.setf(ios::scientific,ios::floatfield);

  string filename_EVNT=dirname+"/"+"DATA_EVNT"+filename;
  ofstream tofile_evnt(filename_EVNT);
  tofile_evnt.precision(5);
  tofile_evnt.setf(ios::scientific,ios::floatfield);

  string filename_CLUS=dirname+"/"+"DATA_CLUS"+filename;
  ofstream tofile_clus(filename_CLUS);
  tofile_clus.precision(5);
  tofile_clus.setf(ios::scientific,ios::floatfield);

  string filename_CONF=dirname+"/"+"DATA_CONF"+filename;
  ofstream tofile_conf(filename_CONF);
  tofile_conf.precision(5);
  tofile_conf.setf(ios::scientific,ios::floatfield);

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
  // this vector only has one member and it is the per capita consumption
  vector<double> consumption;
  // vector containing the landscape state
  vector<unsigned int> landscape;
  // vector containing the production of each patch
  vector<double> agriculturalProduction;
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
    initializeSES(landscape,population,consumption,naturalComponents,agriculturalProduction,c0,n,ao0,ai0,r,ys0,yn0,pSD,ori,ini,ess);
  }

  /////////////////////////////////////////////////////////////////////////////
  // BEGIN OF SIMULATION
  /////////////////////////////////////////////////////////////////////////////

  // entering the time loop
  while(t<SimTime){

    ///////////////////////////////////////////////////////////////////////////
    // SAVING DATA
    ///////////////////////////////////////////////////////////////////////////
    if(t>=t_save)
    {

      tofile_popu << t << " " << population[0] << " " << consumption[0] << "\n";
      tofile_land << t << " ";
      tofile_evnt << t << " ";
      tofile_clus << -1 << " " << t << "\n";

      for(i=0 ; i<landscape.size() ; i++){
        tofile_land << landscape[i] << " ";
      }
      for(i=0 ; i<count_events.size() ; i++){
        tofile_evnt << count_events[i] << " ";
      }
      for(i=0 ; i<naturalComponents.size() ; i++){
        tofile_clus << naturalComponents[i][0] << " " << naturalComponents[i][1] << "\n";
      }

      tofile_land << "\n";
      tofile_evnt << "\n";

      t_save+=dtsave;
    }

    ///////////////////////////////////////////////////////////////////////////
    // CALCULATING PROPENSITY VECTOR
    ///////////////////////////////////////////////////////////////////////////
    getPropensityVector(propensityVector,landscape,naturalComponents,agriculturalProduction,population,consumption,n,Tr,Td,w,a,g,Ta,ori,ini,ess,m0,m);
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
        rungeKutta4(population,consumption,agriculturalProduction,dt,r0,kg,kd,cmin);
      }
      else{
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
      if (reaction==0){landscape[patch]=0;count_events[0]+=1;} //recovery
      else if(reaction==1) {landscape[patch]=1;count_events[1]+=1;} //degradation
      else if(reaction==2) {landscape[patch]=2;count_events[2]+=1;} //expansion
      else if(reaction==3) {landscape[patch]=3;count_events[3]+=1;} //intensification
      else if(reaction==4) {landscape[patch]=0;count_events[4]+=1;} //abandonment to natural
      else if(reaction==5) {landscape[patch]=1;count_events[5]+=1;} //abandonment to degraded
      else {cout << "Error: gill_main.cpp reaction " << reaction << " does not exist.\n";}
      // updating natural connected components
      getNaturalConnectedComponents(naturalComponents, n, landscape);
      // updating agricultural production
      getAgriculturalProduction(agriculturalProduction,landscape,naturalComponents,n,ys0,yn0,pSD,ori,ini,ess,r);

      // update the time and timestep for ODE solving
      t+=dtg;
      dt-=dtg;
    }
  }

  // saving CONF file to re start other simulations from this point
  tofile_conf << t << " " << population[0] << consumption[0];
  for(i=0 ; i<landscape.size() ; i++){
    tofile_conf << " " << landscape[i];
  }
  tofile_conf << "\n";

  /* //print landscape to check
  for(i=0  ; i<n ; i++)
  {
    for(j=0 ; j<n ; j++)
    {
      cout << landscape[j+i*n] << " ";
    }
    cout << "\n";
  }
  for(i=0 ; i<naturalComponents.size() ; i++)
  {
    cout << "patch " << naturalComponents[i][0] << " is in cluster " << naturalComponents[i][1] <<  "\n";
    cout << "state of patch " << naturalComponents[i][0] << " is " << landscape[naturalComponents[i][0]] << "\n";
  }
  for(i=0  ; i<landscape.size() ; i++)
  {
    cout << "patch " << i << " is in state " << landscape[i] << "\n";
  }*/
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::minutes>(stop - start);
  cout << "total execution time " << duration.count() << endl;
  cout << "Time out of simulation =" << t <<", last dtg= " << dtg << ", Cumsum= " << propensityVector.back() <<  "\n";
  cout << "looking for nan emplacement... \n";
  for(ix=0;ix<propensityVector.size();ix++){
    if isnan(propensityVector[ix]){
      i=ix;
      break;
    }
  }
  reaction=(int)i/(n*n); //result from euclidian division
  patch=i%(n*n); // remainder from euclidian division
  cout << "Error in reaction "<< reaction << " at patch " << patch <<"\n";

  return 0;
}
