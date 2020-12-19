#include "gill_functions.h"

#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
void getNeighbourMatrix(vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, unsigned int d)
{

  /*fills a vector containing the neighbours indexes for each patch*/

  int ix,jx,xi,yi,xj,yj;
  unsigned int dx,dy;
  unsigned int manhattanDist;

  neighbourMatrix.resize(n*n);

  for (ix=0; ix<neighbourMatrix.size(); ++ix){
    xi = ix%n;
    yi = (int)ix/n;
    for (xj=0; xj<n; xj++){
      dx=abs(xi-xj);
      // calculating cyclic distances to account for periodic borders
      if (dx>n/2){
        dx=n-dx;
      }
      for (yj=0; yj<n; yj++){
        dy=abs(yi-yj);
        // calculating cyclic distances to account for periodic borders
        if (dy>n/2){
          dy=n-dy;
        }
        manhattanDist = dx+dy;
        if (manhattanDist<=d and manhattanDist>0){
          jx = xj + yj*n;
          neighbourMatrix[ix].push_back(jx);
        }
      }
    }
  }


  return;
}

void getNeighbours(vector<unsigned int> &neighboursList, const vector<vector<unsigned int>> &neighbourMatrix, unsigned int i)
{
  neighboursList = neighbourMatrix[i];
  return;
}

void getNeighboursState(vector<unsigned int> &neighboursState, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, unsigned int i, unsigned int state)
{
  /*
  fills the vector neighboursState so that it contains the indexes of all the
  closest neighbours of i in a given state. the landscape is passed as a constant
  reference so that the vector cannot be modified by the function in main
  */

  /*
  getting the neighbours indexes in neighbour_list vector
  */
  vector<unsigned int> neighboursList;
  getNeighbours(neighboursList,neighbourMatrix,i);

  /*
  getting the index of neighbours in the wanted state
  */
  unsigned long ix;
  for (ix=0 ; ix<neighboursList.size() ; ++ix){
    if (landscape[neighboursList[ix]] == state) {
      neighboursState.push_back( neighboursList[ix] );
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// 2- Calculation of Ecosystem Service provision:
//       - getNaturalConnectedComponents
//       - getEcosystemServiceProvision
//       - getAgriculturalProduction
////////////////////////////////////////////////////////////////////////////////

void getNaturalConnectedComponents(vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape)
{
  /*
  fills a vector where each member is a vector containing the indexes of all the
  natural patches belonging to the same cluster
  */

  vector<unsigned int> naturalPatches;
  unsigned int manhattanDist;
  unsigned int i, j;
  int xi, xj, yi, yj;
  unsigned int dx, dy;


  unsigned int n = (unsigned int) sqrt(landscape.size());
  // clearing natural components vector for refilling
  naturalComponents.clear();

  /*
  get the list of natural patches in the landscape
  */
  for(i=0 ; i<landscape.size() ; ++i){
    if (landscape[i]==0){
      naturalPatches.push_back(i);
    }
  }

  /*
  create an undirected graph with the set of natural patches to calculate
  the connected components. to estimate whether two patches are connected
  we calculate the manhattan distance between them
  */

  using namespace boost;
  {
    typedef adjacency_list< vecS, vecS, undirectedS > Graph;
    Graph G;


    for(i=0 ; i<naturalPatches.size() ; ++i){
      add_edge(i,i,G);
      // converting 1-D coordinates to 2-D
      xi=naturalPatches[i]%n;
      yi=(int)naturalPatches[i]/n;
      for(j=0 ; j<i ; ++j){
        // converting 1-D coordinates to 2-D
        xj=naturalPatches[j]%n;
        yj=(int)naturalPatches[j]/n;
        // calculating manhattan distance between points
        dx=abs(xi-xj);
        dy=abs(yi-yj);
        // calculating cyclic distances to account for periodic borders
        if (dx>n/2){
          dx=n-dx;
        }
        if (dy>n/2){
          dy=n-dy;
        }
        //
        manhattanDist=dx+dy;
        if ( manhattanDist<2 ){
          add_edge(i, j, G);
        }
      }
    }

    /*
    initializing the vector containing the components and calculating components
    */
    vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);

    /*
    converting the nodes indexes into actual landscape coordinates
    */
    naturalComponents.resize(num);
    for (i=0 ; i<naturalComponents.size() ; ++i){
      for (j=0; j<component.size(); j++){
        if(component[j]==i){
          naturalComponents[i].push_back(naturalPatches[j]);
        }
      }
    }

  }
  return;
}

void updateNCCadding(vector<vector<int>> &naturalComponents, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, unsigned int i)
{
  vector<unsigned int> neighboursNatural;
  getNeighboursState(neighboursNatural,neighbourMatrix,landscape,i,0); // state 0 is natural

  if(neighboursNatural.size()==0){ //no natural neighbour
      vector<int> newNaturalComponent;
      newNaturalComponent.push_back(i); // create a new natural component
      naturalComponents.push_back(newNaturalComponent); // add it to the list
  }
  else{
    unsigned long ix,jx;
    vector<unsigned int> clusterList;
    for (ix=0;ix<neighboursNatural.size();ix++){
      // for each of the natural neighbours check their cluster membership and
      // store it to determine whether cluster have been merged or not
      for (jx=0; jx<naturalComponents.size(); jx++){
        // check if neighbour belongs to clster jx
        if (find( naturalComponents[jx].begin(),naturalComponents[jx].end(),neighboursNatural[ix]) != naturalComponents[jx].end()){
          //check if cluster jx is not in the clusterlist
          if (find( clusterList.begin(), clusterList.end(), jx) == clusterList.end()){
            // add cluster jx to list
            clusterList.push_back(jx);
            break;
          }
        }
      }
    }

    // now update getNaturalConnectedComponents
    if (clusterList.size()>0){
      // adding the new natural cell to the first cluster in the list
      naturalComponents[clusterList[0]].push_back(i);
      for (ix=1; ix<clusterList.size(); ix++){
        // putting all the members of the merging cluster in the first one of the list
        for(jx=0;naturalComponents[clusterList[ix]].size();jx++){
          naturalComponents[clusterList[0]].push_back(naturalComponents[clusterList[ix]][jx]);
        }
        // deleting this cluster from naturalComponents after the merge
        naturalComponents.erase(naturalComponents.begin()+clusterList[ix]);
      }
    }
    else{
      cout << "Error: clusterList should not be empty in updateNaturalConnectedComponents function.";
    }

  }

  return;
}

void updateNCCremoving(vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape, unsigned int l)
{

  unsigned int n = (unsigned int) sqrt(landscape.size());

  // find cluster of cell l
  unsigned long ix;
  vector<int>::iterator it;
  for (ix=0;ix<naturalComponents.size();ix++){
    it = find( naturalComponents[ix].begin(),naturalComponents[ix].end(),l);
    if (it!=naturalComponents[ix].end()){
      break;
    }
  }

  // erase cell i from naturalComponents
  naturalComponents[ix].erase(it);

  // fill natural patches to get connected natural components
  vector<unsigned int> naturalPatches;
  unsigned long jx;
  for (jx=0;jx<naturalComponents[ix].size();jx++){
    naturalPatches.push_back(naturalComponents[ix][jx]);
  }

  // erase concerned cluster from naturalComponents
  naturalComponents.erase(naturalComponents.begin()+ix);

  //get connected components from the naturalPatches
  unsigned int manhattanDist;
  unsigned int i, j;
  int xi, xj, yi, yj;
  unsigned int dx, dy;

  /*
  create an undirected graph with the set of natural patches to calculate
  the connected components. to estimate whether two patches are connected
  we calculate the manhattan distance between them
  */

  using namespace boost;
  {
    typedef adjacency_list< vecS, vecS, undirectedS > Graph;
    Graph G;


    for(i=0 ; i<naturalPatches.size() ; ++i){
      add_edge(i,i,G);
      // converting 1-D coordinates to 2-D
      xi=naturalPatches[i]%n;
      yi=(int)naturalPatches[i]/n;
      for(j=0 ; j<i ; ++j){
        // converting 1-D coordinates to 2-D
        xj=naturalPatches[j]%n;
        yj=(int)naturalPatches[j]/n;
        // calculating manhattan distance between points
        dx=abs(xi-xj);
        dy=abs(yi-yj);
        // calculating cyclic distances to account for periodic borders
        if (dx>n/2){
          dx=n-dx;
        }
        if (dy>n/2){
          dy=n-dy;
        }
        //
        manhattanDist=dx+dy;
        if ( manhattanDist<2 ){
          add_edge(i, j, G);
        }
      }
    }

    /*
    initializing the vector containing the components and calculating components
    */
    vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);

    // adding the new components to naturalComponents
    vector<int> newComponent;
    for(i=0; i<num; i++){
      for(j=0; j<component.size(); j++){
        if(i==component[j]){
          newComponent.push_back(naturalPatches[j]);
        }
      }
      naturalComponents.push_back(newComponent);
      newComponent.clear();
    }
  }

  return;
}

void getEcosystemServiceProvision(vector<double> &ecosystemServices, const vector<vector<int>> &naturalComponents, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, double sar)
{
  /*
  returns the exposure to the wanted state of patch i. currently it is only used
  for the exposure to nature. the exposure to nature scales with biodiversity
  hence like a SAR, where the area is the total natural area in contact with
  patch i
  */

  double area;
  double ecosystemServiceProvision;
  vector<unsigned int> neighboursState;
  double neighbourNumber = neighbourMatrix[0].size();
  unsigned long i,ix,jx;

  ecosystemServices.clear();

  for(i=0;i<landscape.size();i++){

    ecosystemServiceProvision=0;
    /*
    getting the state neighbours indexes in neighboursState vector
    */

    getNeighboursState(neighboursState,neighbourMatrix,landscape,i,0); // state 0 is natural

    /*
    calculate the area of each of the neighbour's component
    */
    for (ix=0;ix<neighboursState.size();ix++){
      // for each of the natural neighbours check their cluster membership
      for (jx=0; jx<naturalComponents.size(); jx++){
        // check if neighbour belongs to cluster jx
        if (find( naturalComponents[jx].begin(),naturalComponents[jx].end(),neighboursState[ix]) != naturalComponents[jx].end()){
          area=naturalComponents[jx].size()/landscape.size();
          ecosystemServiceProvision+=(1/neighbourNumber)*pow(area,sar);
          break;
        }
      }
    }

    ecosystemServices.push_back(ecosystemServiceProvision);
    neighboursState.clear();
  }

  return;
}

void getAgriculturalProduction(vector<double> &agriculturalProduction, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double ksi)
{
  /*
  returns the total agricultural production for a given "landscape" and
  minimum yield "y". the minimum yield is for a cropped patch with only
  non natural neighbours. natural neighbours raise yield.
  */

  agriculturalProduction.clear();

  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if(landscape[ix]==2){ // cropped patches
      agriculturalProduction.push_back( ksi*ecosystemServices[ix] ) ;
    }
    else if(landscape[ix]==3){ //intense
      agriculturalProduction.push_back( 1 );
    }
    else{
      agriculturalProduction.push_back(0);
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// 3- Calculation of events' propensities:
//       - getConsumptionDeficit
//       - getSpontaneousPropensity
//       - getAgroPropensity
//       - getAbandonmentPropensity
//       - getPropensityVector
////////////////////////////////////////////////////////////////////////////////

double getConsumptionDeficit(const vector<double> &agriculturalProduction, const vector<double> &population)
{
  double totalAgriculturalProduction=0;
  unsigned long ix;
  for(ix=0;ix<agriculturalProduction.size();ix++){
    totalAgriculturalProduction+=agriculturalProduction[ix];
  }
  return population[0] - totalAgriculturalProduction;
}

void getSpontaneousPropensity(vector<double> &recoveryPropensity, vector<double> &degradationPropensity, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double Tr, double Td)
{
  /*
  fills the recoveryPropensity vector and the degradationPropensity vector
  */

  /*
  traverse the landscape and calculate the recovery propensity for the degraded
  patches and the degradation propensity for the natural patches
  */
  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if (landscape[ix] == 1) { // if patch ix is degraded
      recoveryPropensity.push_back( ecosystemServices[ix]/Tr );
      degradationPropensity.push_back(0);
    }
    else if(landscape[ix] == 0){
      recoveryPropensity.push_back( 0 );
      degradationPropensity.push_back(1/Td*(1-ecosystemServices[ix]) );
    }
    else{
      recoveryPropensity.push_back( 0 );
      degradationPropensity.push_back( 0 );
    }
  }

  return;
}

void getAgroPropensity(vector<double> &expansionPropensity, vector<double> &intensePropensity, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, const vector<double> &agriculturalProduction, const vector<double> &population, double w, double a, double Tag)
{
  /*
  fills the cropping_propensity and restoring_propensity vector
  */

  // use them to store the accumulated probability in order to normalize afterwards
  double expansionCumSum = 0;
  double intenseCumSum = 0;

  // fraction of natural patches
  double fracN = 0;

  /*
  first checking if there is a consumption deficit that justifys human action
  */
  double consumptionDeficit = getConsumptionDeficit(agriculturalProduction,population);
  if (consumptionDeficit<0){ // if humans are satisfied...
    unsigned long ix;
    for(ix=0 ; ix<landscape.size() ; ++ix){
      expansionPropensity.push_back( 0 ); // ... they do not transform the landscape
      intensePropensity.push_back( 0 ); // ... they do not transform the landscape
    }
  }
  else{ // if humans are not happy with they current resource access ...

    /*
    in this loop we calculate the probabilities of choosing one patch over other
    for cropping as a function of the clustering parameter w. we also signal
    all the degraded patches that can be restored. the restoration probability
    is uniform over space.
    */
    unsigned int ix;
    vector<unsigned int> organicNeighbours;
    vector<unsigned int> intenseNeighbours;
    for (ix=0 ; ix<landscape.size() ; ++ix){

      if (landscape[ix] == 0){ // if patch is natural, hence can be converted to organic

        fracN+=1;

        // get the indexes of the organic neighbours
        organicNeighbours.clear();
        getNeighboursState(organicNeighbours,neighbourMatrix,landscape,ix, 2); // state 2 is organic
        // cropping probability expression is taken from bart's and dani's paper
        expansionPropensity.push_back( pow( max(0.1 , (double)organicNeighbours.size() ) , w ) );
        expansionCumSum += expansionPropensity.back();
        intensePropensity.push_back( 0 );

      }
      else if (landscape[ix] == 2 ){ // if patch is organic it can be intensifyed
        // get the indexes of the organic neighbours
        intenseNeighbours.clear();
        getNeighboursState(intenseNeighbours,neighbourMatrix,landscape,ix, 3); // state 3 is intense

        intensePropensity.push_back( pow( max(0.1 , (double)intenseNeighbours.size() ) , w ) );
        intenseCumSum += intensePropensity.back();
        expansionPropensity.push_back(0);

      }
      else{
        expansionPropensity.push_back(0);
        intensePropensity.push_back( 0 );
      }
    }

    /*
    now we normalize the previously obtained probabilities and get calculate the
    probability per unit time of action given the consumption deficit
    */
    fracN/=landscape.size();

    if (expansionCumSum>0 && intenseCumSum>0){ // this is to avoid dividing by zero
      for (ix=0; ix<landscape.size() ; ++ix){
        expansionPropensity[ix]=expansionPropensity[ix]/expansionCumSum/Tag*consumptionDeficit*pow(fracN,a);
        intensePropensity[ix]=intensePropensity[ix]/intenseCumSum/Tag*consumptionDeficit*(1-pow(fracN,a));
      }
    }
    else{
      if (expansionCumSum>0){
        for (ix=0; ix<landscape.size() ; ++ix){
          expansionPropensity[ix]=expansionPropensity[ix]/expansionCumSum*(1/Tag)*consumptionDeficit;
        }
      }
      else if (intenseCumSum>0){
        for (ix=0; ix<landscape.size() ; ++ix){
          intensePropensity[ix]=intensePropensity[ix]/intenseCumSum*(1/Tag)*consumptionDeficit;
        }
      }
    }
  }

  return;
}

void getAbandonmentPropensity(vector<double> &naturalAbandonPropensity, vector<double> &degradedAbandonPropensity, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, double Tab)
{
  /*
  fills the abandonment_propensity vector.
  */

  unsigned int ix;

  for (ix=0; ix<landscape.size(); ++ix){
    if (landscape[ix]==2){ // organic patch
      degradedAbandonPropensity.push_back(0);
      naturalAbandonPropensity.push_back(1/Tab*(1-ecosystemServices[ix]));
    }
    else if(landscape[ix]==3){ //intensive patch
      degradedAbandonPropensity.push_back( 1/Tab*(1-ecosystemServices[ix]) );
      naturalAbandonPropensity.push_back( 0 );
    }
    else{ // non cropped patches
      naturalAbandonPropensity.push_back(0);
      degradedAbandonPropensity.push_back(0);
    }
  }

  return;
}

void getPropensityVector(vector<double> &propensityVector, const vector<vector<unsigned int>> &neighbourMatrix, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices, const vector<double> &agriculturalProduction, const vector<double> &population, double Tr, double Td, double w, double a, double Tag, double Tab)
{
  /*
  calls all the functions to calculate the propensity of each event and merges
  them in a single propensity vector that can be used to run the gillespie algo
  */

  vector<double> recoveryPropensity;
  vector<double> degradationPropensity;
  vector<double> expansionPropensity;
  vector<double> intensePropensity;
  vector<double> naturalAbandonPropensity;
  vector<double> degradedAbandonPropensity;

  getSpontaneousPropensity(recoveryPropensity,degradationPropensity,landscape,ecosystemServices,Tr,Td);
  getAgroPropensity(expansionPropensity,intensePropensity,neighbourMatrix,landscape,agriculturalProduction,population,w,a,Tag);
  getAbandonmentPropensity(naturalAbandonPropensity,degradedAbandonPropensity,landscape,ecosystemServices,Tab);

  // clearing the previous propensity vector to refill it
  propensityVector.clear();

  // making sure the vector is not empty to avoid bad behaviour in next loop
  propensityVector.push_back(recoveryPropensity[0]);
  unsigned long ix;
  // over these loops the cumulative sum of each propensity vector is added to
  // the total propensity vector
  for (ix=1 ; ix<recoveryPropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+recoveryPropensity[ix]);
  }
  for (ix=0 ; ix<degradationPropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+degradationPropensity[ix]);
  }
  for (ix=0 ; ix<expansionPropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+expansionPropensity[ix]);
  }
  for (ix=0 ; ix<intensePropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+intensePropensity[ix]);
  }
  for (ix=0 ; ix<naturalAbandonPropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+naturalAbandonPropensity[ix]);
  }
  for (ix=0 ; ix<degradedAbandonPropensity.size() ; ++ix){
    propensityVector.push_back(propensityVector.back()+degradedAbandonPropensity[ix]);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// 4- Initialization functions:
//       - initializeLandscape
//       - initializePopulation
//       - initializeSES
////////////////////////////////////////////////////////////////////////////////

void initializeLandscape( vector<unsigned int> &landscape, vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, double ao0, double ai0, double w, gsl_rng  *r)
{
  /*
  initializes the landscape given a fraction of initial agricultural patches a0
  */

  //unsigned int number_cropped_patches = 1;
  unsigned long ix,jx,lx;
  // this vector contains the indexes of all the natural patches
  vector<double> probConversion;
  vector<double> probIntense;

  unsigned int nao0=(unsigned int) (ao0*n*n);
  unsigned int nai0=(unsigned int) (ai0*n*n);
  unsigned int na0=nao0+nai0;

  // first build a completely natural landscape with n*n patches
  landscape.push_back(0);
  probConversion.push_back(1);
  for (ix=1 ; ix<n*n; ++ix){
    landscape.push_back(0);
    probConversion.push_back(probConversion.back()+1);
    probIntense.push_back(0);
  }

  vector<unsigned int> organicNeighbours;
  for(ix=0;ix<na0;++ix){
    jx=0;
    while (gsl_rng_uniform(r)*probConversion.back()>probConversion[jx]){
      jx++;
    }
    landscape[jx]=2;
    probIntense[jx]=1;

    // recalculating probconversion
    organicNeighbours.clear();
    getNeighboursState(organicNeighbours,neighbourMatrix,landscape,0, 2);
    probConversion[0]=pow( max(0.1 , (double)organicNeighbours.size() ) , w ) ;
    for(lx=1;lx<probConversion.size();lx++){
      organicNeighbours.clear();
      getNeighboursState(organicNeighbours,neighbourMatrix,landscape,lx, 2);
      probConversion[lx]=probConversion[lx-1]+pow( max(0.1 , (double)organicNeighbours.size() ) , w ) ;
    }
    probConversion[jx]=0;
  }

  vector<unsigned int> intenseNeighbours;
  for(ix=0;ix<nai0;++ix){
    jx=0;
    while (gsl_rng_uniform(r)*probIntense.back()>probIntense[jx]){
      jx++;
    }
    landscape[jx]=3;

    // recalculating probconversion
    intenseNeighbours.clear();
    getNeighboursState(intenseNeighbours,neighbourMatrix,landscape,0, 3);
    probIntense[0]=pow( max(0.1 , (double)intenseNeighbours.size() ) , w ) ;
    for(lx=1;lx<probIntense.size();lx++){
      intenseNeighbours.clear();
      getNeighboursState(intenseNeighbours,neighbourMatrix,landscape,lx, 3);
      probIntense[lx]=probIntense[lx-1]+pow( max(0.1 , (double)intenseNeighbours.size() ) , w ) ;
    }
    probIntense[jx]=0;
  }

  return;
}

void initializePopulation( vector<double> &population, const vector<double> &agriculturalProduction)
{
  /*given an agricultural production, it sets the population at an equilibrium
  level
  */
  unsigned int ix;
  double totalAgriculturalProduction=0;
  for(ix=0; ix<agriculturalProduction.size(); ++ix){
    totalAgriculturalProduction+=agriculturalProduction[ix];
  }
  population.push_back(totalAgriculturalProduction);

  return;
}

void initializeSES( vector<unsigned int> &landscape, vector<double> &population, vector<vector<int>> &naturalComponents, vector<double> &agriculturalProduction, vector<double> &ecosystemServices, vector<vector<unsigned int>> &neighbourMatrix, vector<vector<unsigned int>> &neighbourMatrixES, unsigned int n, double ao0, double ai0, double ksi, double sar, unsigned int d, double w,gsl_rng  *r)
{

  getNeighbourMatrix(neighbourMatrix,n,1);
  getNeighbourMatrix(neighbourMatrixES,n,d);
  initializeLandscape(landscape,neighbourMatrix,n,ao0,ai0,w,r);
  getNaturalConnectedComponents(naturalComponents,landscape);
  getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrix,landscape,sar);
  getAgriculturalProduction(agriculturalProduction, landscape, ecosystemServices, ksi);
  initializePopulation(population,agriculturalProduction);

  return;

}

////////////////////////////////////////////////////////////////////////////////
// 5- ODEs and solver:
//       - populationEquation
//       - rungeKutta4
////////////////////////////////////////////////////////////////////////////////

double populationEquation(double population, double agriculturalProduction)
{
  /*
  returns the expression of the population ODE
  */

  return population*(1-population/agriculturalProduction);

}

void rungeKutta4(vector<double> &population, vector<double> &agriculturalProduction, double dt)
{
  /*
  returns the actualized population after solving the ODE with runge kutta 4 method
  */

  double totalAgriculturalProduction=0;
  unsigned long ix;
  for(ix=0;ix<agriculturalProduction.size();ix++){
    totalAgriculturalProduction+=agriculturalProduction[ix];
  }

  double k1p,k2p,k3p,k4p;
  double p1,p2,p3;
  double deltaP;

  k1p=populationEquation(population[0],totalAgriculturalProduction);
  p1=population[0]+0.5*k1p*dt;

  k2p=populationEquation(p1,totalAgriculturalProduction);
  p2=population[0]+0.5*k2p*dt;

  k3p=populationEquation(p2,totalAgriculturalProduction);
  p3=population[0]+k3p*dt;

  k4p=populationEquation(p3,totalAgriculturalProduction);

  deltaP=dt*(k1p+2*k2p+2*k3p+k4p)/6;

  population[0]+=deltaP;

  return;
}
