#include "gill_functions.h"

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
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
void getNeighbourMatrix(vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, double d)
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
//       - updateNCCadding
//       - updateNCCremoving
//       - getEcosystemServiceProvision
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
  vector<int> newNaturalComponent, newNaturalComponent2;
  newNaturalComponent.push_back(i);

  if(neighboursNatural.size()==0){ //no natural neighbour
      naturalComponents.push_back(newNaturalComponent); // add it to the list
  }
  else{

    vector<unsigned int>::iterator it1;
    vector<vector<int>>::iterator it2;
    vector<int>::iterator it3;
    vector<vector<vector<int>>::iterator> toErase;

    /*
    first traversing the natural components and then the neighbours guarantees that
    pointers to components are located in a sorted way. this is key to ensure
    that the erasing process doesn't mess up with the memory.
    */
    unsigned int neighboursFound = 0;
    for(it2=naturalComponents.begin();it2!=naturalComponents.end();it2++){ // traverse all the components
      for(it1=neighboursNatural.begin();it1!=neighboursNatural.end();it1++){ // traverse all the natural neigbhours
        if( find( it2->begin(), it2->end(), *it1) != it2->end() ){ // found a neighbour in this component
          neighboursFound+=1;
          if ( find( toErase.begin(), toErase.end(), it2 ) == toErase.end() ){
            toErase.push_back(it2);
          }
        }
      }
      if(neighboursFound==neighboursNatural.size()){ // end the search if all neighbours were located
        break;
      }
    }


    // if there is a single component in toErase just add i to that component, no need of merging
    if (toErase.size()>0 && toErase.size()<2){
      toErase[0]->push_back(i);
    }
    else if(toErase.size()>1){ // if there are more, erase them and push back the merged ones
      // erasing components that are going to be merged
      long ix;
      for(ix=0;ix<toErase.size();ix++){
        for(it3=toErase[ix]->begin();it3!=toErase[ix]->end();it3++){
          newNaturalComponent.push_back(*it3);
        }
      }
      /*now erase the components: traverse erase vector backwards to be sure of
      addressing the correct bits of memory erasing first the furthest pointers*/
      for(ix=toErase.size()-1;ix>=0;ix--){
        naturalComponents.erase(toErase[ix]);
      }
      /*now add the new natural component*/
      naturalComponents.push_back(newNaturalComponent);
    }
    else{ // toErase is empty, in which case there is an error in the code
      cout << "Error: toErase size is " << toErase.size() << " but toErase cannot be empty\n";
    }

  }

  return;
}

void updateNCCremoving(vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape, int l)
{
  //cout <<"\n";
  unsigned int n = (unsigned int) sqrt(landscape.size());

  // find cluster of cell l
  unsigned long ix;
  // these iterators are to traverse the naturalComponents
  vector<vector<int>>::iterator it1;
  vector<int>::iterator it2;

  // iterating over natural components and checking to which natural component
  // belonged the natural cell that needs to be removed.

  //this iterator is to store the component to which the removed cell belonged
  vector<vector<int>>::iterator itComp;
  // iterate over components
  for(it1=naturalComponents.begin();it1!=naturalComponents.end();it1++){
    // find whether natural cell l is in this natural component
    it2=find(it1->begin(), it1->end(),l);
    // if it is then erase it and store the component
    if(it2!=it1->end()){
      it1->erase(it2);
      itComp =it1;
      break;
    }
  }

  // fill natural patches to get connected natural components
  vector<unsigned int> naturalPatches;
  for (it2=itComp->begin();it2!=itComp->end();it2++){
    naturalPatches.push_back(*it2);
  }


  // erase concerned cluster from naturalComponents
  naturalComponents.erase(itComp);

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
  unsigned long i,ix,jx;
  unsigned int nNeighbours;

  nNeighbours = neighbourMatrix[0].size();

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
          area=(double)naturalComponents[jx].size()/landscape.size();
          ecosystemServiceProvision+=(double) pow(area,sar)/nNeighbours;
          break;
        }
      }
    }

    ecosystemServices.push_back(ecosystemServiceProvision);
    neighboursState.clear();
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// 3- Calculation of events' propensities:
//       - esSaturationFunction
//       - getAgriculturalProduction
//       - getConsumptionDeficit
//       - getSpontaneousPropensity
//       - getAgroPropensity
//       - getAbandonmentPropensity
//       - getPropensityVector
////////////////////////////////////////////////////////////////////////////////

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
      // putting baseline production 0.5 as a test...
      agriculturalProduction.push_back( 0.2 + ecosystemServices[ix] ) ;
    }
    else if(landscape[ix]==3){ //intense
      agriculturalProduction.push_back( ksi );
    }
    else{
      agriculturalProduction.push_back(0);
    }
  }

  return;
}

double getConsumptionDeficit(const vector<double> &agriculturalProduction, const vector<double> &population)
{
  double totalAgriculturalProduction=0;
  double consumptionDeficit;
  unsigned long ix;
  for(ix=0;ix<agriculturalProduction.size();ix++){
    totalAgriculturalProduction+=agriculturalProduction[ix];
  }
  if(population[0]>0){
    consumptionDeficit=population[0] - totalAgriculturalProduction;
  }
  else{
    consumptionDeficit=0;
  }

  return consumptionDeficit;
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
      recoveryPropensity.push_back( 1/Tr * ecosystemServices[ix] );
      degradationPropensity.push_back(0);
    }
    else if(landscape[ix] == 0){
      recoveryPropensity.push_back( 0 );
      degradationPropensity.push_back( 1/Td*( 1 - ecosystemServices[ix] ));
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

    if (expansionCumSum>0 && intenseCumSum>0){ // this is to avoid dividing by zero
      for (ix=0; ix<landscape.size() ; ++ix){
        expansionPropensity[ix] = expansionPropensity[ix] / expansionCumSum / Tag * consumptionDeficit * (1-a);
        intensePropensity[ix] = intensePropensity[ix] / intenseCumSum / Tag * consumptionDeficit * a;
      }
    }
    else{
      if (expansionCumSum>0){
        for (ix=0; ix<landscape.size() ; ++ix){
          expansionPropensity[ix] = expansionPropensity[ix] / expansionCumSum / Tag * consumptionDeficit;
        }
      }
      else if (intenseCumSum>0){
        if(a>0){
          for (ix=0; ix<landscape.size() ; ++ix){
            intensePropensity[ix] = intensePropensity[ix] / intenseCumSum / Tag * consumptionDeficit;
          }
        }
        else{ // this makes that in the case a=0 there is no intensification
          for (ix=0; ix<landscape.size() ; ++ix){
            intensePropensity[ix] = 0;
          }
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
      naturalAbandonPropensity.push_back( 1/Tab * ( 1 - ecosystemServices[ix] ) );
    }
    else if(landscape[ix]==3){ //intensive patch
      degradedAbandonPropensity.push_back( 1/Tab * ( 1 - ecosystemServices[ix] ) );
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

void initializeLandscape( vector<unsigned int> &landscape, const vector<vector<unsigned int>> &neighbourMatrix, unsigned int n, double a0, double d0, double a, double w, gsl_rng  *r)
{
  /*
  initializes the landscape given a fraction of initial agricultural patches a0 and degraded patches d0
  */

  //unsigned int number_cropped_patches = 1;
  unsigned long ix,jx,lx;
  // this vector contains the indexes of all the natural patches
  vector<double> probConversion;
  vector<double> probIntense;
  vector<double> probDegradation;

  unsigned int nao0=(unsigned int) (a0*n*n*(1-a));
  unsigned int nai0=(unsigned int) (a0*n*n*a);
  unsigned int na0=nao0+nai0;

  // first build a completely natural landscape with n*n patches
  landscape.push_back(0);
  probConversion.push_back(1);

  for (ix=1 ; ix<n*n; ++ix){
    landscape.push_back(0);
    probConversion.push_back(probConversion.back()+1);
    probIntense.push_back(0);
  }
  double cumprob;
  vector<unsigned int> organicNeighbours;
  for(ix=0;ix<na0;++ix){
    jx=0;
    cumprob = gsl_rng_uniform(r)*probConversion.back();
    while (cumprob > probConversion[jx]){
      jx++;
    }
    landscape[jx]=2;
    probIntense[jx]=1;

    // recalculating probconversion
    if (landscape[0]==0){
      organicNeighbours.clear();
      getNeighboursState(organicNeighbours,neighbourMatrix,landscape,0, 2);
      probConversion[0]=pow( max(0.1 , (double)organicNeighbours.size() ) , w ) ;
    }
    else{
      probConversion[0]=0;
    }
    for(lx=1;lx<probConversion.size();lx++){
      if (landscape[lx]==0){
        organicNeighbours.clear();
        getNeighboursState(organicNeighbours,neighbourMatrix,landscape,lx, 2);
        probConversion[lx]=probConversion[lx-1]+pow( max(0.1 , (double)organicNeighbours.size() ) , w ) ;
      }
      else{
        probConversion[lx]=probConversion[lx-1];
      }
    }
  }

  // perform the cumulative sum of probIntense
  for(ix=1; ix<probIntense.size(); ix++){
    probIntense[ix]=probIntense[ix]+probIntense[ix-1];
  }
  vector<unsigned int> intenseNeighbours;
  for(ix=0;ix<nai0;++ix){
    jx=0;
    while (gsl_rng_uniform(r)*probIntense.back()>probIntense[jx]){
      jx++;
    }
    landscape[jx]=3;

    // recalculating probintense
    // first initialize the vector to do the cumulative sum
    if (landscape[0]==2){
      intenseNeighbours.clear();
      getNeighboursState(intenseNeighbours,neighbourMatrix,landscape,0, 3);
      probIntense[0]=pow( max(0.1 , (double)intenseNeighbours.size() ) , w ) ;
    }
    else{
      probIntense[0]=0;
    }
    for(lx=1;lx<probIntense.size();lx++){
      if (landscape[lx]==2){
        intenseNeighbours.clear();
        getNeighboursState(intenseNeighbours,neighbourMatrix,landscape,lx, 3);
        probIntense[lx]=probIntense[lx-1]+pow( max(0.1 , (double)intenseNeighbours.size() ) , w ) ;
      }
      else{
        probIntense[lx]=probIntense[lx-1];
      }
    }
  }

  // now deal with degraded patches
  unsigned int nd0=(unsigned int) (d0*n*n);
  // initialize probDegradation
  if (landscape[0]==0){
    probDegradation.push_back(1);
  }
  else{
    probDegradation.push_back(0);
  }

  // fill probDegradation for first pick
  for (ix=1; ix<landscape.size(); ix++){
    if (landscape[ix]==0){
      probDegradation.push_back(probDegradation.back()+1);
    }
    else{
      probDegradation.push_back(probDegradation.back());
    }
  }

  // now start filling the landscape with the degraded cells
  vector<unsigned int> degradedNeighbours;
  for (ix=0; ix<nd0; ++ix){
    // select to be degraded cell
    jx=0;
    cumprob = gsl_rng_uniform(r)*probDegradation.back();
    while (cumprob > probDegradation[jx]){
      jx++;
    }
    landscape[jx]=1;

    // update the probDegradation
    // first initialize the vector
    if (landscape[0]==0){
      degradedNeighbours.clear();
      getNeighboursState(degradedNeighbours,neighbourMatrix,landscape,0, 1);
      probDegradation[0]=pow( max(0.1 , (double)degradedNeighbours.size() ) , w ) ;
    }
    else{
      probDegradation[0]=0;
    }
    for(lx=1;lx<probDegradation.size();lx++){
      if (landscape[lx]==0){
        degradedNeighbours.clear();
        getNeighboursState(degradedNeighbours,neighbourMatrix,landscape,lx, 1);
        probDegradation[lx]=probDegradation[lx-1]+pow( max(0.1 , (double)degradedNeighbours.size() ) , w ) ;
      }
      else{
        // cout << "lx = " << lx << "\n";
        probDegradation[lx]=probDegradation[lx-1];
      }
    }
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

void initializeSES( vector<unsigned int> &landscape, vector<double> &population, vector<vector<int>> &naturalComponents, vector<double> &agriculturalProduction, vector<double> &ecosystemServices, vector<vector<unsigned int>> &neighbourMatrix, vector<vector<unsigned int>> &neighbourMatrixES, unsigned int n, double a0, double d0, double a, double ksi, double sar, double w, gsl_rng  *r)
{

  initializeLandscape(landscape,neighbourMatrix,n,a0,d0,a,w,r);
  getNaturalConnectedComponents(naturalComponents,landscape);
  getEcosystemServiceProvision(ecosystemServices,naturalComponents,neighbourMatrixES,landscape,sar);
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

////////////////////////////////////////////////////////////////////////////////
// 6- Outputs:
//       - saveAggregated
//       - saveLandscape
//       - saveComponents
////////////////////////////////////////////////////////////////////////////////

void saveAggregated(ofstream &file, double t, const vector<double> &population, const vector<unsigned int> &landscape, const vector<double> &agriculturalProduction, const vector<vector<int>> &naturalComponents, const vector<double> &ecosystemServices, unsigned int nn, double ripleyDistance)
{
  unsigned long numComponents = naturalComponents.size();
  unsigned long ix;
  double meanSize=0;
  double squaredMeanSize=0;
  double maxSize=0;
  double stdSize;
  double componentSize;

  for (ix=0;ix<naturalComponents.size();ix++){
    componentSize=(double)naturalComponents[ix].size()/ecosystemServices.size();
    meanSize+=componentSize;
    squaredMeanSize+=componentSize*componentSize;
    if(componentSize>maxSize){
      maxSize=componentSize;
    }
  }
  if (numComponents>0){
    meanSize/=numComponents;
    squaredMeanSize/=numComponents;
  }
  stdSize=squaredMeanSize-meanSize*meanSize;
  if(stdSize<0){
    stdSize=0;
  }
  else{
    stdSize=sqrt(stdSize);
  }

  double meanES=0;
  double squaredMeanES=0;
  double stdES;
  for (ix=0;ix<ecosystemServices.size();ix++){
    meanES+=ecosystemServices[ix];
    squaredMeanES+=ecosystemServices[ix]*ecosystemServices[ix];
  }
  meanES/=ecosystemServices.size();
  squaredMeanES/=ecosystemServices.size();
  stdES=sqrt(squaredMeanES-meanES*meanES);

  /*
  Second part is a copy of saveAggregated
  */
  double n=0,d=0,a0=0,a1=0;

  for(ix=0;ix<landscape.size();ix++){
    if (landscape[ix]==0){
      n+=1;
    }
    else if (landscape[ix]==1){
      d+=1;
    }
    else if (landscape[ix]==2){
      a0+=1;
    }
    else{
      a1+=1;
    }
  }
  n/=landscape.size();d/=landscape.size();a0/=landscape.size();a1/=landscape.size();

  double totalY;
  for(ix=0; ix<agriculturalProduction.size(); ++ix){
    totalY+=agriculturalProduction[ix];
  }


  /*
  Last part is a copy of saveRipley
  */

  // first get the neighbour matrix given a ripley distance, here it is 1
  vector<vector<unsigned int>> neighbourMatrix;
  getNeighbourMatrix(neighbourMatrix,nn,ripleyDistance);

  // then determine the number of points for each type of land
  double nN=0,nD=0,nA0=0,nA1=0;
  for(ix=0;ix<landscape.size();ix++){
    if (landscape[ix]==0){
      nN+=1;
    }
    else if (landscape[ix]==1){
      nD+=1;
    }
    else if (landscape[ix]==2){
      nA0+=1;
    }
    else{
      nA1+=1;
    }
  }

  // now calculate ripley without normalizing
  double ripleyN=0, ripleyD=0, ripleyA0=0, ripleyA1=0;
  vector<unsigned int> stateNeighbours;
  for(ix=0;ix<landscape.size();ix++){
    stateNeighbours.clear();
    if (landscape[ix]==0){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 0);
      ripleyN+=stateNeighbours.size();
    }
    else if (landscape[ix]==1){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 1);
      ripleyD+=stateNeighbours.size();
    }
    else if (landscape[ix]==2){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 2);
      ripleyA0+=stateNeighbours.size();
    }
    else{
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 3);
      ripleyA1+=stateNeighbours.size();
    }
  }

  if (nN>0){ripleyN *= landscape.size()/(nN*nN);}
  if (nD>0){ripleyD *= landscape.size()/(nD*nD);}
  if (nA0>0){ripleyA0 *= landscape.size()/(nA0*nA0);}
  if (nA1>0){ripleyA1 *= landscape.size()/(nA1*nA1);}

  file << t << " " << population[0] << " " << n << " " << d << " " << a0 << " " << a1 << " " << totalY << " " << numComponents << " " << meanSize << " " << stdSize << " " << maxSize << " " << meanES << " " << stdES << " " << ripleyN << " " << ripleyD << " " << ripleyA0 << " " << ripleyA1 << "\n";

  return;
}

void saveLandscape(ofstream &file, double t, const vector<unsigned int> &landscape)
{
  unsigned long ix;

  file << t;
  for (ix=0;ix<landscape.size();ix++){
    file << " " << landscape[ix];
  }
  file << "\n";

  return;
}

void saveComponents(ofstream &file, double t, const vector<unsigned int> &landscape, const vector<vector<int>> &naturalComponents)
{
  unsigned long ix,jx;
  unsigned int test;
  vector<int>::const_iterator it;

  file << t;
  for (ix=0;ix<landscape.size();ix++){
    test=0;
    for (jx=0;jx<naturalComponents.size();jx++){
      it = find( naturalComponents[jx].begin(), naturalComponents[jx].end(), ix );
      if (it!=naturalComponents[jx].end()){
        file << " " << jx;
        test = 1;
        break;
      }
    }
    if (test==0){
      file << " " << nan("");
    }
  }

  file << "\n";

  return;
}

void saveSensitivityOutput(ofstream &file, unsigned int nn, double ripleyDistance, const vector<double> &population, const vector<vector<int>> &naturalComponents, const vector<unsigned int> &landscape, const vector<double> &ecosystemServices)
{

  /*
  First part is a copy of the saveLandMetrics. Of course this can be improved to
  avoid code repetition
  */
  unsigned long numComponents = naturalComponents.size();
  unsigned long ix;
  double meanSize=0;
  double squaredMeanSize=0;
  double maxSize=0;
  double stdSize;
  double componentSize;

  for (ix=0;ix<naturalComponents.size();ix++){
    componentSize=(double)naturalComponents[ix].size()/ecosystemServices.size();
    meanSize+=componentSize;
    squaredMeanSize+=componentSize*componentSize;
    if(componentSize>maxSize){
      maxSize=componentSize;
    }
  }
  meanSize/=numComponents;
  squaredMeanSize/=numComponents;
  stdSize=sqrt(squaredMeanSize-meanSize*meanSize);

  double meanES=0;
  double squaredMeanES=0;
  double stdES;
  for (ix=0;ix<ecosystemServices.size();ix++){
    meanES+=ecosystemServices[ix];
    squaredMeanES+=ecosystemServices[ix]*ecosystemServices[ix];
  }
  meanES/=ecosystemServices.size();
  squaredMeanES/=ecosystemServices.size();
  stdES=sqrt(squaredMeanES-meanES*meanES);

  /*
  Second part is a copy of saveAggregated
  */
  double n=0,d=0,a0=0,a1=0;

  for(ix=0;ix<landscape.size();ix++){
    if (landscape[ix]==0){
      n+=1;
    }
    else if (landscape[ix]==1){
      d+=1;
    }
    else if (landscape[ix]==2){
      a0+=1;
    }
    else{
      a1+=1;
    }
  }
  n/=landscape.size();d/=landscape.size();a0/=landscape.size();a1/=landscape.size();

  /*
  Last part is a copy of saveRipley
  */

  // first get the neighbour matrix given a ripley distance, here it is 1
  vector<vector<unsigned int>> neighbourMatrix;
  getNeighbourMatrix(neighbourMatrix,nn,ripleyDistance);

  // then determine the number of points for each type of land
  double nN=0,nD=0,nA0=0,nA1=0;
  for(ix=0;ix<landscape.size();ix++){
    if (landscape[ix]==0){
      nN+=1;
    }
    else if (landscape[ix]==1){
      nD+=1;
    }
    else if (landscape[ix]==2){
      nA0+=1;
    }
    else{
      nA1+=1;
    }
  }

  // now calculate ripley without normalizing
  double ripleyN=0, ripleyD=0, ripleyA0=0, ripleyA1=0;
  vector<unsigned int> stateNeighbours;
  for(ix=0;ix<landscape.size();ix++){
    stateNeighbours.clear();
    if (landscape[ix]==0){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 0);
      ripleyN+=stateNeighbours.size();
    }
    else if (landscape[ix]==1){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 1);
      ripleyD+=stateNeighbours.size();
    }
    else if (landscape[ix]==2){
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 2);
      ripleyA0+=stateNeighbours.size();
    }
    else{
      getNeighboursState(stateNeighbours, neighbourMatrix, landscape, ix, 3);
      ripleyA1+=stateNeighbours.size();
    }
  }

  if (nN>0){ripleyN *= landscape.size()/(nN*nN);}
  if (nD>0){ripleyD *= landscape.size()/(nD*nD);}
  if (nA0>0){ripleyA0 *= landscape.size()/(nA0*nA0);}
  if (nA1>0){ripleyA1 *= landscape.size()/(nA1*nA1);}

  file << numComponents << " " << meanSize << " " << stdSize << " " << maxSize << " " << meanES << " " << stdES << " " << population[0] << " " << n << " " << d << " " << a0 << " " << a1 << " " << ripleyN << " " << ripleyD << " " << ripleyA0 << " " << ripleyA1 <<"\n";

  return;
}
