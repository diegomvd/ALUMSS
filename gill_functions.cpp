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

using namespace std;

////////////////////////////////////////////////////////////////////////////////

void getNeighbours(unsigned int i, unsigned int n, vector<unsigned int> &list){
  /*
  changes the vector list so that it contains the indexes of all the closest
  von neumann neighbours of cell i and considering periodic border conditions
  for a 2d landscape of size n
  */

  list.clear();

  /*
  going through all the possible cases. note: the 2d landscape is manipulated
  as a flattened 1d array
  */
  if (i==0) {list={n-1,n*(n-1),1,n}; }
  else if(i==n*n-1) {list={n*n-2,(n-1)*n-1,n*(n-1),n-1}; }
  else if(i==n-1) { list={n-2,n*n-1,0,2*n-1}; }
  else if(i==n*(n-1)) {list={n*n-1,n*(n-2),n*(n-1)+1,0}; }
  else if(i<n) {list={i-1,i+(n-1)*n,i+1,i+n}; }
  else if(i%n==0) {list={i+n-1,i-n,i+1,i+n}; }
  else if(i>(n-1)*n-1) {list={i-1,i-n,i+1,i-n*(n-1)};}
  else if( (i+1)%n==0 ) {list={i-1,i-n,i-(n-1),i+n};}
  else {list={i-1,i-n,i+1,i+n};}
  return;
}

////////////////////////////////////////////////////////////////////////////////

void getStateNeighbours(unsigned int i, unsigned int state, unsigned int n, const vector<unsigned int> &landscape, vector<unsigned int> &state_indexes){
  /*
  fills the vector state_indexes so that it contains the indexes of all the
  closest neighbours of i in a given state. the landscape is passed as a constant
  reference so that the vector cannot be modified by the function in main
  */

  /*
  getting the neighbours indexes in neighbour_list vector
  */
  vector<unsigned int> neighbours_list;
  getNeighbours(i, n, neighbours_list);
  /*
  getting the index of neighbours in the wanted state
  */
  unsigned long ix;
  for (ix=0 ; ix<neighbours_list.size() ; ++ix){
    if (landscape[neighbours_list[ix]] == state) {
      state_indexes.push_back( neighbours_list[ix] );
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void getNaturalCluster(vector<vector<int>> &natural_components, unsigned int n, const vector<unsigned int> &landscape){
  /*
  fills a vector where each member is a vector containing the indexes of all the
  natural patches belonging to the same cluster
  */

  vector<unsigned int> natural_patches;
  unsigned int manhattan_dist;
  unsigned int i, j;
  int xi, xj, yi, yj;
  unsigned int dx, dy;

  // clearing natural components vector for refilling
  natural_components.clear();

  /*
  get the list of natural patches in the landscape
  */
  for(i=0 ; i<landscape.size() ; ++i){
    if (landscape[i]==0){
      natural_patches.push_back(i);
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


    for(i=0 ; i<natural_patches.size() ; ++i){
      add_edge(i,i,G);//not sure on how to do this
      // converting 1-D coordinates to 2-D
      xi=natural_patches[i]%n-1;
      yi=(int)natural_patches[i]/n;
      for(j=0 ; j<i ; ++j){
        // converting 1-D coordinates to 2-D
        xj=natural_patches[j]%n-1;
        yj=(int)natural_patches[j]/n;
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
        manhattan_dist=dx+dy;
        if ( manhattan_dist<2 ){
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
    natural_components.resize(component.size());
    for (i=0 ; i<natural_components.size() ; ++i){
      natural_components[i].push_back(natural_patches[i]);
      natural_components[i].push_back(component[i]);
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

double getExposure2Nature(unsigned int i, unsigned int n, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components){
  /*
  returns the exposure to the wanted state of patch i. currently it is only used
  for the exposure to nature. the exposure to nature scales with biodiversity
  hence like a SAR, where the area is the total natural area in contact with
  patch i
  */

  vector<int> cluster_membership;
  vector<double> component_area;
  double state_exposure=0;

  /*
  getting the state neighbours indexes in neighbours_state vector
  */
  vector<unsigned int> neighbours_state;
  getStateNeighbours(i, 0, n, landscape, neighbours_state); // state 0 is natural

  /*
  identify cluster membership of each natural neighbour
  */
  unsigned long ix, jx;
  for (ix=0 ; ix<neighbours_state.size() ; ++ix){
    jx=0;
    while(natural_components[jx][0]!=neighbours_state[ix]){
      jx++;
    }
    cluster_membership.push_back(natural_components[jx][1]);
  }

  /*
  calculate the area of each cluster
  */
  int area;
  for(ix=0 ; ix<cluster_membership.size(); ++ix){
    area=0;
    for(jx=0 ; jx<natural_components.size() ; ++jx){
      if (cluster_membership[ix]==natural_components[jx][1]){
        area++;
      }
    }
    component_area.push_back(area);
  }

  for (ix=0 ; ix<component_area.size(); ++ix){
    state_exposure+=0.25*pow(component_area[ix]/(n*n),0.25);
  }

  return state_exposure;
}

////////////////////////////////////////////////////////////////////////////////

void getRecoveryPropensity(unsigned int n, double Tr, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &recovery_propensity){
  /*
  fills the recovery_propensity vector given the mean time to recovery Tr, and
  the landscape
  */

  /*
  traverse the landscape and calculate the recovery propensity for the degraded
  patches
  */
  double exposure_to_nature;
  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if (landscape[ix] == 1) { // if patch ix is degraded
      exposure_to_nature=getExposure2Nature(ix,n,landscape,natural_components);
      recovery_propensity.push_back( exposure_to_nature/Tr );
    }
    else{
      recovery_propensity.push_back( 0 );
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void getDegradationPropensity(unsigned int n, double Td, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, vector<double> &degradation_propensity){
  /*
  fills the degradation_propensity vector given the mean time to degradation Td,
  and the landscape
  */

  /*
  traverse the landscape and calculate the degradation propensity for the natural
  patches
  */
  double exposure_to_nature;
  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if (landscape[ix] == 0) { // if patch ix is natural
      exposure_to_nature=getExposure2Nature(ix,n,landscape,natural_components);
      degradation_propensity.push_back( 1/Td*(1-exposure_to_nature) );
    }
    else{
      degradation_propensity.push_back( 0 );
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void getActionPropensity(unsigned int n, double w, double a, double g, double consumption_deficit,const vector<unsigned int> &landscape, vector<double> &organic_propensity, vector<double> &intense_propensity){
  /*
  fills the cropping_propensity and restoring_propensity vector
  */

  // use them to store the accumulated probability in order to normalize afterwards
  int org_cum_sum = 0;
  int int_cum_sum = 0;

  /*
  first checking if there is a consumption deficit that justifys human action
  */
  if (consumption_deficit<0){ // if humans are satisfied...
    unsigned long ix;
    for(ix=0 ; ix<landscape.size() ; ++ix){
      organic_propensity.push_back( 0 ); // ... they do not transform the landscape
    }
    for(ix=0 ; ix<landscape.size() ; ++ix){
      intense_propensity.push_back( 0 ); // ... they do not transform the landscape
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

    for (ix=0 ; ix<landscape.size() ; ++ix){

      if (landscape[ix] == 0){ // if patch is natural, hence can be converted to organic

        // get the indexes of the organic neighbours
        vector<unsigned int> organic_neighbours;
        getStateNeighbours(ix, 2, n, landscape, organic_neighbours); // state 2 is organic
        // cropping probability expression is taken from bart's and dani's paper
        organic_propensity.push_back( pow( max(0.1 , (double)organic_neighbours.size() ) , w ) );
        org_cum_sum += organic_propensity.back();
        intense_propensity.push_back( 0 );

      }
      else if (landscape[ix] == 2 ){ // if patch is organic it can be intensifyed

        // get the indexes of the organic neighbours
        vector<unsigned int> intense_neighbours;
        getStateNeighbours(ix, 3, n, landscape, intense_neighbours); // state 2 is organic
        intense_propensity.push_back( pow( max(0.1 , (double)intense_neighbours.size() ) , w ) );
        int_cum_sum += intense_propensity.back();
        organic_propensity.push_back(0);

      }
      else if (landscape[ix] == 1 || landscape[ix] == 3){
        organic_propensity.push_back(0);
        intense_propensity.push_back( 0 );

      }
      else cout << "Error: getActionPropensity: landscape contains an unrecognized value\n";
    }

    /*
    now we normalize the previously obtained probabilities and get calculate the
    probability per unit time of action given the consumption deficit
    */

    if (org_cum_sum>0 && int_cum_sum>0){ // this is to avoid dividing by zero
      for (ix=0; ix<landscape.size() ; ++ix){
        organic_propensity[ix]=organic_propensity[ix]/org_cum_sum*g*consumption_deficit*a;
        intense_propensity[ix]=intense_propensity[ix]/int_cum_sum*g*consumption_deficit*(1-a);
      }
    }
    else{
      if (org_cum_sum>0){
        for (ix=0; ix<landscape.size() ; ++ix){
          organic_propensity[ix]=organic_propensity[ix]/org_cum_sum*g*consumption_deficit;
        }
      }
      else if (int_cum_sum>0){
        for (ix=0; ix<landscape.size() ; ++ix){
          intense_propensity[ix]=intense_propensity[ix]/int_cum_sum*g*consumption_deficit;
        }
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

void getAbandonmentPropensity(double Ta, const vector<unsigned int> &landscape, const vector<double> &agricultural_production, const vector<double> &maintenance_costs, vector<double> &abandonmentO_propensity, vector<double> &abandonmentI_propensity){
  /*
  fills the abandonment_propensity vector. the probability per unit time is
  proportional to the maintenance deficit, which is the difference between
  a patch's maintenance cost and its production
  */

  unsigned int local=0;

  unsigned int ix;
  double maintenance_deficit=0;

  if (local==1){

    for (ix=0; ix<landscape.size(); ++ix){
      maintenance_deficit=maintenance_costs[ix] - agricultural_production[ix];
      if (landscape[ix]==2){ // organic patch
        abandonmentI_propensity.push_back(0);
        if (maintenance_deficit>0){
          abandonmentO_propensity.push_back(maintenance_deficit/Ta);
        }
        else{
          abandonmentO_propensity.push_back(0);
        }
      }
      else if(landscape[ix]==3){ //intensive patch
        abandonmentO_propensity.push_back(0);
        if (maintenance_deficit>0){
          abandonmentI_propensity.push_back(maintenance_deficit/Ta);
        }
        else{
          abandonmentI_propensity.push_back(0);
        }
      }
      else{ // non cropped patches
        abandonmentO_propensity.push_back(0);
        abandonmentI_propensity.push_back(0);
      }
    }
  }
  else{

    double total_maintenance_cost=0;
    double total_agricultural_production=0;
    vector<double> patch_maintenance_deficit;

    // calculate global maintenance deficit
    for (ix=0; ix<landscape.size(); ++ix){
      total_maintenance_cost+=maintenance_costs[ix];
      total_agricultural_production+=agricultural_production[ix];
      patch_maintenance_deficit.push_back(maintenance_costs[ix]-agricultural_production[ix])
    }
    maintenance_deficit=total_maintenance_cost-total_agricultural_production;
    if (maintenance_deficit>0){
      for(ix=0; ix<landscape.size(); ++ix){
        if (landscape[ix]==2){ // organic
          abandonmentO_propensity.push_back(max(0.0,patch_maintenance_deficit[ix]/Ta));
          abandonmentI_propensity.push_back(0);
        }
        else if (landscape[ix]==3){//intense
          abandonmentO_propensity.push_back(0);
          abandonmentI_propensity.push_back(max(0.0,patch_maintenance_deficit[ix]/Ta));
        }
        else{//non cropped
          abandonmentO_propensity.push_back(0);
          abandonmentI_propensity.push_back(0);
        }
      }
    }
    else{ // no maintenance deficit
      for(ix=0; ix<landscape.size(); ++ix){
        abandonmentO_propensity.push_back(0);
        abandonmentI_propensity.push_back(0);
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

void getPropensityVector(vector<double> &propensity_vector, unsigned int n, double Tr, double Td, double w, double a, double g, double Ta, double consumption_deficit, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components, const vector<double> &agricultural_production, const vector<double> &maintenance_costs){
  /*
  calls all the functions to calculate the propensity of each event and merges
  them in a single propensity vector that can be used to run the gillespie algo
  */

  vector<double> recovery_propensity;
  vector<double> degradation_propensity;
  vector<double> organic_propensity;
  vector<double> intense_propensity;
  vector<double> abandonmentO_propensity;
  vector<double> abandonmentI_propensity;

  getRecoveryPropensity(n,Tr,landscape,natural_components,recovery_propensity);
  getDegradationPropensity(n,Td,landscape,natural_components,degradation_propensity);
  getActionPropensity(n, w, a, g,consumption_deficit,landscape,organic_propensity, intense_propensity);
  getAbandonmentPropensity(Ta,landscape,agricultural_production,maintenance_costs,abandonmentO_propensity,abandonmentI_propensity);

  // clearing the previous propensity vector to refill it
  propensity_vector.clear();

  // making sure the vector is not empty to avoid bad behaviour in next loop
  propensity_vector.push_back(recovery_propensity[0]);
  unsigned long ix;
  // over these loops the cumulative sum of each propensity vector is added to
  // the total propensity vector
  for (ix=1 ; ix<recovery_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+recovery_propensity[ix]);
  }
  for (ix=0 ; ix<degradation_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+degradation_propensity[ix]);
  }
  for (ix=0 ; ix<organic_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+organic_propensity[ix]);
  }
  for (ix=0 ; ix<intense_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+intense_propensity[ix]);
  }
  for (ix=0 ; ix<abandonmentO_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+abandonmentO_propensity[ix]);
  }
  for (ix=0 ; ix<abandonmentI_propensity.size() ; ++ix){
    propensity_vector.push_back(propensity_vector.back()+abandonmentI_propensity[ix]);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void initializeLandscape( vector<unsigned int> &landscape, unsigned int n, double ao0, double ai0, gsl_rng  *r){
  /*
  initializes the landscape given a fraction of initial agricultural patches a0
  */

  //unsigned int number_cropped_patches = 1;
  unsigned long ix;
  unsigned int ix_land4crop;
  unsigned int ix_landscape;
  // this vector contains the indexes of all the natural patches
  vector<unsigned int> land4crop;

  // first build a completely natural landscape with n*n patches
  for (ix=0 ; ix<n*n; ++ix){
    landscape.push_back(0);
    land4crop.push_back(ix);
  }

  unsigned int nao0=(unsigned int) (ao0*n*n);
  unsigned int nai0=(unsigned int) (ai0*n*n);
  for (ix=0 ; ix<nao0 ; ++ix){
    // selecting random location for initial cropped patches
    ix_land4crop = gsl_rng_uniform_int(r, land4crop.size() );
    ix_landscape = land4crop[ix_land4crop];
    landscape[ix_landscape]=2; // cropping selected patch with organic agriculture
    land4crop.erase(land4crop.begin()+ix_land4crop); // removing the selected patch so it cannot be selected again
  }
  for (ix=0 ; ix<nai0 ; ++ix){
    // selecting random location for initial cropped patches
    ix_land4crop = gsl_rng_uniform_int(r, land4crop.size() );
    ix_landscape = land4crop[ix_land4crop];
    landscape[ix_landscape]=3; // cropping selected patch with intense agriculture
    land4crop.erase(land4crop.begin()+ix_land4crop); // removing the selected patch so it cannot be selected again
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void initializePopulation(vector<double> &population, const vector<double> &consumption, const vector<double> &agricultural_production){
  /*given an agricultural production, it sets the population at an equilibrium
  level
  */
  unsigned int ix;
  double total_agricultural_production=0;
  for(ix=0; ix<agricultural_production.size(); ++ix){
    total_agricultural_production+=agricultural_production[ix];
  }

  population.push_back(total_agricultural_production/consumption[0]);

  return;
}

////////////////////////////////////////////////////////////////////////////////

void initializeMaintenanceCosts(vector<double> &maintenance_costs, double y0, double m, const vector<unsigned int> &landscape){
  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if (landscape[ix]==2 || landscape[ix]==3){ //cropped patches
      maintenance_costs.push_back(m*y0);
    }
    else{
      maintenance_costs.push_back(0);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void getAgriculturalProduction(vector<double> &agricultural_production, unsigned int n, double y0, double k, double phi, const vector<unsigned int> &landscape, const vector<vector<int>> &natural_components){
  /*
  returns the total agricultural production for a given "landscape" and
  minimum yield "y". the minimum yield is for a cropped patch with only
  non natural neighbours. natural neighbours raise yield.
  */

  agricultural_production.clear();

  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if(landscape[ix]==2){ // cropped patches
      agricultural_production.push_back(y0*(1+phi*getExposure2Nature(ix,n,landscape,natural_components)));
    }
    else if(landscape[ix]==3){ //intense
      agricultural_production.push_back(y0*k);
    }
    else{
      agricultural_production.push_back(0);
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////

void updateMaintenanceIntense(vector<double> &maintenance_costs, double dtp, double y0, double m, double k, double Ti, const vector<unsigned int> &landscape){
  /*
  the maintenance cost of intense patches is incremented
  */

  unsigned long ix;
  for (ix=0 ; ix<landscape.size() ; ++ix){
    if(landscape[ix]==3){ // intense
      maintenance_costs[ix]+= y0*(k-m)/Ti*dtp;
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////

double populationEquation(double r0, double agricultural_production, double population, double consumption){
  /*
  returns the expression of the population ODE
  */

  return r0*population*(1-consumption*population/agricultural_production);

}

////////////////////////////////////////////////////////////////////////////////

double consumptionEquation(double kg, double kd, double minimum_consumption, double agricultural_production, double population, double consumption){
  /*
  returns the expression of the consumption ODE, the growing and decreasing rates
  are different.
  */
  double res;

  if( consumption < agricultural_production/population){
    res = (agricultural_production/population - consumption)*kg ;
  }
  else if(consumption> agricultural_production/population && consumption>minimum_consumption){
    res = (agricultural_production/population - consumption)*kd*kg ;
  }
  else{
    res = 0;
  }

  return res;
}

////////////////////////////////////////////////////////////////////////////////

void RungeKutta4(double r0, double kg, double kd, double minimum_consumption, double dt, double agricultural_production, vector<double> &population, vector<double> &consumption){
  /*
  returns the actualized population after solving the ODE with runge kutta 4 method
  */

  double k1p,k2p,k3p,k4p,k1c,k2c,k3c,k4c;
  double p1,p2,p3,c1,c2,c3;
  double deltaP, deltaC;

  k1p=populationEquation(r0,agricultural_production,population[0],consumption[0]);
  k1c=consumptionEquation(kg,kd,minimum_consumption,agricultural_production,population[0],consumption[0]);
  p1=population[0]+0.5*k1p*dt;
  c1=consumption[0]+0.5*k1c*dt;

  k2p=populationEquation(r0,agricultural_production,p1,c1);
  k2c=consumptionEquation(kg,kd,minimum_consumption,agricultural_production,p1,c1);
  p2=population[0]+0.5*k2p*dt;
  c2=consumption[0]+0.5*k2c*dt;

  k3p=populationEquation(r0,agricultural_production,p2,c2);
  k3c=consumptionEquation(kg,kd,minimum_consumption,agricultural_production,p2,c2);
  p3=population[0]+k3p*dt;
  c3=consumption[0]+k3c*dt;

  k4p=populationEquation(r0,agricultural_production,p3,c3);
  k4c=consumptionEquation(kg,kd,minimum_consumption,agricultural_production,p3,c3);

  deltaP=dt*(k1p+2*k2p+2*k3p+k4p)/6;
  deltaC=dt*(k1c+2*k2c+2*k3c+k4c)/6;

  population[0]+=deltaP;
  consumption[0]+=deltaC;

  return;
}
