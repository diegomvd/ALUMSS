
  vector<unsigned int> natural_patches;
  vector< vector<unsigned int> > cluster_vector;
  Node * ahc_result;
  double cutoff_distance=1.1;
  unsigned int i, j;
  unsigned int xi, xj, yi, yj;

  /*
  get the list of natural patches in the landscape
  */
  for(i=0 ; i<landscape.size() ; ++i){
    if (landscape==0){
      natural_patches.push_back(i);
    }
  }

  /*
  allocate memory for the distance matrix and build the double pointer structure
  */
  // allocate memory for each row of the distance matrix
  double **distmat=(double **)malloc(natural_patches.size()*sizeof(double *));
  // allocate memory for each column knowing that the matrix is symetrical
  for (i=0 ; i<natural_patches.size() ; ++i){
    distmat[i]=(double *)malloc(i*sizeof(double));
  }

  /*
  fill the distance matrix with the values
  */
  for(i=0 ; i<natural_patches.size() ; ++i){
    // converting 1-D coordinates to 2-D
    xi=i%n-1;
    yi=(int)i/n;
    for(j=0 ; j<i ; ++j){
      // converting 1-D coordinates to 2-D
      xj=j%n-1;
      yj=(int)j/n;
      // entering the values in the matrix
      distmat[i][j]=sqrt( (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) )
    }
  }

  /*
  performing agglomerative hierarchical clustering
  */
  ahc_result=treecluster(natural_patches.size(),2,0,0,0,0,0,'s',distmat);

  /*
  freeing the distance matrix
  */
  for (i=0 ; i<natural_patches.size() ; ++i){
    free(distmat[i]);
  }
  free distmat;

  /*
  treating cluster_data to cut clusters at distance cutoff
  */
  // first point to the end cluster_data. cluster_data points to the first element
  // of an array of pointers to a Node struct with natural_patches.size()-1 elements

  for (i=natural_patches.size()-2 ; i>=0 ; --i){
    // clusters formed by distances bigger than the max distance we chose aren't considered
    if (ahc_result[i].distance < cutoff_distance){
      j=i;
      // we traverse the tree from bigger clusters to smaller ones by examining
      // every merge
      while( ahc_result[j].left > natural_patches.size()-1 ){
        // by construction of the clustering algorithm, j is assigned to the
        // location where cluster ahc_result[j] was created
        j=ahc_result[j].left - natural_patches.size();
      }
      // here we got to a point cluster that we want to store in a vector of
      // vectors containing the indexes of patches in the same natural cluster
      cluster_vector[].push_back(ahc_resul[j].left)
      cluster_vector[].push_back(ahc_resul[j].right)

      // now we do the same for departing from the initial right subnode
      j=i;
      while( ahc_result[j].right > natural_patches.size()-1 ){
        j=ahc_result[j].right - natural_patches.size();
      }
      cluster_vector[].push_back(ahc_resul[j])

      // go to subnode left and create a vector containing all clusters with id<neleemts
      // go to subnode right and put the same info in the previously created vector

    }
  }
