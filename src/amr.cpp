
#include "amr.hpp"

// *********************************************************************

amr::Amr::Amr(sc_MPI_Comm mpicomm, ftl::section_t configuration, 
		size_t user_data_size, Context& context){

  ftl::getKeyword(configuration, "[FOREST]", "n_trees_x", n_trees_x);
  ftl::getKeyword(configuration, "[FOREST]", "n_trees_y", n_trees_y);

  // setup the connectivity between the trees in the forest
  connectivity = p4est_connectivity_new_brick(n_trees_x, n_trees_y, periodic_x, periodic_y);

  // divide the forest over the mpi processes
  forest = p4est_new(mpicomm, connectivity, user_data_size, NULL, &context );

}

// *********************************************************************

amr::Amr::~Amr(){

  // Destroy the p4est and the connectivity structure. 
  p4est_destroy(forest);
  p4est_connectivity_destroy(connectivity);

}

// *********************************************************************

void amr::Amr::refineForest(int refine_recursive, 
		       p4est_refine_t refine_fn, p4est_init_t init_fn){

  p4est_refine(forest, refine_recursive, refine_fn, init_fn);

}

// *********************************************************************
