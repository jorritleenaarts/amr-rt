
#include "amr.hpp"

// *********************************************************************

amr::Amr::Amr(sc_MPI_Comm mpicomm, ftl::section_t configuration, 
		size_t user_data_size, Context& context){

  ftl::getKeyword(configuration, "[FOREST]", "n_trees_x", n_trees_x);
  ftl::getKeyword(configuration, "[FOREST]", "n_trees_y", n_trees_y);

  // setup the connectivity between the trees in the forest
  connectivity = p4est_connectivity_new_brick(n_trees_x, n_trees_y, periodic_x, periodic_y);

  // divide the forest over the mpi processes
  forest = p4est_new(mpicomm, connectivity, user_data_size, nullptr, &context );

}

// *********************************************************************

amr::Amr::~Amr(){


  // free all patches
  p4est_iterate(forest, nullptr, nullptr, , nullptr,
		nullptr);
  
  // Destroy the p4est and the connectivity structure.  
  p4est_destroy(forest);
  p4est_connectivity_destroy(connectivity);

}

// *********************************************************************

void amr::Amr::partitionForest(int allow_for_coarsening,
			       p4est_weight_t weight_fn){

  p4est_partition (forest, allow_for_coarsening, weight_fn);
  
}
  
// *********************************************************************

void amr::Amr::iterateForest(p4est_ghost_t * ghost_layer,
		       void *user_data,
		       p4est_iter_volume_t iter_volume,
		       p4est_iter_face_t iter_face,
			     p4est_iter_corner_t iter_corner){

  p4est_iterate(forest, ghost_layer, user_data, iter_volume, iter_face,
		iter_corner);
};

// *********************************************************************

void amr::Amr::refineForest(int refine_recursive, int maxlevel,
		   p4est_refine_t refine_fn,
		   p4est_init_t init_fn,
		   p4est_replace_t replace_fn){

  p4est_refine_ext(forest, refine_recursive, maxlevel, refine_fn,
		   init_fn, replace_fn);
};

// *********************************************************************

void amr::Amr::coarsenForest(int coarsen_recursive,
		       int callback_orphans,
		       p4est_coarsen_t coarsen_fn,
		       p4est_init_t init_fn,
		       p4est_replace_t replace_fn){
  
  p4est_coarsen_ext(forest, coarsen_recursive, callback_orphans,
		    coarsen_fn, init_fn, replace_fn);
   };

// *********************************************************************
