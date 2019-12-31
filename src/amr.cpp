
#include "amr.hpp"

// *********************************************************************

void amr::allocate_patch(p4est_quadrant_t* q){

  if (q->p.user_data != nullptr){
    printf("q->p.user_data not null, stopping.\n");
    abort();
  }

  // allocate patch and copy into user data pointer.
  patch::Patch*  p = new patch::Patch();
  q->p.user_data = p;
 
}


// *********************************************************************

void amr::deallocate_patch(p4est_quadrant_t* q){

  if (q->p.user_data == nullptr){
    printf("q->p.user_data == null, stopping.\n");
    abort();
  }

  patch::Patch* p=static_cast<patch::Patch*>(q->p.user_data);
  delete p;
  q->p.user_data = nullptr;
}

// *********************************************************************

void amr::allocate_patches(p4est_iter_volume_info_t* info, void* user_data)
{
  allocate_patch(info->quad);  
}

// *********************************************************************

void amr::deallocate_patches(p4est_iter_volume_info_t * info,
			void *user_data){

  deallocate_patch(info->quad);

}

// *********************************************************************

void amr::nullify_user_data(p4est_t * p4est,
			    p4est_topidx_t which_tree,
			    p4est_quadrant_t * quadrant){

  quadrant->p.user_data = nullptr;

}

// *********************************************************************

amr::Forest::Forest(sc_MPI_Comm mpicomm, ftl::section_t configuration, 
		size_t user_data_size, Context& context){

  ftl::getKeyword(configuration, "[FOREST]", "n_trees_x", n_trees_x);
  ftl::getKeyword(configuration, "[FOREST]", "n_trees_y", n_trees_y);

  // setup the connectivity between the trees in the forest
  connectivity = p4est_connectivity_new_brick(n_trees_x, n_trees_y, 
					      periodic_x, periodic_y);

  // divide the forest over the mpi processes
  forest = p4est_new(mpicomm, connectivity, user_data_size, 
		     nullify_user_data, &context );

  // initialize ghost layer
  ghost = nullptr;

}

// *********************************************************************

amr::Forest::~Forest(){

  // free all patches
  p4est_iterate(forest, nullptr, nullptr, deallocate_patches, nullptr,
		nullptr);
  
  // Destroy the p4est and the connectivity structure.  
  p4est_destroy(forest);
  p4est_connectivity_destroy(connectivity);

}

// *********************************************************************

void amr::Forest::partition(int allow_for_coarsening,
			    p4est_weight_t weight_fn){

  p4est_partition(forest, allow_for_coarsening, weight_fn);
  
}
  
// *********************************************************************

void amr::Forest::iterate(p4est_ghost_t * ghost_layer,
			  void *user_data,
			  p4est_iter_volume_t iter_volume,
			  p4est_iter_face_t iter_face,
			  p4est_iter_corner_t iter_corner){

  p4est_iterate(forest, ghost_layer, user_data, iter_volume, iter_face,
		iter_corner);
};

// *********************************************************************

void amr::Forest::iterateVolume(p4est_iter_volume_t iter_volume,
				p4est_ghost_t* ghost_layer,
				void* user_data){

  p4est_iterate(forest, ghost_layer, user_data, iter_volume, nullptr,
		nullptr);
};

// *********************************************************************

void amr::Forest::refine(int refine_recursive, int maxlevel,
		   p4est_refine_t refine_fn,
		   p4est_replace_t replace_fn){

  p4est_refine_ext(forest, refine_recursive, maxlevel, refine_fn,
		   nullify_user_data, replace_fn);
};

// *********************************************************************

void amr::Forest::coarsen(int coarsen_recursive,
			  p4est_coarsen_t coarsen_fn,
			  p4est_replace_t replace_fn){
  
  const int callback_orphans = 0;
  p4est_coarsen_ext(forest, coarsen_recursive, callback_orphans, 
		    coarsen_fn, nullify_user_data, replace_fn);
};

// *********************************************************************

void amr::Forest::balance(p4est_connect_type_t btype,
			  p4est_replace_t replace_fn) {

  p4est_balance_ext(forest,btype, nullify_user_data, replace_fn);
}

// *********************************************************************

void amr::Forest::writeVTKFile(){

  p4est_vtk_write_file (forest, NULL, P4EST_STRING "_step1");
}

// *********************************************************************

void amr::Forest::create_ghost_layer(p4est_connect_type_t btype){

  if (ghost != nullptr) {
    printf("ghost not null, stopping.\n");
    abort();
  }

  ghost = p4est_ghost_new(forest, btype);
  
}

// *********************************************************************

void amr::Forest::destroy_ghost_layer(){

  if (ghost == nullptr) {
    printf("ghost == null, stopping.\n");
    abort();
  }

  p4est_ghost_destroy(ghost);
  
}

// *********************************************************************
