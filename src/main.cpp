
// standard libraries
#include <iostream>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <cmath>
#include <cassert>

// p4est libraries

// local includes
#include "amr.hpp"
#include "patch.hpp"
#include "ftools.hpp"
#include "context.hpp"
#include "hdf.hpp"
#include "math.hpp"

// *********************************************************************

void initialize_patch(p4est_quadrant_t* q, Context* context){
  /*
   * Allocate space to hold the user data in the patch connected to
   * quadrant q.
   */

  patch::Patch* p = static_cast<patch::Patch*>(q->p.user_data);

   //compute pixel coordinates of lower left corner in input atmos
  int xmin = (q->x * context->nx_in_atmos) / P4EST_ROOT_LEN;
  int ymin = (q->y * context->ny_in_atmos) / P4EST_ROOT_LEN;

  int xmax = xmin +
    (P4EST_QUADRANT_LEN(q->level) * context->nx_in_atmos)/P4EST_ROOT_LEN - 1;
   int ymax = ymin +
    (P4EST_QUADRANT_LEN(q->level) * context->ny_in_atmos)/P4EST_ROOT_LEN - 1;

   int ng = context->n_guard_zones;
   
   p->xs = xmin;
   p->xe = xmax;
   p->ys = ymin;
   p->ye = ymax;

   p->xr = context->nx_per_patch;
   p->yr = context->ny_per_patch;

   p->xrb = context->nx_per_patch + 2 * ng;
   p->yrb = context->ny_per_patch + 2 * ng;
 
   p->density.resize(-ng,p->xr,-ng,p->yr);
   p->density = 0.0;
   
}

// *********************************************************************

int initial_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{
    return 1;
}

// *********************************************************************

int density_refine_condition(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{

  Context* context = static_cast<Context*>(p4est->user_pointer);
  patch::Patch* p = static_cast<patch::Patch*>(quadrant->p.user_data);

  float minimum_density = 1e-13;
  float maximum_density = 1e-6;
  

  // first check whether patch falls below minimum density requirement
  float rmax = 1e-30;
    for (int j = 0; j < p->yr; ++j ){
      for (int i = 0; i < p->xr; ++i ){
	rmax = max(rmax, p->density(j,i));
      }
    }

    if (rmax < minimum_density) {
      // density is too low in patch so we don't want to refine
      return 0;
    }

    // then check whether all children fall above maximum density requirement
    float rmin = 1e30;
    for (int j = 0; j < p->yr; ++j ){
      for (int i = 0; i < p->xr; ++i ){
	rmin = min(rmin, p->density(j,i));
      }
    }
  if (rmin > maximum_density) {
    // density is too high in patch so we don't want to refine
    return 0;
  }



  const float density_gradient_threshold = 0.2;

  mem::Array<float,3> grad  = math::grid_gradient(p->density);

  float g0 = 0.0;
  float g1 = 0.0;
  for (int j = 0; j < p->yr; ++j ){
    for (int i = 0; i < p->xr; ++i ){
      g0 =  fabs(grad(j,i,0) / p->density(j,i));
      g1 =  fabs(grad(j,i,1) / p->density(j,i));
          if ( g0 >= density_gradient_threshold || 
	       g1 >= density_gradient_threshold ){	 
	return 1;
      }
    }
  }

  return 0;

}


// *********************************************************************

int density_coarsen_condition(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrants[])
{

  Context* context = static_cast<Context*>(p4est->user_pointer);
  const float density_gradient_threshold = 0.02;


  float minimum_density = 1e-13;
  float maximum_density = 1e-6;
  
  // first check whether all children fall below minimum density requirement
  float rmax = 1e-30;
  for (int k = 0; k < P4EST_CHILDREN; k++) {
    patch::Patch* p = static_cast<patch::Patch*>(quadrants[k]->p.user_data);
    for (int j = 0; j < p->yr; ++j ){
      for (int i = 0; i < p->xr; ++i ){
	rmax = max(rmax, p->density(j,i));
      }
    }
  }
  if (rmax < minimum_density) {
    patch::Patch* p = static_cast<patch::Patch*>(quadrants[0]->p.user_data);
    //    printf("level %i: (%i, %i) coarsened (min density).\n", 
    //	   quadrants[0]->level, p->xs, p->ys);
    return 1;
  }

  // then check whether all children fall above maximum density requirement
  float rmin = 1e30;
  for (int k = 0; k < P4EST_CHILDREN; k++) {
    patch::Patch* p = static_cast<patch::Patch*>(quadrants[k]->p.user_data);
    for (int j = 0; j < p->yr; ++j ){
     for (int i = 0; i < p->xr; ++i ){
       rmin = min(rmin, p->density(j,i));
     }
   }
  }
  if (rmin > maximum_density) {
    patch::Patch* p = static_cast<patch::Patch*>(quadrants[0]->p.user_data);
    // printf("level %i: (%i, %i) coarsened (max density).\n", 
    //	   quadrants[0]->level, p->xs, p->ys);
    return 1;
  }

  // minmax did not flag for coarsening. now check for gradients.
  for (int k = 0; k < P4EST_CHILDREN; k++) {

    patch::Patch* p = static_cast<patch::Patch*>(quadrants[k]->p.user_data);
    mem::Array<float,3> grad  = math::grid_gradient(p->density);

    float g0 = 0.0;
    float g1 = 0.0;
    for (int j = 0; j < p->yr; ++j ){
      for (int i = 0; i < p->xr; ++i ){
	g0 =  fabs(grad(j,i,0) / p->density(j,i));
	g1 =  fabs(grad(j,i,1) / p->density(j,i));
	if ( g0 >= density_gradient_threshold || 
	     g1 >= density_gradient_threshold ){
	  // printf("quadrant at level %i is retained as is.\n", quadrants[0]->level);	
	  return 0;
	}
      }
    }  
  } 
  
  patch::Patch* p = static_cast<patch::Patch*>(quadrants[0]->p.user_data);
  //printf("level %i: (%i, %i) coarsened (gradient).\n",
  //	 quadrants[0]->level, p->xs, p->ys);
  return 1;

}

// *********************************************************************

void refine_patch_function(p4est_t * p4est,
			    p4est_topidx_t which_tree,
			    int num_outgoing,
			    p4est_quadrant_t * outgoing[],
			    int num_incoming,
			    p4est_quadrant_t * incoming[]){

  Context* context = static_cast<Context*>(p4est->user_pointer);

  // coarse user data
  patch::Patch* oldp = static_cast<patch::Patch*>(outgoing[0]->p.user_data);
  
  //allocate and initialize child patches
  assert(P4EST_CHILDREN == num_incoming);
  patch::Patch* newp;
  for (int i = 0; i < P4EST_CHILDREN; i++) {
    amr::allocate_patch(incoming[i]);
    initialize_patch(incoming[i], context);
    newp = static_cast<patch::Patch*>(incoming[i]->p.user_data);
    newp->inject_and_interpolate_from_coarse(oldp);
    
  }

  // deallocate old patch
  amr::deallocate_patch(outgoing[0]);

}

// *********************************************************************

void coarsen_patch_function(p4est_t * p4est,
			    p4est_topidx_t which_tree,
			    int num_outgoing,
			    p4est_quadrant_t * outgoing[],
			    int num_incoming,
			    p4est_quadrant_t * incoming[]){

  Context* context = static_cast<Context*>(p4est->user_pointer);

  patch::Patch* oldp[4];

  // allocate coarse user data
  assert(num_incoming == 1);
  amr::allocate_patch(incoming[0]);
  initialize_patch(incoming[0], context);
  patch::Patch* newp = static_cast<patch::Patch*>(incoming[0]->p.user_data);

  assert(num_outgoing == P4EST_CHILDREN);

  for (int i = 0; i < P4EST_CHILDREN; i++) { 
    patch::Patch* oldp = static_cast<patch::Patch*>(outgoing[i]->p.user_data);
    newp->average_from_fine(oldp);    
  }

  for (int i = 0; i < P4EST_CHILDREN; i++) {
    amr::deallocate_patch(outgoing[i]);
  }

}

// *********************************************************************

void replace_patch(p4est_t * p4est,
		   p4est_topidx_t which_tree,
		   int num_outgoing,
		   p4est_quadrant_t * outgoing[],
		   int num_incoming,
		   p4est_quadrant_t * incoming[]){

 if (num_outgoing > 1) {
   // coarsening
   coarsen_patch_function(p4est, which_tree, num_outgoing, outgoing, 
			  num_incoming, incoming);
 } else {
   //refinement
   refine_patch_function(p4est, which_tree, num_outgoing, outgoing, 
			 num_incoming, incoming);
 }
}

// *********************************************************************

void initialize_new_patches(p4est_iter_volume_info_t* info, void* user_data)
  // callback function following p4est_iter_volume_t prototype
{

  p4est_quadrant_t* q = info->quad;
  Context* context = static_cast<Context *>(info->p4est->user_pointer);
  patch::Patch* p = static_cast<patch::Patch*>(q->p.user_data);

  initialize_patch(q, context);
 
}

// *********************************************************************

void fill_boundaries(p4est_iter_volume_info_t* info, void* user_data)
  // callback function following p4est_iter_volume_t prototype
{

  p4est_quadrant_t* q = info->quad;
  patch::Patch* p = static_cast<patch::Patch*>(q->p.user_data);

  // left face
  for (int j = 0; j < p->yr; ++j)  {
    p->density(j, -1) = p->density(j, 0);
  }
  // right face
  for (int j = 0; j < p->yr; ++j) {
      p->density(j, p->xr) = p->density(j, p->xr-1);
  }
  // bottom face
  for (int i = 0; i < p->xr; ++i) {
    p->density(-1, i) = p->density(0, i);
  }
  //top face
  for (int i = 0; i < p->xr; ++i)  {
    p->density(p->yr, i) = p->density(p->yr-1, i);
  }

  //and the four corners
  p->density(-1, -1) = p->density(0, 0);
  p->density(-1, p->xr) = p->density(0, p->xr-1);
  p->density(p->yr, p->xr) = p->density(p->yr-1, p->xr-1);
  p->density(p->yr, -1) = p->density(p->yr-1, 0);

}

// *********************************************************************

void read_atmos(p4est_iter_volume_info_t* info, void* user_data)
  // callback function following p4est_iter_volume_t prototype
{

  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;
  Context* context = static_cast<Context *>(info->p4est->user_pointer);
  patch::Patch* p = static_cast<patch::Patch*>(q->p.user_data);
  hdf::Hdf* f = static_cast<hdf::Hdf*>(user_data);

  mem::Array<float,2> tmp2d(p->yr, p->xr);
  const std::vector<int> start = {p->xs, p->ys};
  const std::vector<int> block = {p->xr, p->yr};
    
  f->readParallel("mass_density", "float", tmp2d.getDataBlock(), start, block);
  for (int j = 0; j < p->yr; ++j){
    for (int i = 0; i < p->xr; ++i){
      p->density(j, i) = tmp2d(j, i);
    }
   }
 
}

// *********************************************************************

void face_iteration_test(p4est_iter_face_info_t * info, void *user_data){

  p4est_iter_face_side_t* side[2];
  sc_array_t* sides = &(info->sides);

  /*
   * has this face leaves on both sides or just one (which means a
   * boundary face)
   */
 
  // every face has one side
  side[0] = p4est_iter_fside_array_index_int(sides, 0);
  //internal faces have two
  if (sides->elem_count == 2){
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
  }

  
  for (int i=0; i < sides->elem_count; ++i){
    printf("%i %i %i \n",
	   side[i]->face,
	   side[i]->is_hanging,
	   side[i]->is.full.is_ghost);
  }
  printf("\n");
  
    
}

// *********************************************************************

int main(int argc, char **argv)
{

  /*
   * initialize MPI using the MPI manager sc which is needed for 
   */
  int mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  int myid;
  sc_MPI_Comm_rank(mpicomm, &myid);

  /* 
   * These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. 
   */
  sc_init (mpicomm, 1, 1, nullptr, SC_LP_ESSENTIAL);
  p4est_init(nullptr, SC_LP_DEFAULT);

  hdf::init();

  /*
   * read the input
   */
  const std::string config_file = "test.cfg";
  ftl::section_t configuration = ftl::readConfigFile(config_file);

  /*
   * create global context
   */
  Context context(configuration);

  /* 
   * create forest, put in extra  { ... } block to ensure destructor
   * is called before MPI is finalized
   */
  {

    amr::Forest forest(mpicomm, configuration, 0 , context);

    /*
     * refine the initial forest, and redistribute evenly over all
     * processes
     */
    int recursive = 1;
    forest.refine(recursive, context.initial_refinement_level,
		  initial_refine_fn,
		  nullptr);
    forest.partition();
   
    /*
     * allocate a patch at each leaf, and allocate space for user data
     * in each patch
     */
  
    forest.iterateVolume(amr::allocate_patches);
    forest.iterateVolume(initialize_new_patches);
    
    // read atmosphere
    /*
    std::string atmosfile;
    ftl::getKeyword(configuration, "[INPUT_ATMOS]", "f", atmosfile);
    hdf::Hdf f;
    f.open(atmosfile, "r", mpicomm);
    forest.iterateVolume(read_atmos, nullptr, &f);
    f.close();
    forest.iterateVolume(fill_boundaries);
    */

    

    /* refine */

    /*
    recursive = 0;
    int allowed_level = context.initial_refinement_level + 2;
    for (int i = 0; i < 2; ++i){
      forest.refine(recursive, allowed_level,
		    density_refine_condition,
		    refine_patch_function);
      forest.balance(P4EST_CONNECT_FACE, replace_patch);
    } 
    */

    /*coarsen*/

    /*
    recursive = 0;
    for (int i = 0; i < 3; ++i){
      forest.coarsen(recursive, density_coarsen_condition, 
		     coarsen_patch_function);
      forest.balance(P4EST_CONNECT_FACE, replace_patch);
    } 

    forest.writeVTKFile(); 
    */ 
        
    /* create the ghost quadrants */
    forest.create_ghost_layer(P4EST_CONNECT_FULL);

    /*
      note that the following bracket calls the destructor on the forest
      object and deallocates all patches
    */
  }
  
  // finalize hdf
  hdf::finalize();
  
  //finalize MPI
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;

}
// *********************************************************************
