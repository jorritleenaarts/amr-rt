
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

int initial_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{

  Context* context = static_cast<Context *>(p4est->user_pointer);
 
  if (quadrant->level < context->initial_refinement_level) {
    return 1;
  } else {
    return 0;
  }

}

// *********************************************************************

int density_refine_condition(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{

  Context* context = static_cast<Context *>(p4est->user_pointer);
  patch::Patch* p = static_cast<patch::Patch*>(quadrant->p.user_data);

  const float density_gradient_threshold = 0.2;


  mem::Array<float,3> grad = math::grid_gradient(p->density);
  std::cout << "refine = 1" << std::endl;
   return 1;
   
  // compute gradient of density
 
  mem::Array<float,2> absgrad = math::abs_array(grad);

 
  
  for (int j = p->ys; j < p->ye; ++j ){
    for (int i = p->xs; i < p->xe; ++i ){
      if (absgrad(j,i) /= p->density(j,i) > density_gradient_threshold){
	std::cout << "refine = 1" << std::endl;
	return 1;
      }
    }
  }
  std::cout << "refine = 0" << std::endl;
  return 0;

}

// *********************************************************************

void refine_patch_function(p4est_t * p4est,
			    p4est_topidx_t which_tree,
			    int num_outgoing,
			    p4est_quadrant_t * outgoing[],
			    int num_incoming,
			    p4est_quadrant_t * incoming[]){

			   }

// *********************************************************************

void nullify_user_data_ptr(p4est_t* p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t* quadrant)
// callback function following p4est_init_t prototype
{
  // CHRASHES WHEN RUNNING WITH MORE THAN ONE MPI PROCESS
  // quadrant->p.user_data = nullptr;
}

// *********************************************************************

void allocate_new_patches(p4est_iter_volume_info_t* info, void* user_data)
  // callback function following p4est_iter_volume_t prototype
{

  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;

  
  if (q->p.user_data != nullptr){
    printf("q->p.user_data not null at tree %i with refinement level %i\n",
	   info->treeid, q->level);
    printf("This might be a memory leak, or the leaves are not properly initialized.\n");
    //  abort();
  }
  
  
  // allocate patch and copy into user data pointer.
  std::cout << "making new patch on" << std::endl;
  patch::Patch*  p = new patch::Patch();
  q->p.user_data = p;
   
}

// *********************************************************************

void initialize_new_patches(p4est_iter_volume_info_t* info, void* user_data)
  // callback function following p4est_iter_volume_t prototype
{

  p4est_quadrant_t* q = info->quad;
  p4est_topidx_t which_tree = info->treeid;
  Context* context = static_cast<Context *>(info->p4est->user_pointer);
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
   p->xr = xmax - xmin + 1;
   p->yr = ymax - ymin + 1;
   
   p->xsb = xmin - ng;
   p->xeb = xmax + ng;
   p->ysb = ymin - ng;
   p->yeb = ymax + ng;
   p->xrb = p->xr + 2 * ng;
   p->yrb = p->yr + 2 * ng;

   p->density.resize(-ng,p->xr,-ng,p->yr);
   std::cout << "here: " <<  p->density.size() << std::endl;   
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
  //printf("density(0,0)=%e\n",p->density(0,0));
  
}

// *********************************************************************

void face_iteration_test(p4est_iter_face_info_t * info, void *user_data){

  p4est_iter_face_side_t* side[2];
  sc_array_t* sides = &(info->sides);

  // has this face leaves on both sides or just one (which means a boundary face)
  //std::cout << sides->elem_count << " "<< std::endl;

  // every face has one side
  side[0] = p4est_iter_fside_array_index_int(sides, 0);
  //internal faces have two
  if (sides->elem_count == 2){
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
  }

  /*
  for (int i=0; i < sides->elem_count; ++i){
    printf("%i %i %i \n",
	   side[i]->face,
	   side[i]->is_hanging,
	   side[i]->is.full.is_ghost);
  }
  printf("\n");
  */
    
}

// *********************************************************************

int main(int argc, char **argv)
{

  // *********************************************************************
  // initialize MPI using the MPI manager sc which is needed for p4est
  // *********************************************************************
  
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

  // *********************************************************************
  // read the input
  // *********************************************************************

  const std::string config_file = "test.cfg";
  ftl::section_t configuration = ftl::readConfigFile(config_file);


  /*
   * create global context
   */
  Context context;

 
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "nx", context.nx_in_atmos);
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "ny", context.ny_in_atmos);
  
  ftl::getKeyword(configuration, "[PATCH]", "nx_per_patch",
		  context.nx_per_patch);
  ftl::getKeyword(configuration, "[PATCH]", "ny_per_patch",
		  context.ny_per_patch);

  // compute initial refinement
  // for now a crappy implementation 
  int x_ref_lvl = log2(context.nx_in_atmos / context.nx_per_patch);
  context.initial_refinement_level = 1; // x_ref_lvl; 

  /* 
   * create forest, put in extra  { ... } block to ensure destructor
   * is called before MPI is finalized
   */
  {

    amr::Amr AMR(mpicomm, configuration, sizeof(patch::Patch), context);

    // initialize hdf
    hdf::init();

    
    // refine the initial forest
    //
    int recursive = 1;
    AMR.refineForest(recursive,
    		     initial_refine_fn,
		     nullify_user_data_ptr);

    /* 
     * refine may lead to uneven number of leaves per MPI process,
     * so redistribute
     */
    AMR.partitionForest(0, nullptr);

    /* 
     * iterate over all leaves and read initial atmosphere
     */
    AMR.iterateForest(nullptr, &myid, allocate_new_patches, nullptr, nullptr);
    AMR.iterateForest(nullptr, &myid, initialize_new_patches, nullptr, nullptr);
    
    std::string atmosfile;
    ftl::getKeyword(configuration, "[INPUT_ATMOS]", "f", atmosfile);

     hdf::Hdf f;
    f.open(atmosfile, "r", mpicomm);
    AMR.iterateForest(nullptr, &f, read_atmos, nullptr, nullptr);
    f.close();


    /* adapt */
    recursive = 1;
    int allowed_level = context.initial_refinement_level + 2;
    //    AMR.refine(recursive, allowed_level, density_refine_condition,
    //	       nullptr, refine_patch_function);
    AMR.refine(recursive, allowed_level, density_refine_condition,
	       nullptr, nullptr);

    
    /*
    recursive = 1
    AMR::coarsen(recursive, callbackorphans,
      step3_coarsen_err_estimate, NULL,
      step3_replace_quads);
	p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL,
	step3_replace_quads);
      */    
	

    /* create the ghost quadrants */
    //    p4est_ghost_t* ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  
    // AMR.iterateForest(nullptr, &myid, nullptr, face_iteration_test, nullptr);
    
  }

  // finalize hdf
  hdf::finalize();
  
  // *********************************************************************
  // finalize MPI
  // *********************************************************************

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;

}
// *********************************************************************
