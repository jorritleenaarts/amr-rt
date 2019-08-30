
// standard libraries
#include <iostream>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <cmath>

// p4est libraries
#include "p4est.h"
#include "p4est_iterate.h"

// local includes
#include "amr.hpp"
#include "patch.hpp"
#include "ftools.hpp"
#include "context.hpp"

// *********************************************************************

int initial_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant){

  Context* context = static_cast<Context *>(p4est->user_pointer);
  // cast int8_t as int so that cout knows to interpret it as an int and not a char.
  std::cout << int(quadrant->level) << "/" << context->initial_refinement_level << std::endl;

  if (quadrant->level < context->initial_refinement_level) {
    return 1;
  } else {
    return 0;
  }

}

// *********************************************************************

void initialize_quadrant(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant){

  Context* context = static_cast<Context *>(p4est->user_pointer);
  patch::Patch* lpatch = static_cast<patch::Patch *>(quadrant->p.user_data);



}

// *********************************************************************

void iterate_fn(p4est_iter_volume_info_t* info, void* user_data)
{

  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;

  int* myid_ptr = (int *) user_data;
  printf("myid = %d; tree = %d; level = %d \n",*myid_ptr, which_tree ,q->level);
  
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
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init(NULL, SC_LP_DEFAULT);


  // *********************************************************************
  // read the input
  // *********************************************************************

  const std::string config_file = "test.cfg";
  ftl::section_t configuration = ftl::readConfigFile(config_file);


  /*
   * create global context
   */
  Context context;

  // compute initial refinement
  int nx;
  int ny;
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "nx", nx);
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "ny", ny);
  
  ftl::getKeyword(configuration, "[PATCH]", "nx_per_patch", context.nx_per_patch);
  ftl::getKeyword(configuration, "[PATCH]", "ny_per_patch", context.ny_per_patch);
  
  int x_ref_lvl = log2(nx / context.nx_per_patch);
  std::cout <<nx << " " << context.nx_per_patch << " " << x_ref_lvl << std::endl;
  context.initial_refinement_level = 1;// x_ref_lvl; // for now a crappy implementation 

  /* 
   * create forest, put in extra  { ... } block to ensure destructor
   * is called before MPI is finalized
   */
  {

    std::cout <<  "patch size:" << sizeof(patch::Patch) << std::endl;
    amr::Amr AMR(mpicomm, configuration, sizeof(patch::Patch), context);
    
    // refine the initial forest
    int recursive = 1;
    AMR.refineForest(recursive, initial_refine_fn, NULL);

  }

  /*
  // refine lead to uneven number of leaves per MPI process, so redistribute
  p4est_partition (forest, 0, NULL);

  // iterate over all leaves
  p4est_iterate(forest, NULL, &myid, iterate_fn , NULL, NULL);
  
  //  std::cout << "Hello World from rank " << myid << std::endl;

  // Destroy the p4est and the connectivity structure. 
  p4est_destroy (forest);
  p4est_connectivity_destroy (connectivity);

  */


  // *********************************************************************
  // finalize MPI
  // *********************************************************************

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;

}
// *********************************************************************
