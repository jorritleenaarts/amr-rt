
// standard libraries
#include <iostream>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <cmath>

// p4est libraries


// local includes
#include "amr.hpp"
#include "patch.hpp"
#include "ftools.hpp"
#include "context.hpp"

// *********************************************************************

int initial_refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{

  Context* context = static_cast<Context *>(p4est->user_pointer);
  // cast int8_t as int so that cout knows to interpret it as an
  // int and not a char.
  //std::cout << int(quadrant->level) << "/" << context->initial_refinement_level
  //	    << std::endl;

  if (quadrant->level < context->initial_refinement_level) {
    return 1;
  } else {
    return 0;
  }

}

// *********************************************************************

void initialize_quadrant(p4est_t * p4est, p4est_topidx_t which_tree,
	      p4est_quadrant_t * quadrant)
{

  Context* context = static_cast<Context *>(p4est->user_pointer);
  patch::Patch* lpatch = static_cast<patch::Patch *>(quadrant->p.user_data);

}

// *********************************************************************

void iterate_fn(p4est_iter_volume_info_t* info, void* user_data)
{

  p4est_quadrant_t* q = info->quad;
  patch::Patch* p = static_cast<patch::Patch*>(q->p.user_data);
  p4est_topidx_t which_tree = info->treeid;
  Context* context = static_cast<Context *>(info->p4est->user_pointer);
  

  //compute pixel coordinatse of lower left corner in input atmos
  int xmin = (q->x * context->nx_in_atmos)/P4EST_ROOT_LEN;
  int ymin = (q->y * context->ny_in_atmos)/P4EST_ROOT_LEN;

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

   //std::cout << p->xrb << " " << p->yrb << std::endl;
   //p->density.resize(-ng, p->xr + 1,-ng, p->yr + 1 );
   //p->dummy.resize(p->xrb*p->yrb);

   p->dummy2 = (int *) malloc(sizeof(int)*p->xrb*p->yrb);

   
     /*for (int j = 0; j < p->yr ; ++j){
      for (int i = 0; i < p->xr; ++i){
	p->density(j,i) =  10000 * (xmin + i) + (ymin + j);
     }
     } */
   // printf("%d  %08d \n",p->density(-1,-1), p->density(0,0));
 
  
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

  // compute initial refinement
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "nx", context.nx_in_atmos);
  ftl::getKeyword(configuration, "[INPUT_ATMOS]", "ny", context.ny_in_atmos);
  
  ftl::getKeyword(configuration, "[PATCH]", "nx_per_patch",
		  context.nx_per_patch);
  ftl::getKeyword(configuration, "[PATCH]", "ny_per_patch",
		  context.ny_per_patch);
  
  int x_ref_lvl = log2(context.nx_in_atmos / context.nx_per_patch);
  // for now a crappy implementation 
  context.initial_refinement_level = 7; //x_ref_lvl; 

  /* 
   * create forest, put in extra  { ... } block to ensure destructor
   * is called before MPI is finalized
   */
  {

    std::cout <<  "patch size:" << sizeof(patch::Patch) << std::endl;
    amr::Amr AMR(mpicomm, configuration, sizeof(patch::Patch), context);
    
    // refine the initial forest
    int recursive = 1;
    AMR.refineForest(recursive, initial_refine_fn, nullptr);

    /* 
     * refine may lead to uneven number of leaves per MPI process,
     * so redistribute
     */
    AMR.partitionForest(0, nullptr);

    /* 
     * iterate over all leaves and read initial atmosphere
     */
    AMR.iterateForest(nullptr, &myid, iterate_fn , nullptr, nullptr);

  }

  // *********************************************************************
  // finalize MPI
  // *********************************************************************

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;

}
// *********************************************************************
