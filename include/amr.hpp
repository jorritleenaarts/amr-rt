#ifndef AMR_HPP
#define AMR_HPP

#include "p4est.h"
#include "p4est_iterate.h"
#include "p4est_extended.h"

#include "ftools.hpp"
#include "context.hpp"

namespace amr{

  class Amr{
/*
 * wrapper around the p4est C library
 */

  private:

    // hard coded parameters
    const int periodic_x = 0; // not periodic in x
    const int periodic_y = 0; // not periodic in y

    // input parameters
    int n_trees_x = 1; //number of trees in x
    int n_trees_y = 1; // number of trees in y

    // the p4est objects that
    p4est_connectivity_t* connectivity;
    p4est_t* forest;

  public:

    Amr(sc_MPI_Comm mpicomm, ftl::section_t configuration, 
	size_t user_data_size, Context& context);

    ~Amr();

    void refineForest(int refine_recursive, p4est_refine_t refine_fn, p4est_init_t init_fn);

    void partitionForest(int allow_for_coarsening, p4est_weight_t weight_fn);

    void iterateForest(p4est_ghost_t * ghost_layer,
		       void *user_data,
		       p4est_iter_volume_t iter_volume,
		       p4est_iter_face_t iter_face,
		       p4est_iter_corner_t iter_corner);

    void refine(int refine_recursive, int maxlevel,
		   p4est_refine_t refine_fn,
		   p4est_init_t init_fn,
		p4est_replace_t replace_fn);

    void coarsen(int coarsen_recursive,
		 int callback_orphans,
		 p4est_coarsen_t coarsen_fn,
		 p4est_init_t init_fn,
		 p4est_replace_t replace_fn);

  };

}

#endif // AMR_HPP
