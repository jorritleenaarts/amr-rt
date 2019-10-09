#ifndef AMR_HPP
#define AMR_HPP

#include "p4est.h"
#include "p4est_iterate.h"
#include "p4est_extended.h"
#include <p4est_vtk.h>

#include "patch.hpp"
#include "ftools.hpp"
#include "context.hpp"

namespace amr{

  /*
   * allocates the patch connected to the quadrant
   */
  void allocate_patch(p4est_quadrant_t* q);

  /*
   * deallocates the patch connected to the quadrant, set user_data =
   * nullptr
   */
  void deallocate_patch(p4est_quadrant_t* q);

  /*
   * helper functions that conform to the p4est_iter_volume_t
   * function prototype
   */

  /*
   * deallocate the patch objects in each leaf of the forest
   */
  void deallocate_patches(p4est_iter_volume_info_t * info,
			  void *user_data);

  /*
   * allocate a patch object in each leaf of the forest
   */
  void allocate_patches(p4est_iter_volume_info_t * info,
			void *user_data);

  /*
   * helper functions that conform to the p4est_init_t function
   * prototype
   */

  /*
   * this function should be called whenever a new quadrant is
   * created using the p4est_[new,refine,coarsen] functions.
   */
  void nullify_user_data(p4est_t * p4est,
			 p4est_topidx_t which_tree,
			 p4est_quadrant_t * quadrant);


  class Forest{
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

    Forest(sc_MPI_Comm mpicomm, ftl::section_t configuration, 
	   size_t user_data_size, Context& context);

    ~Forest();

    void partition(int allow_for_coarsening = 0,
		   p4est_weight_t weight_fn = nullptr);

    void iterate(p4est_ghost_t * ghost_layer,
		 void *user_data,
		 p4est_iter_volume_t iter_volume,
		 p4est_iter_face_t iter_face,
		 p4est_iter_corner_t iter_corner);

    void iterateVolume(p4est_iter_volume_t iter_volume,
		       p4est_ghost_t* ghost_layer = nullptr,
		       void* user_data = nullptr);

    void refine(int refine_recursive, int maxlevel,
		p4est_refine_t refine_fn,
		p4est_replace_t replace_fn);
    
    void coarsen(int coarsen_recursive,
		 p4est_coarsen_t coarsen_fn,
		 p4est_replace_t replace_fn);

    void balance(p4est_connect_type_t btype,
		 p4est_replace_t replace_fn);

    void writeVTKFile();

  };

}

#endif // AMR_HPP
