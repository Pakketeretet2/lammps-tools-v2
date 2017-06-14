#ifndef LT_TRANSFORMATIONS_H
#define LT_TRANSFORMATIONS_H

#include "lt_block_data.h"

extern "C" {

void lt_transformations_rotate_all( lt_block_data_handle *bdh,
                                    const double *axis,
                                    const double *origin, double angle );

void lt_transformations_shift_all( lt_block_data_handle *bdh,
                                   const double *delta );

void lt_transformations_center_box_on( lt_block_data_handle *bdh,
                                       const double *o );

void lt_transformations_shift( lt_block_data_handle *bdh, const double *delta,
                               size_t id_size, const int *ids );

void lt_transformations_rotate( lt_block_data_handle *bdh, const double *axis,
                                const double *origin, double angle,
                                size_t id_size, const int *ids );

void lt_transformations_unfold_mols( lt_block_data_handle *bdh );


} // extern "C"

// Because pybind11 supports vectors but not arrays.
void lt_transformations_rotate_all_pb( lt_block_data_handle *bdh,
                                       const std::vector<double> &axis,
                                       const std::vector<double> &origin,
                                       double angle )
{
	lt_transformations_rotate_all( bdh, axis.data(), origin.data(), angle );
}

void lt_transformations_shift_all_pb( lt_block_data_handle *bdh,
                                      const std::vector<double> &delta )
{
	lt_transformations_shift_all( bdh, delta.data() );
}

void lt_transformations_center_box_on_pb( lt_block_data_handle *bdh,
                                          const std::vector<double> &o )
{
	lt_transformations_center_box_on( bdh, o.data() );
}

void lt_transformations_shift_pb( lt_block_data_handle *bdh,
                                  const std::vector<double> &delta,
                                  const std::vector<int> &ids )
{
	lt_transformations_shift( bdh, delta.data(), ids.size(), ids.data() );
}

void lt_transformations_rotate_pb( lt_block_data_handle *bdh,
                                   const std::vector<double> &axis,
                                   const std::vector<double> &origin,
                                   double angle,
                                   const std::vector<int> &ids )
{
	lt_transformations_rotate( bdh, axis.data(), origin.data(), angle,
	                           ids.size(), ids.data() );
}



#endif // LT_TRANSFORMATIONS_H
