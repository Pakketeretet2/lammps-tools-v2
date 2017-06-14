import transformations_

def rotate_all( b, axis, origin, angle ):
    transformations_.rotate_all( b.handle, axis, origin, angle )

def shift_all( b, delta ):
    transformations_.shift_all( b.handle, delta )

def center_box_on( b, origin ):
    transformations_.center_box_on( b.handle, origin )

def rotate( b, axis, origin, angle, ids ):
    transformations_.rotate( b.handle, axis, origin, angle, ids )

def shift( b, delta, ids ):
    transformations_.shift( b.handle, delta, ids )

def unfold_mols( b ):
    transformations_.unfold_mols( b.handle )
