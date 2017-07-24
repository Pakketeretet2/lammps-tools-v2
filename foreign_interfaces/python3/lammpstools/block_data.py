import block_data_
import data_field_
import data_field

import numpy as np
import sys


SPECIAL_COLS_ID   = block_data_.SPECIAL_COLS.ID
SPECIAL_COLS_TYPE = block_data_.SPECIAL_COLS.TYPE
SPECIAL_COLS_MOL  = block_data_.SPECIAL_COLS.MOL
SPECIAL_COLS_X    = block_data_.SPECIAL_COLS.X
SPECIAL_COLS_Y    = block_data_.SPECIAL_COLS.Y
SPECIAL_COLS_Z    = block_data_.SPECIAL_COLS.Z
SPECIAL_COLS_VX   = block_data_.SPECIAL_COLS.VX
SPECIAL_COLS_VY   = block_data_.SPECIAL_COLS.VY
SPECIAL_COLS_VZ   = block_data_.SPECIAL_COLS.VZ
SPECIAL_COLS_IX   = block_data_.SPECIAL_COLS.IX
SPECIAL_COLS_IY   = block_data_.SPECIAL_COLS.IY
SPECIAL_COLS_IZ   = block_data_.SPECIAL_COLS.IZ


class xyz_array_acessor:
    """ Provides an abstraction to indexing the xyz array in block data. """
    def __init__(self, x_arr, y_arr, z_arr):
        """ Initialises references to the arrays. """
        if len(x_arr) != len(z_arr) or len(x_arr) != len(y_arr):
            raise RuntimeError("Arrays not of equal lengths!")

        self.x = x_arr
        self.y = y_arr
        self.z = z_arr

    def __getitem__(self, atom_index):
        """ Returns a !copy! of atom coordinates. """
        return np.array( [self.x[atom_index], self.y[atom_index], self.z[atom_index] ])

    def __setitem__(self, atom_index, new_values):
        """ Returns a !copy! of atom coordinates. """
        self.x[atom_index] = new_values[0]
        self.y[atom_index] = new_values[1]
        self.z[atom_index] = new_values[2]

    def __len__(self):
        """ Returns the length of the arrays. """
        return len(self.x)


class domain_data:
    """ Some basic info about the simulation domain """
    def __init__( self, xlo, xhi, periodic ):
        self.xlo = xlo
        self.xhi = xhi
        self.periodic = periodic

    def distance(self, xi, xj):
        """ Calculates the Euclidian distance respecting the PBC """
        Lx = self.xhi[0] - self.xlo[0]
        Ly = self.xhi[1] - self.xlo[1]
        Lz = self.xhi[2] - self.xlo[2]
        r  = xi - xj
        if self.periodic & 1:
            if r[0] >  0.5*Lx:
                r[0] -= Lx
            elif r[0] <= -0.5*Lx:
                r[0] += Lx

        if self.periodic & 2:
            if r[1] >  0.5*Ly:
                r[1] -= Ly
            elif r[1] <= -0.5*Ly:
                r[1] += Ly

        if self.periodic & 4:
            if r[2] >  0.5*Lz:
                r[2] -= Lz
            elif r[2] <= -0.5*Lz:
                r[2] += Lz

        dist = np.linalg.norm(r)

        return dist, r


class block_meta:
    """ Some block metadata """
    def __init__(self,t,N,domain):
        self.t = t
        self.N = N
        self.domain = domain
        self.atom_style = "atomic"


class block_data:
    """ A class representing a block of data. """
    def __init__(self, handle = None):
        """ Initialiser. Empty."""
        self.stored_name_map = False
        if handle is None:
            # self.handle = new_block_data_handle()
            self.handle = block_data_.block_data_handle()
        else:
            self.handle = handle

        self.name_to_col = None
        self.data_types  = None
        self.data = None

    def store_name_mapping(self):
        # Grab all names for mapping:
        N = block_data_.n_data_fields(self.handle)
        self.name_to_col = {}
        self.data = []
        self.data_types = np.zeros( N, dtype = int )

        for i in range(0, N):
            df = block_data_.data_by_index(self.handle,i)
            ttype = data_field_.get_type( df )
            self.data_types[i] = ttype

            if ttype == data_field_.TYPES.DOUBLE:
                raw_data = data_field_.as_float(df)
            elif ttype == data_field_.TYPES.INT:
                raw_data = data_field_.as_int(df)
            else:
                raise RuntimeError("Unkown data type encountered!")
            self.data.append(raw_data)
            name = data_field_.get_name(df)
            self.name_to_col[ name ] = i
        self.stored_name_map = True

    def n_data_fields(self):
        """ Returns the number of data fields. """
        return block_data_.n_data_fields(self.handle)

    def name_of_data(self, index):
        """ Returns the name of given data index. """
        if index < 0 or index >= self.n_data_fields():
            raise RuntimeError("Index out of bounds!")

        df = block_data_.data_by_index(self.handle, index)
        return data_field_.get_name(df)

    def data_by_name(self, name):
        """ Grab raw data by name. """
        # Lazy create name map:
        if self.stored_name_map is False:
            self.store_name_mapping()

        raw_data = self.data[ self.name_to_col[ name ] ]
        return raw_data

    def overwrite_data_by_name(self, name, data):
        """ Overwrites given data with new data. """
        if not self.has_data(name):
            raise RuntimeError("Data named", name, "not in block_data!")

        col_idx = self.name_to_col[name]
        data_type = self.data_types[col_idx]



    def has_data(self, name):
        """ Checks if block data has named data. """
        try:
            df = self.data_by_name( name )
            return True
        except KeyError:
            # Apparently named data is not  there.
            return False

    def get_ref_(self):
        """ Returns a ref to the pointer contained. This eases some foreign
            function calling but shouldn't be used inside Python too much! """
        return self.handle.get_const_ref()

    def get_ptr_(self):
        """ Returns a ref to the pointer contained. This eases some foreign
            function calling but shouldn't be used inside Python too much! """
        return self.handle.get_ptr()

    def remove_data_field(self, name):
        """ Removes given data field from block data. """
        block_data_.remove_field( self.handle, name )

    def replace_data_field(self, name, new_data_field):
        """ Replaces given data field from block data. """
        block_data_.swap_fields( self.handle, name, new_data_field.handle )


    def add_data_field(self, d, special_field_type = None):
        """ Adds given data_field to this block. """
        if special_field_type is None:
            block_data_.add_data_field( self.handle, d.handle )
        else:
            block_data_.add_special_field( self.handle, d.handle,
                                           special_field_type )

        # Check if the data field is properly added:
        test_data_field = block_data_.data_by_name( self.handle, d.name() )
        if data_field.data_field(test_data_field).name() != d.name():
            raise RuntimeError("Data field mismatch!")

        print("I just added a field named", d.name())

        # Reset name_mapping to force this to be in sync:
        self.store_name_mapping()

    def set_domain(self, xlo = None, xhi = None, periodic_bits = None):
        """ Sets domain info """
        if not xlo is None:
            self.meta.domain.xlo[0] = xlo[0]
            self.meta.domain.xlo[1] = xlo[1]
            self.meta.domain.xlo[2] = xlo[2]

        if not xhi is None:
            self.meta.domain.xhi[0] = xhi[0]
            self.meta.domain.xhi[1] = xhi[1]
            self.meta.domain.xhi[2] = xhi[2]

        if not periodic_bits is None:
            self.meta.domain.periodic = periodic_bits

        block_data_.set_domain(self.handle,
                               self.meta.domain.xlo[0], self.meta.domain.xhi[0],
                               self.meta.domain.xlo[1], self.meta.domain.xhi[1],
                               self.meta.domain.xlo[2], self.meta.domain.xhi[2],
                               self.meta.domain.periodic)



class block_data_custom(block_data):
    """ The actual block data for atoms: """
    def __init__(self, meta, ids, types, x, mol = None, handle = None):
        """ Initialiser. """
        super(block_data_custom,self).__init__(handle)
        if handle is None:
            self.init_from_arrays(meta,ids,types,x,mol, True)
        else:
            self.init_from_arrays(meta,ids,types,x,mol, False)

    @classmethod
    def init_from_handle(cls, handle, no_copy = False):
        """ Initialises from a block_data_handle. """
        dom = domain_data( np.array( [0.0,0.0,0.0] ),
                           np.array( [0.0,0.0,0.0] ), 0 )
        meta = block_meta( handle.time_step(), handle.n_atoms(), dom )

        N = handle.n_atoms()
        ids   = block_data_.special_field_int( handle, 0 )
        types = block_data_.special_field_int( handle, 1 )
        mol   = None
        x = block_data_.special_field_double( handle, 3 )
        y = block_data_.special_field_double( handle, 4 )
        z = block_data_.special_field_double( handle, 5 )

        if no_copy:
            X = xyz_array_acessor( x, y, z )
        else:
            X = np.empty( [ N, 3], dtype = float )
            X[:,0] = x
            X[:,1] = y
            X[:,2] = z

        if( block_data_.has_special_field( handle, 2 ) ):
            mol = block_data_.special_field_int( handle, 2 )
        return cls(meta, ids, types, X, mol, handle)


    def init_from_arrays(self,meta,ids,types, x, mol = None,
                         add_fields_to_handle = False):
        """ Initialises from np arrays """
        # Check lengths:
        if not all(meta.N == length for length in
                   [ len(ids), len(types), len(x) ] ):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")

        if (not mol is None) and (len(mol) != meta.N):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")

        self.meta  = meta

        if add_fields_to_handle:
            self.handle.set_n_atoms( meta.N )

            atom_style_int = 0
            if self.meta.atom_style == "molecular": atom_style_int = 1

            block_data_.set_meta( self.handle, meta.t, meta.N,
                                  meta.domain.xlo[0], meta.domain.xhi[0],
                                  meta.domain.xlo[1], meta.domain.xhi[1],
                                  meta.domain.xlo[2], meta.domain.xhi[2],
                                  meta.domain.periodic, atom_style_int )

            max_type = np.max( types )
            self.handle.set_n_types( max_type )

            int_type = data_field.DATA_TYPE_INT
            float_type = data_field.DATA_TYPE_DOUBLE

            d_id   = data_field.new_data_field(   "id", int_type, meta.N )
            d_type = data_field.new_data_field( "type", int_type, meta.N )
            d_x    = data_field.new_data_field(    "x", float_type, meta.N )
            d_y    = data_field.new_data_field(    "y", float_type, meta.N )
            d_z    = data_field.new_data_field(    "z", float_type, meta.N )

            for i in range(0,meta.N):
                d_id[i] = ids[i]
                d_type[i] = types[i]
                d_x[i] = x[i,0]
                d_y[i] = x[i,1]
                d_z[i] = x[i,2]

            self.add_data_field(   d_id, SPECIAL_COLS_ID )
            self.add_data_field( d_type, SPECIAL_COLS_TYPE )
            self.add_data_field(    d_x, SPECIAL_COLS_X )
            self.add_data_field(    d_y, SPECIAL_COLS_Y )
            self.add_data_field(    d_z, SPECIAL_COLS_Z )

            if mol is not None:
                d_mol = data_field.new_data_field( "mol", int_type, meta.N )
                data_field.copy_to_data_field( mol, d_mol )
                self.add_data_field( d_mol, SPECIAL_COLS_MOL )
                self.meta.atom_style = "molecular"

        # This needs to happen no matter what:
        self.ids   = block_data_.special_field_int( self.handle, 0 )
        self.types = block_data_.special_field_int( self.handle, 1 )

        xf = block_data_.special_field_double( self.handle, 3 )
        yf = block_data_.special_field_double( self.handle, 4 )
        zf = block_data_.special_field_double( self.handle, 5 )

        self.x = xyz_array_acessor( xf, yf, zf )

        if mol is None:
            self.mol = None
        else:
            self.mol = block_data_.special_field_int( self.handle, 2 )
            self.handle.set_atom_style( 1 )

class block_data_local(block_data):
    """ The block data for a dump local: """

    def __init__(self, handle):
        super(block_data_local,self).__init__(handle)

        self.dom = domain_data( np.array( [0,0,0] ), np.array( [0,0,0] ), 0 )
        self.meta = block_meta( handle.time_step(), handle.n_atoms(), self.dom )
        self.store_name_mapping()
