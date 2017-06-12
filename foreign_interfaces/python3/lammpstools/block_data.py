import block_data_
import data_field_
import data_field

import numpy as np

class xyz_array_acessor:
    """ Provides an abstraction to indexing the xyz array in block data. """
    def __init__(self, x_arr, y_arr, z_arr):
        self.x = x_arr
        self.y = y_arr
        self.z = z_arr

    def __getitem__(self, atom_index):
        return (self.x[atom_index], self.y[atom_index], self.z[atom_index])

    
    # def __getitem__(self, atom_index, coord_index):
    #     if coord_index == 0:
    #         return self.x[atom_index]
    #     elif coord_index == 1:
    #         return self.y[atom_index]
    #     elif coord_index == 2:
    #         return self.z[atom_index]
    #     else:
    #         raise RuntimeError("coordinate index",coord_index,"out of bounds!")

    def __len__(self):
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


    def data_by_name(self, name):
        """ Grab data field by name. """
        # Lazy create name map:
        if self.stored_name_map is False:
            self.store_name_mapping()

        raw_data = self.data[ self.name_to_col[ name ] ]
        return raw_data

        # if df is None:
        #     raise RuntimeError("Failed to find data field of name", name)
        # # Resolve type here:
        # if data_field_.get_type( df ) == data_field_.TYPES.INT:
        #     return data_field_.as_int( df )
        # elif data_field_.get_type( df ) == data_field_.TYPES.DOUBLE:
        #     return data_field_.as_float( df )
        # else:
        #     raise RuntimeError("Failed to resolve data type for data field",
        #                        data_field_.get_name(df))


    def get_ref(self):
        """ Returns a ref to the pointer contained. This eases some foreign
            function calling but shouldn't be used inside Python too much! """
        return self.handle.get_const_ref()



class block_data_custom(block_data):
    """ The actual block data for atoms: """
    def __init__(self, meta, ids, types, x, mol = None, handle = None, no_copy = False):
        """ Initialiser. """
        super(block_data_custom,self).__init__(handle)
        self.no_copy = no_copy
        self.init_from_arrays(meta,ids,types,x,mol)

    @classmethod
    def init_from_handle(cls, handle, no_copy = False):
        """ Initialises from a block_data_handle. """
        dom = domain_data( np.array( [0,0,0] ), np.array( [0,0,0] ), 0 )
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
        return cls(meta, ids, types, X, mol, handle, no_copy)


    def init_from_arrays(self,meta,ids,types, x, mol = None, no_copy = False):
        """ Initialises from np arrays """
        # Check lengths:
        if not all(meta.N == length for length in
                   [ len(ids), len(types), len(x) ] ):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")

        if (not mol is None) and (len(mol) != meta.N):
            raise RuntimeError("Length mismatch in arrays passed to block_data!")

        self.meta  = meta
        self.ids   = ids
        self.types = types
        self.x     = x
        if mol is None:
            self.mol = None
        else:
            self.mol = mol
            self.meta.atom_style = "molecular"


class block_data_local(block_data):
    """ The block data for a dump local: """

    def __init__(self, handle):
        super(block_data_local,self).__init__(handle)

        self.dom = domain_data( np.array( [0,0,0] ), np.array( [0,0,0] ), 0 )
        self.meta = block_meta( handle.time_step(), handle.n_atoms(), self.dom )
        self.store_name_mapping()
