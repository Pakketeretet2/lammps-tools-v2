import block_data_
import data_field_
import data_field

import numpy as np



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
    """ The actual block data for atoms: """
    def __init__(self,meta,ids,types, x, mol = None):
        self.init_from_arrays(meta,ids,types,x,mol)

    @classmethod
    def init_from_handle(cls, handle):
        """ Initialises from a block_data_handle. """
        dom = domain_data( np.array( [0,0,0] ), np.array( [0,0,0] ), 0 )
        meta = block_meta( handle.time_step(), handle.n_atoms(), dom )

        ids   = block_data_.special_field_int( handle, 0 )
        types = block_data_.special_field_int( handle, 1 )
        mol   = None
        x = block_data_.special_field_double( handle, 3 )
        y = block_data_.special_field_double( handle, 4 )
        z = block_data_.special_field_double( handle, 5 )

        N = handle.n_atoms()
        X = np.empty( [ N, 3], dtype = float )
        X[:,0] = x
        X[:,1] = y
        X[:,2] = z

        # TODO: Domain data
        # TODO: Neatly check for mol.

        return cls(meta, ids, types, X, mol)


    def init_from_arrays(self,meta,ids,types, x, mol = None):
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


class block_data_local:
    """ The block data for a dump local: """
    def __init__(self, handle):

        self.dom = domain_data( np.array( [0,0,0] ), np.array( [0,0,0] ), 0 )
        self.meta = block_meta( handle.time_step(), handle.n_atoms(), self.dom )
        self.data = []
        self.data_names = []
        self.data_types = []

        """ Init block data from handle. """
        for i in range(0, block_data_.n_data_fields(handle)):
            df = block_data_.data_by_index(handle,i)
            self.data_names.append( data_field_.get_name( df ) )
            self.data_types.append( data_field_.get_type( df ) )

            # Grab the data as the proper underlying type:
            type = data_field_.get_type( df )
            if type == data_field_.TYPES.DOUBLE:
                tmp = data_field_.as_float( block_data_.data_by_index( handle, i ) )
            elif type == data_field_.TYPES.INT:
                tmp = data_field_.as_int( block_data_.data_by_index( handle, i ) )
            else:
                    raise RuntimeError("Unknown data type encountered in block_data")

            if self.data_types[i] == data_field_.TYPES.DOUBLE:
                tmp_data = np.zeros( self.meta.N, dtype = float )
            else: # Assume int.
                tmp_data = np.zeros( self.meta.N, dtype = int )

            for j in range(0, self.meta.N):
                tmp_data[j] = tmp[j]
            self.data.append( tmp_data )
