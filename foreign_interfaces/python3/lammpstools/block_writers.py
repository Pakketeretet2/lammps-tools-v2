import block_writers_

def to_lammps_data( fname, wmode, b ):
    """ Writes given block_data to lammps data file. """
    return block_writers_.to_lammps_data( fname, wmode, b.handle )

def to_lammps_dump( fname, wmode, b ):
    """ Writes given block_data to lammps dump file. """
    return block_writers_.to_lammps_dump( fname, wmode, b.handle )
