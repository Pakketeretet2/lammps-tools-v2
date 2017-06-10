import dump_reader_
import data_field_
import block_data_
import block_data

import sys

class dump_reader:
    """
    Interface to the C++ dump reader. This is the preferred method of
    reading dump files from Python due to its enhanced performance and better
    flexibility.

    It wraps the C handle in an __init__/__del__ pair and add some features.

    Initialise with a dump name, a file format and a dump format.
    Recognized file formats are: "PLAIN", "GZIP" and "BIN"
    Recognized dump formats are: "LAMMPS"

    If the dump format is "LAMMPS", you can pass an optional arg to indicate
    the dump file comes from a dump local instead.

    Examples:
    d = dump_reader.dump_reader( "dump_file.dump.bin", file_format = "BIN",
                                 dump_format = "LAMMPS" )

    """

    def __init__(self, fname, file_format, dump_format, is_local = False):
        """ Initialises dump reader. """
        self.local  = is_local
        self.handle = None

        fformat = dump_reader_.FILE_FORMATS.UNSET
        dformat = dump_reader_.DUMP_FORMATS.UNSET


        if file_format == "PLAIN":
            fformat = dump_reader_.FILE_FORMATS.PLAIN
        elif file_format == "GZIP":
            fformat = dump_reader_.FILE_FORMATS.GZIP
        elif file_format == "BIN":
            fformat = dump_reader_.FILE_FORMATS.BIN

        if dump_format == "LAMMPS":
            dformat = dump_reader_.DUMP_FORMATS.LAMMPS

        if (fformat == dump_reader_.FILE_FORMATS.UNSET or
            dformat == dump_reader_.FILE_FORMATS.UNSET):
            raise RuntimeError("Dump or file format not recognised!")

        if is_local:
            if fformat != dump_reader_.FILE_FORMATS.PLAIN:
                raise RuntimeError("Dump style local only works for " +
                                   "plain text format!")

        if self.local:
            self.handle = dump_reader_.new_local( fname, fformat, dformat )
        else:
            self.handle = dump_reader_.new( fname, fformat, dformat )

        status = dump_reader_.status(self.handle)
        if status !=  dump_reader_.DUMP_READER_STATUS.IS_GOOD:
            print("Error opening dump file ", fname, " for file format ",
                  fformat, " and dump format ", dformat, "!")
            print("Got status",status)
            sys.exit(-1)

    def __del__(self):
        """ Deletes the dump reader handle. """
        if self.handle != None:
            dump_reader_.delete( self.handle )

    def __iter__(self):
        """ Makes an iterator: """
        return self

    def __next__(self):
        if not self.at_eof():
            b = self.next_block()
            if self.status() == dump_reader_.DUMP_READER_STATUS.IS_GOOD:
                return b
        raise StopIteration

    def status(self):
        """ Returns the status of the internal dump reader. """
        status = dump_reader_.status( self.handle )
        return status

    def at_eof(self):
        """ Returns True if the dump file is at EOF. """
        status = self.status()
        return status == dump_reader_.DUMP_READER_STATUS.AT_EOF

    def next_block(self):
        """ Returns the next block_data. """
        bh = block_data_.block_data_handle()
        status = dump_reader_.get_next_block(self.handle, bh)
        if status == 0:
            if self.local:
                b = block_data.block_data_local( bh )
            else:
                b = block_data.block_data_custom.init_from_handle( bh )
            return b
        else:
            return None

    def set_column_headers(self, headers):
        """ Sets the column headers for the dump file. """
        if self.local:
            self.set_column_headers_local(headers)
        else:
            self.set_column_headers_custom(headers)

    def set_column_headers_local(self, headers):
        """ Sets the column headers for a dump file with local info. """

    def set_column_headers_custom(self, headers):
        """ Sets the column headers for a dump file containing atoms. """
        for i, h in zip( range(0,len(headers)), headers):
            dump_reader_.set_column_header( self.handle, i, h )

        header_set = [ False ]*6

        # Pretty names for the special columns:
        special_col_names = [ "ID", "MOL", "TYPE", "X", "Y", "Z" ]
        special_col_enum  = [ block_data_.SPECIAL_COLS.ID,
                              block_data_.SPECIAL_COLS.MOL,
                              block_data_.SPECIAL_COLS.TYPE,
                              block_data_.SPECIAL_COLS.X,
                              block_data_.SPECIAL_COLS.Y,
                              block_data_.SPECIAL_COLS.Z ]
        for h in headers:
            # Assume some defaults:
            if   h == "id":   i = 0
            elif h == "mol":  i = 1
            elif h == "type": i = 2
            elif h == "x":    i = 3
            elif h == "y":    i = 4
            elif h == "z":    i = 5
            # Add more here:
            # Bail if not guessable:
            else: continue

            name = special_col_names[i]

            if( header_set[i]  ):
                print( "Warning: Overwriting", name, "header.",
                       file = sys.stderr )
            else:
                print( "Guessing header", h, "\tis for", name,
                       file = sys.stderr )
            self.set_special_column( h, special_col_enum[i] )
            header_set[i] = True



    def set_special_column(self, header, special_field):
        """ Sets the given column header as special field. """
        dump_reader_.set_special_column( self.handle, header, special_field )


    def set_column_type(self, header, type):
        """ Sets the type of given column name. """
        dump_reader_.set_column_type( self.handle, header, type )

    def get_column_type(self, header):
        """ Returns the type of given column name. """
        return dump_reader_.get_column_type( self.handle, header )

    def set_default_column_type(self,type):
        """ Sets the default type for all columns in dump. """
        if type == "int":
            proper_type = data_field_.TYPES.INT
        elif type == "float" or type == "double":
            proper_type = data_field_.TYPES.DOUBLE
        else:
            raise RuntimeError( "Type " + type + " not recognized!" );

        dump_reader_.set_default_column_type( self.handle, proper_type )
