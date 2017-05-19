import dump_reader_
import block_data_
import block_data

class dump_reader:
    """
    Interface to the C++ dump reader. This is the preferred method of
    reading dump files from Python due to its enhanced performance and better
    flexibility.

    It wraps the C handle in an __init__/__del__ pair and add some features.
    """

    def __init__(self, fname, fformat, dformat):
        """ Initialises dump reader. """
        print("Called dump_reader.__init__")
        self.handle = dump_reader_.new_dump_reader( fname, fformat, dformat )

    def __del__(self):
        """ Deletes the dump reader handle. """
        dump_reader_.delete_dump_reader( self.handle )

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
        status = dump_reader_.dump_reader_status( self.handle )
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
            b = block_data.block_data.init_from_handle( bh )
            return b
        else:
            return None

    def set_column_headers(self, headers):
        """ Sets the column headers for the dump file. """
        for i, h in zip( range(0,len(headers)), headers):
            dump_reader_.set_column_header( self.handle, i, h )


    def set_special_column(self, header, special_field):
        """ Sets the given column header as special field. """
        dump_reader_.set_special_column( self.handle, header, special_field )
