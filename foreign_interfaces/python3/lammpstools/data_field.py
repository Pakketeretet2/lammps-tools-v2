import data_field_

DATA_TYPE_DOUBLE = data_field_.TYPES.DOUBLE
DATA_TYPE_INT = data_field_.TYPES.INT

class data_field:
    """ Wrapper around data_field """
    def __init__(self, handle, data_type = DATA_TYPE_DOUBLE):
        self.handle = handle
        self.data_type = data_type

    def size(self):
        return data_field_.get_size( self.handle )

    def __len__(self):
        return self.size()

    def __iter__(self):
        self.internal_iter_counter = 0
        return self

    def __next__(self):
        if self.internal_iter_counter >= self.size():
            raise StopIteration
        
        value = self.__getitem__(self.internal_iter_counter)
        self.internal_iter_counter += 1
        return value

    def type(self):
        return data_field_.get_type( self.handle )

    def name(self):
        return data_field_.get_name( self.handle )

    def set_name(self, name):
        data_field_.set_name( self.handle, name )

    def set_size(self, size):
        data_field_.set_size( self.handle, size )

    def get_data(self):
        """ Returns the raw data interpreted as the proper type. """
        if self.data_type == DATA_TYPE_DOUBLE:
            return data_field_.as_float(self.handle)
        elif self.data_type  == data_field_.TYPES.INT:
            return data_field_.as_int(self.handle)
        else:
            raise RuntimeError("Unkown data type encountered!")

    def __getitem__(self, index):
        """ Indexes internal data. """
        if index < 0 or index >= self.size():
            raise RuntimeError("Index out of bounds!")

        if self.data_type == DATA_TYPE_INT:
            return data_field_.get_indexed_int_data(self.handle, index)
        elif self.data_type == DATA_TYPE_DOUBLE:
            return data_field_.get_indexed_double_data(self.handle, index)
        else:
            raise RuntimeError("Unknown data type encountered!")

    def __setitem__(self, index, value):
        """ Indexes internal data. """
        if index < 0 or index >= self.size():
            raise RuntimeError("Index out of bounds!")
        s = 0
        if self.data_type == DATA_TYPE_INT:
            s = data_field_.set_indexed_int_data(self.handle, index, value)
        elif self.data_type == DATA_TYPE_DOUBLE:
            s = data_field_.set_indexed_double_data(self.handle, index, value)
        else:
            raise RuntimeError("Unknown data type encountered!")
        if s != 0:
            raise RuntimeError("Attempted to modify data in const data_field!")

    def __iadd__(self, index, number):
        """ In/decrements internal data. """
        if index < 0 or index >= self.size():
            raise RuntimeError("Index out of bounds!")

        if self.data_type == DATA_TYPE_INT:
            old_val = data_field_.get_indexed_int_data(self.handle, index)
            new_val = old_val + number
            data_field_.set_indexed_int_data(self.handle, index, new_val)
        elif self.data_type == DATA_TYPE_DOUBLE:
            old_val = data_field_.get_indexed_double_data(self.handle, index)
            new_val = old_val + number
            data_field_.set_indexed_double_data(self.handle, index, new_val)
        else:
            raise RuntimeError("Unknown data type encountered!")

def new_data_field(name, dtype, size):
    """ Makes a freshly instantiated data_field and returns a handle. """
    return data_field( data_field_.new_data_field( name, dtype, size ), dtype )
