import data_field_


class data_field:
    """ Wrapper around data_field """
    def __init__(self, handle):
        self.handle = handle

    def size(self):
        return data_field_.get_size( self.handle )

    def type(self):
        return data_field_.get_type( self.handle )

    def name(self):
        return data_field_.get_name( self.handle )
