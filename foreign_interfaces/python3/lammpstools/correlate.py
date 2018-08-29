import correlation_
import data_field

def correlate( b, data, data_type, x0, x1, dx, dims ):
    """ Correlates data over distances. """

    if data_type == data_field.DATA_TYPE_INT:
        return correlation_.correlate_int( b.get_ref_(), data, x0, x1, dx, dims )
    elif data_type == data_field.DATA_TYPE_DOUBLE:
        return correlation_.correlate_double( b.get_ref_(), data, x0, x1, dx, dims )
    else:
        print( "Unknown data type", data_type, file = sys.stderr )
        return None


def correlate( b, df, x0, x1, dx, dims ):
    if df.type() == data_field.DATA_TYPE_INT:
        return correlation_.correlate_int( b.get_ref_(), df.get_data(), x0, x1, dx, dims )
    elif df.type() == data_field.DATA_TYPE_DOUBLE:
        return correlation_.correlate_double( b.get_ref_(), df.get_data(), x0, x1, dx, dims )
    else:
        print( "Unknown data type for data field", df.name(), file = sys.stderr )
        return None
