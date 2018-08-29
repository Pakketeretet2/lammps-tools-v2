import histogram_, block_data_


def make_histogram( data, y0, y1, N_bins ):
    """ Make histogram of given data in interval [y0, y1) with N_bins bins. """

    if type(data) == block_data_.VectorDouble:
        return histogram_.make_histogram_double( data, y0, y1, N_bins )
        pass
    elif type(data) == block_data_.VectorInt:
        return histogram_.make_histogram_int( data, y0, y1, N_bins )
        pass
    else:
        print("Cannot histogram given data!", file = sys.stderr)
        return None
