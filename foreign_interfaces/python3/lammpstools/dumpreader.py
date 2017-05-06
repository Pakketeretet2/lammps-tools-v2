import dumpreader
import block_data

if __name__ == "__main__":

    print("Testing dumpreader")
    dr = dumpreader.new_dump_reader( "../../../test/lammps_dump_file_test.dump.bin", 2, 0)

    headers = [ 'id', 'type', 'x', 'y', 'z', 'c_pe' ]
    dumpreader.print_list_of_strings( headers );
    dumpreader.set_column_headers( dr, 6, headers )
    
    bd = block_data.block_data_handle()
    status = 0

    while dumpreader.get_next_block( dr, bd ) == 0:
        print("At t = ", bd.time_step(), ", N = ", bd.n_atoms())
        size = 0
        t = 0
        xf = block_data.data_as_double( bd, "x", size, t )
        print("size is ", size,", t = ", t)

    dumpreader.delete_dump_reader(dr)

