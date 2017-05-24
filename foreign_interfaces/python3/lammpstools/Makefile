CC = clang++
# INC = -I../../../cpp_lib/
INC = -I../../../c_interface/ -I../../../internal_include/
LNK = -llammpstools
FLAGS = -O3 -std=c++11 -shared -fPIC `python-config --cflags --ldflags`
FLAGS += $(INC)
FLAGS += $(LNK)

OBJS = block_data_.so dump_reader_.so data_field_.so

.PHONY = clean all

all: $(OBJS) $(LIB)

%_.so: %_python.cpp
	$(CC) $(FLAGS) $^ -o $@

clean:
	rm -f $(OBJS)
