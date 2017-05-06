CC = g++
INC = -I../../../c_interface/
LNK = -llammpstools
FLAGS = -O3 -std=c++11 -shared -fPIC -llammpstools `python-config --cflags --ldflags`
FLAGS += $(INC)
FLAGS += $(LNK)

OBJS = dumpreader.so block_data.so

.PHONY = clean all

all: $(OBJS) $(LIB)

%.so: %_python.cpp
	$(CC) $(FLAGS) $^ -o $@

clean:
	rm -f $(OBJS)
