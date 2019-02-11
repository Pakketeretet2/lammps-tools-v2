CC = g++
LINK = $(CC) -llammpstools

LAMMPS_TOOLS_DIR = $(HOME)/projects/lammps-tools-v2/cpp_lib
CLARA_DIR = $(HOME)/projects/Clara

FLAGS = -O3 -g -I$(LAMMPS_TOOLS_DIR) -I$(CLARA_DIR)/single_include/

.PHONY = all

all: get_completion

get_completion: get_completion.o
	$(LINK) $^ -o $@

get_completion.o:
	$(CC) $(FLAGS) -c get_completion.cpp -o get_completion.o
