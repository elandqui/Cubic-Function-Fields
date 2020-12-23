# This is the makefile to create the library libhbounds.a
# type 'make' or make hbounds'
#
#-----------------
# COMPILER FLAGS:
#-----------------

CC = g++
# -g   debugging output (gdb)

INCLUDE_DIR = -I. -I/home/research/elandqui/Software/include
# include directories
LIBS_DIR = -L/home/research/elandqui/Software/lib
# directories containing the libraries

LIBS = -lntl -lgmp -lm
# library includes 

OBJ= invariants.o
# object files for library

SRC= invariants.cc 
# Source files for library

CFLAGS = -c  -O3 -Wall -Wno-deprecated $(INCLUDE_DIR)
# -g   debugging output (gdb)

LDFLAGS = $(INCLUDE_DIR) 
# -g   debugging output (gdb)

AR=ar
# command to make a library

ARFLAGS=cr 
# arguments for AR

libinvariants: $(OBJ)  
	$(AR) $(ARFLAGS) libinvariants.a $(OBJ)
	ranlib libinvariants.a

# first line is dependency list
# second line is compile command
invariants.o: invariants.cc invariants.h 
	$(CC) $(CFLAGS) invariants.cc 

clean: 
	rm *.o 

