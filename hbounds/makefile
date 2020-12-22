# This is the makefile to create the library libhbounds.a
# type 'make' or make hbounds'
#
#-----------------
# COMPILER FLAGS:
#-----------------

CC = g++
# -g   debugging output (gdb)

INCLUDE_DIR = -I. -I/home/elandqui/include/
# include directories
LIBS_DIR = -L../Invariants -L/home/elandqui/lib/  
# directories containing the libraries

LIBS = -linvariants -lntl -lgmp -lm
# library includes 

OBJ= hbounds.o
# object files for library

SRC= hbounds.cc 
# Source files for library

CFLAGS = -c  -O3 -Wall -Wno-deprecated $(INCLUDE_DIR)
# -g   debugging output (gdb)

LDFLAGS = $(INCLUDE_DIR) 
# -g   debugging output (gdb)

AR=ar
# command to make a library

ARFLAGS=cr 
# arguments for AR

libhbounds: $(OBJ)  
	$(AR) $(ARFLAGS) libhbounds.a $(OBJ)
	ranlib libhbounds.a

# first line is dependency list
# second line is compile command
hbounds.o: hbounds.cc hbounds.h 
	$(CC) $(CFLAGS) hbounds.cc 

clean: 
	rm *.o 