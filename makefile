# This is the makefile to create the library libcubic.a
# type either 'make'  or 'make libcubic' 
# After making a new version of the library the libcubic.a and ideal.h should
# be copied to /u2/mbauer/code/lib and /u2/mbauer/code/include
# To test the library 
# type 'make test' and then './test.out'
# or
# type 'make testsq' and then './testsq.out' 
#
# If there is an error: 'ake: Fatal error: Don't know how to make target
# run dos2unix < makefile > makefile.new
# mv makefile.new makefile
# Don't worry about US keyboard errors and such.
# 
#-----------------
# COMPILER FLAGS:
#-----------------

CC = gcc
# -g   debugging output (gdb)

INCLUDE_DIR = -I. -I/home/research/elandqui/Software/include
# include directories
LIBS_DIR = -L../Invariants -L/home/research/elandqui/Software/lib
# directories containing the libraries

LIBS = -linvariants -lntl -lgmp -lm
# library includes 

OBJ= cubic_ideal.o
# object files for library

SRC= cubic_ideal.cpp 
# Source files for library

CFLAGS = -c  -O3 -Wall -Wno-deprecated $(INCLUDE_DIR)
# -g   debugging output (gdb)

LDFLAGS = $(INCLUDE_DIR) 
# -g   debugging output (gdb)

AR=ar
# command to make a library

ARFLAGS=cr 
# arguments for AR

libcubic: $(OBJ)  
	$(AR) $(ARFLAGS) libcubic.a $(OBJ)
	ranlib libcubic.a

# first line is dependency list
# second line is compile command

cubic_ideal.o: cubic_ideal.cpp cubic_ideal.h 
	$(CC) $(CFLAGS) cubic_ideal.cpp 

clean: 
	rm *.o 
