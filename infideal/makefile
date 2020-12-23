CC = g++
# -g   debugging output (gdb)

INCLUDE_DIR = -I. -I.. -I/home/research/elandqui/Software/include
# include directories
LIBS_DIR = -L../CubicIdeal -L../Invariants -L/home/research/elandqui/Software/lib
# directories containing the libraries

LIBS = -libcubic -linvariants -lntl -lgmp -lm
# library includes 

OBJ= infrastructure_ideal.o
# object files for library

SRC= infrastructure_ideal.cc 
# Source files for library

CFLAGS = -c -O3 -Wall -Wno-deprecated $(INCLUDE_DIR)
# -g   debugging output (gdb)

LDFLAGS = $(INCLUDE_DIR) 
# -g   debugging output (gdb)

AR=ar
# command to make a library

ARFLAGS=cr 
# arguments for AR

libinfideal: $(OBJ)  
	$(AR) $(ARFLAGS) libinfideal.a $(OBJ)
	ranlib libinfideal.a

# first line is dependency list
# second line is compile command
infrastructure_ideal.o: infrastructure_ideal.cc infrastructure_ideal.h 
	$(CC) $(CFLAGS) infrastructure_ideal.cc 

clean: 
	rm *.o 
