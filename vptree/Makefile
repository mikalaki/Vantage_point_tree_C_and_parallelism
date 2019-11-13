# ####################################################################
#
#			  	   C Makefile
#	   Parallels and Distributed Systems
#				     Exercise 1
#
# Author: Michael Karatzas <mikalaki@ece.auth.gr>
# 						   <mikalaki@it.auth.gr>
# Adapted from:
#  Elearning Examples and edit it to fit the needs of the exercise.
#
# ####################################################################
#
# 'make'        build executable file 'main'
# 'make lib'	build the libraries .a
# 'make clean'  removes all .o and executable files
#
# # call everytime
# .PRECIOUS: %.a
# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7

# define compile-time flags #-Wall
CFLAGS = -Wall -O2 -pthread -fcilkplus -fopenmp
# CFLAGS =   -pthread -fcilkplus -fopenmp
# define any directories containing header files
INCLUDES = -I./inc

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:
LDFLAGS =-L./lib

# define any libraries to link into executable:
#   To ling libusb-1.0 :
#   LIBS = -lusb-1.0
LIBS = -lm

# define the source file for the library
SRCLIB = src/vptree

#define
SRC= lib/vptree



# define the different possible executables
TYPES = sequential openmp cilk pthreads

# # define the executable file  name
MAIN = main

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above
#

# # call everytime
# .PRECIOUS: %.a
#
all: lib $(addprefix $(MAIN)_, $(TYPES))

lib: $(addsuffix .a, $(addprefix $(SRCLIB)_, $(TYPES)))

$(MAIN)_%: $(MAIN).c $(SRC)_%.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file)

#.o.a == (%.a : %.o )
.o.a:
	ar rcs $@ $< && mv src/vptree_*.a lib

# .c.o == (%.o : %.c )
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	$(RM) *.o *~ $(addprefix $(MAIN)_, $(TYPES)) lib/vptree_*.a
