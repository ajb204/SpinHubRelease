#
# Makefile for the CATIA cpmg fitting class functions
#
#           (C) D.F. Hansen
#           University of Toronto, Canada
#           flemming@pound.med.utoronto.ca

#
# Add the environment varaible CATIA_DIR to your source
# ie.
#      setenv CATIA_DIR "${HOME}/Programs/Catia"
# 
# NOTE: requires 3 important (standard) system libraries:
#       lapack
#       blas
#       fortran (such as gfortran or libg2c)
# All can be installed via standard repositories
#

FC=gfortran
RANLIB = ranlib 
AR=ar 
PLATFORM= $(shell uname -s)
MACHINETYPE=$(shell uname -p)
RM=rm -f
SRC=src/
BIN ?= ../../bin/
LIB ?= ../../extern/libs/

# Dependency locations are overrideable. Linux normally finds these in the
# compiler defaults; macOS defaults are discovered from Homebrew when present.
BOOST_INC ?= $(shell if command -v brew >/dev/null 2>&1; then echo $$(brew --prefix boost 2>/dev/null)/include; else echo /usr/include; fi)
LAPACK_LIB ?= $(shell if command -v brew >/dev/null 2>&1; then echo $$(brew --prefix lapack 2>/dev/null)/lib; fi)
GFORT_LIB ?= $(shell if command -v brew >/dev/null 2>&1; then dirname $$($(FC) -print-file-name=libgfortran.dylib 2>/dev/null); fi)
LIBS ?= $(LIB)

#
FFLAGS= -O3 -funroll-loops -Ofast -Wno-argument-mismatch -Wno-aggressive-loop-optimizations 
CFLAGS= -O3 -std=c++11 -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP -funroll-loops  -ffast-math

CFLAG_PARSER= -O3 -Iinclude -std=c++11 -ffast-math -I$(BOOST_INC)

INCLUDE= -Iinclude/  -Iinclude/lapack -I$(BOOST_INC)

ifeq ($(PLATFORM),Darwin)
	ifeq ($(MACHINETYPE),arm)
		LLVM_PATH= /opt/homebrew/opt/llvm/
	else
		LLVM_PATH= /usr/local/opt/llvm
	endif
	CC = $(LLVM_PATH)/bin/clang++

	#CXXFLAGS := -Wall -static -static-libgcc -O3 -DMAC_OS_X

	CFLAGS += -DMAX_OS_X
	CFLAGS += -Wno-return-type-c-linkage -Wno-uninitialized -Wno-missing-braces -Wno-unused-variable -Wno-return-type -Wno-c++11-compat-deprecated-writable-strings -DNDEBUG -Wno-unused-but-set-variable -Wno-deprecated-declarations -Wno-reorder-ctor -Wno-self-assign-overloaded -Wno-vla-cxx-extension
	CFLAG_PARSER+= -DMAX_OS_X -Wno-deprecated-declarations


	#will need to likely edit these in the future: need lapack, a fortran library and blas on the system.


endif
ifeq ($(PLATFORM),Linux)
	CC=g++
	CFLAGS := -static -static-libgcc -DNDEBUG -DLINUX -O3
	CFLAG_PARSER=-Iinclude -DLINUX -std=c++11 -Wno-deprecated-declarations
endif

#
OBJS := $(patsubst $(SRC)%.cpp,lib/%_$(PLATFORM).o,$(wildcard $(SRC)*.cpp))
OBJF := $(patsubst include/lapackexpo/%.f,lib/%_$(PLATFORM).o,$(wildcard include/lapackexpo/*.f))


ifeq ($(PLATFORM),Darwin)
LDFLAG_PARSER=-llapack $(if $(GFORT_LIB),-L$(GFORT_LIB)) -lgfortran -lblas -lm
else
LDFLAG_PARSER=$(if $(LAPACK_LIB),-L$(LAPACK_LIB)) -llapack $(if $(GFORT_LIB),-L$(GFORT_LIB)) -lgfortran -lblas -lm
endif

#
##################################################################
#
####
all: lib $(OBJS) $(OBJF) lib/CATIA.o
	@echo "================================================================="
	@echo "   CATIA (CPMG And Trosy Intelligent Analsysis) - version 0.60   "
	@echo " "
	@echo "   The catia c++ class functions were compiled succesfully "
	@echo " "
	@echo "   If you use this object class in your work, please refer to:"
	@echo " " 
	@echo "                    (C) CATIA"
	@echo "                        D. Flemming Hansen"
	@echo "                        University of Toronto, Canada            "
	@echo "                        flemming@pound.med.utoronto.ca           "
	@echo "================================================================="
lib:
	mkdir lib
#
clean: 
	$(RM) $(OBJS)
	$(RM) lib/catia_${PLATFORM}_$(MACHINETYPE).a
	$(RM) catia_${PLATFORM}_${MACHINETYPE}
	$(RM) $(BIN)catia_${PLATFORM}_${MACHINETYPE}

#
clean-dist: clean
	$(RM) catia_$(PLATFORM).a

#
cleanInst: 
	$(RM) $(OBJS) $(OBJF) lib/CATIA.o

#catia library
$(LIBS)catia_$(PLATFORM)_$(MACHINETYPE).a: $(OBJS) $(OBJF)
	$(AR) -ru $@ $(OBJS) $(OBJF)
	$(RANLIB) $@

# ordinary files (c++)
lib/%_$(PLATFORM).o: src/%.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# ordinary files (fortran)
lib/%_$(PLATFORM).o: include/lapackexpo/%.f
	$(FC) $(FFLAGS) -c $< -o $@

#create an object for the parser.
lib/CATIA.o: CATIA.cpp
	$(CC) $(CFLAG_PARSER)  -c $<  -o $@

#create catia parser
$(BIN)catia_${PLATFORM}_${MACHINETYPE}: $(LIBS)catia_$(PLATFORM)_$(MACHINETYPE).a lib/CATIA.o
	$(CC)  lib/CATIA.o  $< $(LDFLAG_PARSER) -o $@

#assemble the catia library, compile the parser
install:  $(BIN)catia_${PLATFORM}_${MACHINETYPE} $(LIBS)catia_$(PLATFORM)_$(MACHINETYPE).a
	@echo "created: " $<


