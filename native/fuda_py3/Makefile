######################################
# Makefile for FUDA
#
# 1) compile minpack
# 2) compile fuda with minpack
# 3) compile nmrPipe modules
# 4) compile fuda dataIO with nmrPipe
# 5) transfer results into fuda folder
#
# A.Baldwin 14th Sept 2022
# Streamlined 25th March 2024
######################################

# Set compiler for fuda and dataIO (C++)
CC = g++
#fortran compiler  (FORTRAN)
FC=gfortran      
#Compiler for nmrPipe  (C)
CN = gcc
#executable run to start python.
PYTHON_EXE ?= python3

# Setup variables constructed from platform.
OBJ_PATH = py/
SRC = src/
BIN ?= ../../bin/
IO = $(SRC)io/

#####  Final objects, created by this script #####
#####  transferred to install_path           #####
PYTHON_INST    = fuda.py
PYTHON_IO_INST = fudaIO.py
PYTHON_SO      = $(OBJ_PATH)fudalib.so
PYTHON_SO_IO   = $(OBJ_PATH)dataIO.so
FUDA_LIB       = $(OBJ_PATH)libfuda.a 

###### distribution NMRPipe locations. ########
NMRPIPE_PATH = $(SRC)nmrPipe_format

#####  Python location and version. ########
PYTHON_PREFIX = $(shell $(PYTHON_EXE) -c 'import sys; print (sys.prefix)')
PYTHON_VERSION = $(shell $(PYTHON_EXE) -c 'import sys; print (sys.version.split(".")[0]+"."+sys.version.split(".")[1])') 
PYTHON_INCLUDE = $(shell $(PYTHON_EXE) -c 'import sysconfig; print(sysconfig.get_path("include"))')
NUMPY_INCLUDE = $(shell $(PYTHON_EXE) -c 'import numpy; print(numpy.get_include())')
PYTHONLDFLAG = -L$(PYTHON_PREFIX)/lib 
$(info PYTHON_PREFIX: $(PYTHON_PREFIX))
$(info PYTHON_VERSION: $(PYTHON_VERSION))
$(info PYTHON_INCLUDE: $(PYTHON_INCLUDE))
$(info NUMPY_INCLUDE: $(NUMPY_INCLUDE))
$(info PYTHONLDFLAG: $(PYTHONLDFLAG))

### numpy/arrayobject.h header ########
# have copied in the headers to include: t
# these will be linked to the python version
# cross your fingers that this works.
##DISABLED
#PYTHON_INCLUDE2 = include

##### SYSTEM SETTINGS ########################
PLATFORM = $(shell uname -s)
CPU = $(shell uname -p)

#### MINPACK SETTINGS #################
#folder for minpack
FFLAGS=-O3
LIBMINPACK = $(SRC)libminpack.a
OBJF := $(patsubst $(SRC)minpack/%.f,$(SRC)minpack/%.o,$(wildcard $(SRC)minpack/*.f))

######## Fuda library.  ########################
DATAIO_LIB   = $(IO)libdataIO.a # dataIO library
OBJ := $(patsubst $(SRC)C/%.C,$(SRC)%.o,$(wildcard $(SRC)C/*.C))

############# NMRPIPE INSTALL ###################
#stripped down pipe to get enough for showhdr
NMRPIPE_EXTRA_OBJ := $(patsubst $(NMRPIPE_PATH)/%.c,$(NMRPIPE_PATH)/%.o,$(wildcard $(NMRPIPE_PATH)/*.c))

##########################################################
#Objects for fuda data IO
NMRPIPE_OBJ  = $(IO)dataIO_nmrpipe.o
SPARKY_OBJ = $(IO)dataIO_sparky.o
DATAIO_OBJ = $(IO)dataIO.o
OBJ_IO = $(DATAIO_OBJ) $(NMRPIPE_OBJ) $(SPARKY_OBJ)

####### FIGURE OUT FLAGS FOR FUDA COMPILE #############################

# Setup platform dependent variables.
ifeq ($(PLATFORM),Linux)
#	ifeq($(CPU),i686)
		CFLAGS = -pthread -fno-strict-aliasing \
			 -DNDEBUG -g -O3 -Wall -fPIC
		LDFLAGS = -pthread $(FUDA_LIB) $(LIBMINPACK) -lm
		CFLAGS_PYTHON = -pthread -fno-strict-aliasing \
			 -DNDEBUG -g -O3 -Wall -fPIC \
			 -I. -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -DHAVE_CONFIG_H
		LDFLAGS_PYTHON = -pthread -shared $(FUDA_LIB) $(LIBMINPACK) -lm
#	endif
#	ifeq($(CPU),x86_64)
#		CFLAGS = -pthread -fno-strict-aliasing \
#		         -DNDEBUG -g -O3 -Wall -fPIC
#		LDFLAGS = -pthread $(FUDA_LIB) $(LIBMINPACK) -lm
#		CFLAGS_PYTHON = -pthread -fno-strict-aliasing \
#	  	          -DNDEBUG -g -O3 -Wall -fPIC \
#		          -I. -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -DHAVE_CONFIG_H
#		LDFLAGS_PYTHON = -pthread -shared $(FUDA_LIB) $(LIBMINPACK) -lm
#	endif
endif

#LLVM_PATH= $(CONDA_PREFIX)
ifeq ($(PLATFORM),Darwin)
	ifeq ($(CPU),i386)
		#LLVM_PATH= /usr/local/opt/llvm
		CFLAGS = \
			-O3 -Wall -Wstrict-prototypes  -Wno-deprecated-declarations -Wno-constant-logical-operand -Wno-unused-but-set-variable -Wno-string-plus-int -Wno-format -Wno-c++11-compat-deprecated-writable-strings -Wno-\#warnings -Wno-dangling-else -Wno-unused-function -Wno-return-type
		#CFLAGS += -DNDEBUG -g 
		#	-fno-strict-aliasing  \
		#	-no-cpp-precomp  \

		LDFLAGS = $(FUDA_LIB) $(LIBMINPACK) -lm
		CFLAGS_PYTHON = \
			-O3 -Wall -Wstrict-prototypes   -Wno-deprecated-declarations -Wno-constant-logical-operand -Wno-unused-but-set-variable -Wno-string-plus-int -Wno-format -Wno-c++11-compat-deprecated-writable-strings -Wno-\#warnings -Wno-dangling-else -Wno-unused-function -Wno-return-type -Wno-writable-strings  -Wno-vla-cxx-extension \
			-I.src -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)
		#-I$(PYTHON_INCLUDE2)
		#CFLAGS_PYTHON += -DNDEBUG -g 
		#			-fno-strict-aliasing  \
		#	-no-cpp-precomp  
		LDFLAGS_PYTHON = -bundle -undefined dynamic_lookup \
			 -ldl -framework CoreFoundation $(FUDA_LIB) $(LIBMINPACK) -lm
	endif
	ifeq ($(CPU),arm)
		#LLVM_PATH= /opt/homebrew/opt/llvm/
		CXXFLAGS = -mmacosx-version-min=11.1
		CFLAGS = \
			-fno-strict-aliasing \
			-DNDEBUG -O3 -Wall 
		LDFLAGS = $(FUDA_LIB) $(LIBMINPACK) -lm
		CFLAGS_PYTHON = -fno-strict-aliasing -fno-common -DNDEBUG -O3 -I. -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -DTHREADING -DNMR64
		#remove include2
		LDFLAGS_PYTHON = -bundle -undefined dynamic_lookup -ldl -framework CoreFoundation $(FUDA_LIB) $(LIBMINPACK) -lm
	endif
	#CC = $(LLVM_PATH)/bin/clang++
	#CC = g++

endif
ifeq ($(PLATFORM),IRIX)
	CFLAGS = -xansi
	LDFLAGS =
endif
ifeq ($(PLATFORM),IRIX64)
	CFLAGS = -xansi
	LDFLAGS =
endif
ifeq ($(PLATFORM),OSF1)
	CFLAGS = -std1
	LDFLAGS =
endif




####### FIGURE OUT FLAGS FOR DATAIO COMPILE #############################
# Setup platform dependent variables.

ifeq ($(PLATFORM),Linux)
	CFLAGS_IO = -g
	LDFLAGS_IO = 
	CFLAGS_PYTHON_IO = -fPIC -g -O2 -I. -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -DHAVE_CONFIG_H
	LDFLAGS_PYTHON_IO = -shared
endif
ifeq ($(PLATFORM),Darwin)
	ifeq ($(CPU),i386)
		CFLAGS_IO = -O3
		LDFLAGS_IO = 
		#CFLAGS_PYTHON_IO = -fno-strict-aliasing -no-cpp-precomp -fno-common -dynamic -DNDEBUG -g -O3  -I. -I/sw/include -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)  -I$(PYTHON_INCLUDE2)  -DTHREADING
		#CFLAGS_PYTHON_IO = -DMAC_OSX -fno-strict-aliasing -no-cpp-precomp -fno-common -dynamic -DNDEBUG -g -O3  -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)  -I$(PYTHON_INCLUDE2)  -DTHREADING -Wno-deprecated-declarations -Wno-constant-logical-operand -Wno-unused-but-set-variable -Wno-string-plus-int -Wno-format -Wno-c++11-compat-deprecated-writable-strings -Wno-\#warnings -Wno-dangling-else -Wno-unused-function -Wno-return-type -Wno-delete-abstract-non-virtual-dtor
		CFLAGS_PYTHON_IO = -DMAC_OSX -fno-strict-aliasing -no-cpp-precomp -fno-common -O3  -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)  -DTHREADING -Wno-deprecated-declarations -Wno-constant-logical-operand -Wno-unused-but-set-variable -Wno-string-plus-int -Wno-format -Wno-c++11-compat-deprecated-writable-strings -Wno-\#warnings -Wno-dangling-else -Wno-unused-function -Wno-return-type -Wno-delete-abstract-non-virtual-dtor
		#removed include2  -I$(PYTHON_INCLUDE2)
		#LDFLAGS_PYTHON_IO = -Wl,-F. -bundle -framework Python -L/sw/lib
		#LDFLAGS_PYTHON_IO = -Wl,-F. -bundle $(PYTHONLDFLAG) -lpython2.7
		LDFLAGS_PYTHON_IO = -bundle -undefined dynamic_lookup
	endif
	ifeq ($(CPU),arm)
		CFLAGS_IO = -g  
		LDFLAGS_IO = 
		CFLAGS_PYTHON_IO = -fno-strict-aliasing -fno-common -DNDEBUG -O3 -I. -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) -DTHREADING -DNMR64
		# LDFLAGS_PYTHON_IO = -Wl,-rpath,$(PYTHON_PREFIX)/lib -F. -bundle  $(PYTHONLDFLAG) -lpython$(PYTHON_VERSION) -ldl -framework CoreFoundation -flat_namespace -undefined suppress 
		LDFLAGS_PYTHON_IO = -bundle -undefined dynamic_lookup -ldl -framework CoreFoundation
	endif
endif

ifeq ($(PLATFORM),IRIX)
	CFLAGS_IO = -xansi
	LDFLAGS_IO =
endif
ifeq ($(PLATFORM),IRIX64)
	CFLAGS_IO = -xansi
	LDFLAGS_IO =
endif
ifeq ($(PLATFORM),OSF1)
	CFLAGS_IO = -std1
	LDFLAGS_IO =
endif


########## FLAGS FOR NMRPIPE #############################################

CPFLAGS_PIPE    = -DNO_DMX_SRC -DMAC_OSX -DUSE_NMR_PLUGIN -D__NO_STRING_INLINES -Wno-comment -Wno-parentheses -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -I$(NMRPIPE_PATH) -fwritable-strings -O3



#########################################################################
# Here goes the rules.


#compile the objects: fuda, fudaIO, minPack and nmrPipe
#this can be done with multiple cpus
all : $(OBJ) $(OBJ_IO) $(OBJF)  $(NMRPIPE_EXTRA_OBJ) $(SRC)fuda_python.o $(IO)dataIO_python.o

#libdataIO.a
#must be done in series
$(DATAIO_LIB) : $(OBJ_IO) 
	$(info Creating $@)
	ar ru $@ $?
	ranlib $@

#libfuda.a
#must be done in series
$(FUDA_LIB) : $(OBJ) | $(OBJ_PATH)
	$(info Creating $@)
	ar ru $@ $?
	ranlib $@

#libminpack.a
#must be done in series
$(LIBMINPACK): $(OBJF)
	$(info Creating $@)
	ar ru $@ $?
	ranlib $@

#fudalib.so
#compile static libraries
$(PYTHON_SO) : | $(OBJ_PATH)
$(PYTHON_SO) : $(FUDA_LIB) $(LIBMINPACK) $(SRC)fuda_python.o
	$(info Creating $@)
	$(CC) $(SRC)fuda_python.o $(LDFLAGS_PYTHON) -o $@

#io/dataIO.so
$(PYTHON_SO_IO) : | $(OBJ_PATH)
$(PYTHON_SO_IO) :  $(IO)dataIO_python.o  $(DATAIO_LIB) $(NMRPIPE_EXTRA_OBJS)
	$(CC) $(LDFLAGS_PYTHON_IO) $(IO)dataIO_python.o \
		$(DATAIO_LIB) $(NMRPIPE_EXTRA_OBJ) -lm -o $@

#compile showhdr
$(BIN)showhdr: $(NMRPIPE_EXTRA_OBJ)
	$(info Compiling $@)
	$(CN) $(NMRPIPE_EXTRA_OBJ) -o $@

#make objects: io/dataIO_python
$(IO)dataIO_python.o : $(IO)dataIO_python.C  
	$(CC) $(CFLAGS_PYTHON_IO) -c $< -o $@

#make fuda_python object
$(SRC)fuda_python.o : $(SRC)C/fuda_python.C 
	$(CC) $(CFLAGS_PYTHON) -c $< -o $@

#make objects: io/dataIO.C dataIO_nmrpipe.C and dataIO_sparky.C to get io/dataIO*.o 
$(IO)data%.o: $(IO)data%.C
	$(CC) $(CFLAGS_IO) -c $< -o $@

#make objects: compile C/fuda*.C to get fuda*.o 
$(SRC)fuda_%.o: $(SRC)C/fuda_%.C 
	$(CC) $(CFLAGS) -c $< -o $@ 

#make objects: compile nmrPipe objects
$(NMRPIPE_EXTRA_OBJ) : 
	$(CN) $(CPFLAGS_PIPE)  -c $(patsubst %.o,%.c,$@)  -o $@

# Clean rules.
clean :
	rm -f $(SRC)*.o $(SRC)*.a $(SRC)*.so $(SRC)*/*.o $(SRC)*/*.a $(SRC)*/*.so $(OBJ_PATH)* $(BIN)showhdr $(BIN)nmrPipeFit.py

# Clean rules.
cleanInst : 
	rm -f $(SRC)*.o $(SRC)*.a $(SRC)*.so $(SRC)*/*.o $(SRC)*/*.a $(SRC)*/*.so

# Output directories must be recreated after distclean.
$(OBJ_PATH):
	mkdir -p $@

$(BIN):
	mkdir -p $@

# Install rules. Copy results into place, assemble objects/compile as needed
install : $(BIN)showhdr $(FUDA_LIB) $(PYTHON_SO) $(PYTHON_SO_IO) $(OBJ_PATH)$(PYTHON_INST) $(OBJ_PATH)$(PYTHON_IO_INST) $(BIN)nmrPipeFit.py

$(OBJ_PATH)$(PYTHON_INST) : $(SRC)$(PYTHON_INST) | $(OBJ_PATH)
	cp $< $@
	#done

$(OBJ_PATH)$(PYTHON_IO_INST) : $(IO)$(PYTHON_IO_INST) | $(OBJ_PATH)
	cp $< $@
	#done

$(BIN)nmrPipeFit.py : bin/nmrPipeFit.py | $(BIN)
	cp $< $@
	#done



