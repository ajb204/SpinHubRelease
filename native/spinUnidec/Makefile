#
# Makefile for decon
#
#           (C) A.Baldwin
#           University of Oxford
#           andrew.baldwin@chem.ox.ac.uk
#

CC = g++
MPICC = g++

# OpenMP link mode.  On macOS, STATIC=True means only libomp is linked
# statically; libc++ and libSystem remain normal macOS system dylibs.
# OMP_LINK can be overridden directly: make OMP_LINK=static|dynamic
STATIC ?= False
OMP_LINK ?= $(if $(filter True true TRUE 1 yes YES,$(STATIC)),static,dynamic)
PARALLEL = True
DOUBLE = False
FLOAT = True

VERSION = 1.0

PLATFORM = $(shell uname -p)
MACHINETYPE = $(shell uname -s)

SRC=src/
EXTERN_LIBS ?= ../../extern/libs
BIN ?= ../../bin/

$(info Platform: $(PLATFORM))
$(info Arch    : $(MACHINETYPE))
$(info Version : $(VERSION))
$(info Compiler: $(CC))
$(info OpenMP : $(OMP_LINK))

CCFLAGS = -O3 -std=c++11 -Isrc -Wno-vla-cxx-extension -fopenmp

# FFTW runtime libraries available in this environment are versioned only.
# Use exact shared object names so the build works without development symlinks.
ifeq ($(FLOAT),True)
  LDFLAGS += $(EXTERN_LIBS)/libfftw3f.a
  ifeq ($(PARALLEL),True)
    LDFLAGS += $(EXTERN_LIBS)/libfftw3f_threads.a
  endif
endif

ifeq ($(DOUBLE),True)
  CCFLAGS += -DDOUBLE2D -DDOUBLE3D -DDOUBLE1D -DDOUBLE4D
  LDFLAGS +=  $(EXTERN_LIBS)/libfftw3.a
  ifeq ($(PARALLEL),True)
    LDFLAGS += $(EXTERN_LIBS)/libfftw3_threads.a
  endif
endif

ifeq ($(MACHINETYPE),Darwin)
	CCFLAGS += -DMAC_OS_X -Wno-switch-bool -Wno-return-type -Wno-vla-cxx-extension
	ifeq ($(PLATFORM),arm)
		LLVM_PATH= /opt/homebrew/opt/llvm/
	else
		LLVM_PATH= /usr/local/opt/llvm
	endif
	CC = $(LLVM_PATH)/bin/clang++

	# Homebrew libomp ships both libomp.dylib and libomp.a.  Using the
	# archive by full path prevents ld from selecting the dylib.
	LIBOMP_PREFIX ?= /opt/homebrew/opt/libomp
	LIBOMP_STATIC ?= $(LIBOMP_PREFIX)/lib/libomp.a
	ifeq ($(OMP_LINK),static)
		LDFLAGS += $(LIBOMP_STATIC)
	else ifeq ($(OMP_LINK),dynamic)
		LDFLAGS += -L$(LIBOMP_PREFIX)/lib -lomp
	else
		$(error OMP_LINK must be either 'static' or 'dynamic')
	endif

else
    STATIC = True
    ifeq ($(MACHINETYPE),Linux)
		CCFLAGS += -DLINUX
    else
	   	UNAME_S := "Windows"
	   	CCFLAGS += -DWINDOWS 
	   	LDFLAGS := -lfftw3 -lfftw3f -lpthread -lm -static-libgcc -static-libstdc++ -o
    endif
endif

# Preserve the historical full-static behaviour away from macOS.
# macOS does not support a fully static executable, so Darwin deliberately
# never receives -static.
ifneq ($(MACHINETYPE),Darwin)
  ifneq ($(filter True true TRUE 1 yes YES,$(STATIC)),)
    CCFLAGS += -static
  endif
endif


OBJS = \
	lib/general_${PLATFORM}.o \
	lib/slice1DUSTA_$(PLATFORM).o  \
	lib/slice1D_$(PLATFORM).o  \
	lib/slice2D_$(PLATFORM).o  \
	lib/slice2Dfit_$(PLATFORM).o  \
	lib/slice3D_$(PLATFORM).o \
	lib/slice4D_$(PLATFORM).o    \
	lib/pipeClass_$(PLATFORM).o \
	lib/deconParse_$(PLATFORM).o \
	lib/deconMain_$(PLATFORM).o   

#lib/decon_parallel_$(PLATFORM).o   
#lib/decon_$(PLATFORM).o   \




TARGETS=$(BIN)decon_parallel_$(MACHINETYPE)_$(PLATFORM) $(BIN)decon_$(MACHINETYPE)_$(PLATFORM)
#$(BIN)decon_parallel_$(UNAME_S)_$(UNAME_P)

.PHONY: all static dynamic clean install

# Convenience entry points.  These produce the same binary names; use one
# mode at a time.  Equivalent forms are `make OMP_LINK=static` and
# `make OMP_LINK=dynamic`.
static:
	$(MAKE) clean
	$(MAKE) OMP_LINK=static STATIC=True all

dynamic:
	$(MAKE) clean
	$(MAKE) OMP_LINK=dynamic STATIC=False all

all: $(BIN) $(TARGETS)
	@echo "================================================================="
	@echo "   unidecNMRdecon - version "$(VERSION) 
	@echo "   Compiled successfully "	
	@echo "                    (C) decon"
	@echo "                        A.Baldwin"
	@echo "                        C.Buchanan"
	@echo "                        University of Oxford  "          
	@echo "                        andrew.baldwin@chem.ox.ac.uk "         
	@echo "================================================================="

$(BIN):
	mkdir -p $(BIN)

#$(BIN)decon_$(UNAME_S)_$(UNAME_P): $(SRC)decon.cpp $(SRC)deconMain.hpp $(SRC)deconParse.hpp $(SRC)slice.hpp $(SRC)pipeClass.hpp $(SRC)slice1D.hpp $(SRC)slice2D.hpp $(SRC)slice3D.hpp $(SRC)slice4D.hpp $(SRC)fftw3.h | $(BIN)
#	$(CC) $(CCFLAGS) $< $(LDFLAGS) -o $@

#$(BIN)decon_parallel_$(UNAME_S)_$(UNAME_P): $(SRC)decon_parallel.cpp $(SRC)deconMain.hpp $(SRC)deconParse.hpp $(SRC)slice.hpp $(SRC)pipeClass.hpp $(SRC)slice1D.hpp $(SRC)slice2D.hpp $(SRC)slice3D.hpp $(SRC)slice4D.hpp $(SRC)fftw3.h | $(BIN)
#	$(MPICC) $(CCFLAGS) $< $(LDFLAGS) -o $@


# ordinary files (c++)
lib/%_$(PLATFORM).o: src/%.cpp
	$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@



#do final linking.
$(BIN)decon_parallel_$(MACHINETYPE)_${PLATFORM}: $(OBJS) lib/decon_parallel_${PLATFORM}.o 
	$(CC)  $^ ${LDFLAGS} -o  $@


#do final linking.
$(BIN)decon_$(MACHINETYPE)_${PLATFORM}: $(OBJS) lib/decon_${PLATFORM}.o 
	$(CC)  $^ ${LDFLAGS} -o  $@




clean:
	@echo "Removing objects"
	rm -f lib/*.o
	rm -f $(TARGETS)

install: $(TARGETS)
	@echo "Compiling binaries" $< 
