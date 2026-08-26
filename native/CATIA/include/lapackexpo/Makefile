FC	= gfortran
FFLAGS	= -O3

RM=rm -f
# Among the 3 possibilities below, uncomment the appropriate
# case for your environment and comment the others.

# case 1: works when LAPACK and BLAS are installed.
OBJLIBS =
LIBS    =  -llapack -lblas
RANLIB = ranlib 
AR=ar 
PLATFORM= $(shell uname -s)
MACHINETYPE=$(shell uname -p)


# case 2: works when LAPACK is not installed but BLAS is.
#LIBS    = -lblas
#OBJLIBS = lapack.o

# case 3: works when neither LAPACK nor BLAS are installed.
#OBJLIBS = blas.o lapack.o
#LIBS	=

#############
#############

OBJS   = expokit.o

all:
	$(FC) -c expokit.f -o $(OBJS)_$(PLATFORM)
	$(AR) -r libexpokit_$(PLATFORM).a $(OBJS)
	$(RANLIB) libexpokit_$(PLATFORM).a

clean:
	$(RM) $(EXPOKIT)
	$(RM) *.a


