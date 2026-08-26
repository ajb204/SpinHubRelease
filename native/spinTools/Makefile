#
# Makefile for spinTools
#
#           (C) A.Baldwin
#           University of Oxford
#           andrew.baldwin@chem.ox.ac.uk
#


CC = g++

LDFLAGS = -lm
CFLAGS = -O3
RM = rm -f

#location of source code
SRC := ./src/
#final destination for executables
BIN ?= ../../bin/

#assemble lists of objects and executables automatically
OBJS := $(patsubst $(SRC)%.cpp,$(SRC)%.o,$(wildcard $(SRC)*.cpp))
EXES := $(patsubst $(SRC)%.cpp,$(BIN)%.out,$(wildcard $(SRC)*.cpp))


#.PHONY: prog install
all: $(OBJS) 
	@echo "Compling spinhub Tool objects"

install: $(EXES)

$(BIN)%.out: $(SRC)%.o
	$(CC)  $<  $(LDFLAGS) -o $@

$(SRC)%.o: $(SRC)%.cpp
	$(CC)  -c $(CFLAGS) $< -o $@

cleanO:
	$(RM) $(OBJS)

clean:
	$(RM) $(OBJS)
	$(RM) $(EXES)



