#------------------------------------------------------------------------------
# ROOT definitions
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lMathMore  -lMinuit2

#------------------------------------------------------------------------------
# g++ options

CXX = g++
CXXFLAGS = -g -O -fPIC -Wall -Wno-deprecated  
LD = g++

SOFLAGS = -shared #-dynamiclib #-single_module #-undefined dynamic_lookup

CXXFLAGS     += $(ROOTCFLAGS) 

LIBS         = ${ROOTGLIBS} ${ROOTLIBS} -L/home-pbs/nchanon/CMSSW_8_0_20/src/IPHCNtuple/MEM/Minimizer/SubGradient_h.so 
#------------------------------------------------------------------------------
# Compilation

all : testMinimizer

clean:
	rm *.o
	rm test*

objects= TestMinimizer.o

testMinimizer: $(objects) 
	$(CXX) -o testMinimizer TestMinimizer.o ${LIBS}  $(CXXFLAGS) #${SOFLAGS}
#	$(CXX) -o test $(objects) ${ROOTLIBS} ${ROOTGLIBS} $(CXXFLAGS) ${SOFLAGS} 
#/home-pbs/nchanon/CMSSW_8_0_20/src/IPHCNtuple/MEM/Minimizer/libSubGradient.so

TestMinimizer.o: TestMinimizer.cpp Functions.h SubGradient.h
	$(CXX) ${LIBS} ${CXXFLAGS} -c TestMinimizer.cpp 

#SubGradient_h.so: SubGradient.h
#	$(CXX) ${LIBS} ${CXXFLAGS} -o $@ ${SOFLAGS}

#SubGradient.o: SubGradient.h
#	$(CXX) ${LIBS} ${CXXFLAGS} -c SubGradient.h
