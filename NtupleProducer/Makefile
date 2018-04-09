CMSINC      = -I$(CMSSW_BASE)/src

#EGAMLIB      = -L$(CMSSW_BASE)/lib/slc6_amd64_gcc481 -lEgammaAnalysisElectronTools
#FORMLIB      = -L$(CMSSW_BASE)/lib/slc6_amd64_gcc481 -lEgammaAnalysisElectronTools

JETLIBS     = -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc630 -lCondFormatsJetMETObjects -lFWCoreMessageLogger
JETINC      = -I$(CMSSW_RELEASE_BASE)/src

BOOSTSYS       = $(shell scram tool tag boost BOOST_BASE)
BOOSTINC       = -I$(BOOSTSYS)/include

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++ -std=c++11 -fPIC -Wno-write-strings -Wno-writable-strings -Wno-attributes 
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS       = $(ROOTCFLAGS)
INCLUDE_FLAGS  = $(JETINC) $(BOOSTINC)
LDLIBS         = $(JETLIBS) $(ROOTLIBS)
GLIBS          = $(ROOTGLIBS)

EXE           = test/NtupleProducer

INC 	      = include/NtupleProducer.h include/Tree.h include/Base.h \
                include/Electron.h include/Muon.h include/Ntuple.h \
		include/Event.h include/Jet.h include/Truth.h include/Tau.h include/GenJet.h \
		include/TriggerObj.h

SRC	      = src/NtupleProducer.cxx src/Tree.cxx src/Base.cxx \
                src/Electron.cxx src/Muon.cxx src/Ntuple.cxx \
		src/Event.cxx src/Jet.cxx src/Truth.cxx src/Tau.cxx src/GenJet.cxx \
		src/TriggerObj.cxx

OBJS          = obj/NtupleProducer.o obj/Tree.o obj/Base.o \
                obj/Electron.o obj/Muon.o obj/Ntuple.o \
		obj/Event.o obj/Jet.o obj/Truth.o obj/Tau.o obj/GenJet.o \
		obj/TriggerObj.o

LIB           = libNtuple.so

all: 	      $(LIB) $(EXE)

$(LIB):	      $(INC) $(SRC)
	      @echo "####### Generating dictionary"
	      @rootcint -f obj/NtupleProducerDict.cxx -c -p $(CXXFLAGS) \
	      $(INCLUDE_FLAGS) -I. $(INC) include/LinkDef.h

	      @echo "####### Building library $(LIB)"
	      @$(CXX) $(SOFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDE_FLAGS) -I. $(SRC) \
	      obj/NtupleProducerDict.cxx -o $(LIB) $(ROOTLIBS) 
	      
	      @echo  "####### Removing generated dictionary"
	      @rm -f obj/NtupleProducerDict.cxx obj/NtupleProducerDict.h
	      @rm -f *.o

$(EXE):	      $(LIB) src/NtupleProducer.cxx
	      @echo "####### Building object file for executable"
	      @$(CXX) -c $(CXXFLAGS) $(INCLUDE_FLAGS) src/NtupleProducer.cxx -o obj/NtupleProducer.o
	      @echo "####### Building executable"
	      @$(CXX) $(CXXFLAGS) obj/NtupleProducer.o -o $(EXE) $(LDLIBS) $(LIB)

clean:
	      @rm -f $(OBJS) $(EXE) obj/NtupleProducerDict.cxx obj/NtupleProducerDict.h $(LIB)
