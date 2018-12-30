CMSINC      = -I$(CMSSW_BASE)/src

#EGAMLIB      = -L$(CMSSW_BASE)/lib/slc6_amd64_gcc481 -lEgammaAnalysisElectronTools
#FORMLIB      = -L$(CMSSW_BASE)/lib/slc6_amd64_gcc481 -lEgammaAnalysisElectronTools

JETLIBS     = -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc630 -lCondFormatsJetMETObjects -lJetMETCorrectionsModules -lFWCoreMessageLogger
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

PILEUP	      = test/Make_PU_Histograms.exe

EXE           = test/NtupleProducer

INC 	      = include/NtupleProducer.h include/Tree.h include/Base.h include/Helper.h \
                include/Electron.h include/Muon.h include/Ntuple.h \
		include/ElectronExt.h include/MuonExt.h \
		include/Event.h include/Jet.h include/Truth.h include/Tau.h include/GenJet.h \
		include/EventExt.h include/JetExt.h include/TruthExt.h include/TauExt.h include/GenJetExt.h \
		include/TriggerObj.h include/TriggerObjExt.h include/Sync.h

SRC	      = src/NtupleProducer.cxx src/Tree.cxx src/Base.cxx src/Helper.cxx \
                src/Electron.cxx src/Muon.cxx src/Ntuple.cxx \
		src/ElectronExt.cxx src/MuonExt.cxx \
		src/Event.cxx src/Jet.cxx src/Truth.cxx src/Tau.cxx src/GenJet.cxx \
		src/EventExt.cxx src/JetExt.cxx src/TruthExt.cxx src/TauExt.cxx src/GenJetExt.cxx \
		src/TriggerObj.cxx src/TriggerObjExt.cxx src/Sync.cxx

OBJS          = obj/NtupleProducer.o obj/Tree.o obj/Base.o obj/Helper.o \
                obj/Electron.o obj/Muon.o obj/Ntuple.o \
		obj/ElectronExt.o obj/MuonExt.o \
		obj/Event.o obj/Jet.o obj/Truth.o obj/Tau.o obj/GenJet.o \
		obj/EventExt.o obj/JetExt.o obj/TruthExt.o obj/TauExt.o obj/GenJetExt.o \
		obj/TriggerObj.o obj/TriggerObjExt.o obj/Sync.o

LIB           = libNtuple.so

all: 	      $(LIB) $(EXE) $(PILEUP)

$(LIB):	      $(INC) $(SRC)
	      @echo "####### Generating dictionary"
	      @rootcint -f obj/NtupleProducerDict.cxx -c -p $(CXXFLAGS) \
	      $(INCLUDE_FLAGS) -I. $(INC) include/LinkDef.h

	      @echo "####### Building library $(LIB)"
	      @$(CXX) $(SOFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDE_FLAGS) -I. $(SRC) \
	      obj/NtupleProducerDict.cxx $(LDLIBS) -o $(LIB) $(ROOTLIBS) 
	      
	      @echo  "####### Removing generated dictionary"
	      @rm -f obj/NtupleProducerDict.cxx obj/NtupleProducerDict.h
	      @rm -f *.o

$(EXE):	      $(LIB) src/NtupleProducer.cxx
	      @echo "####### Building object file for executable"
	      @$(CXX) -c $(CXXFLAGS) $(INCLUDE_FLAGS) src/NtupleProducer.cxx -o obj/NtupleProducer.o
	      @echo "####### Building executable"
	      @$(CXX) $(CXXFLAGS) obj/NtupleProducer.o -o $(EXE) $(LDLIBS) $(LIB)
	      
$(PILEUP):   src/Make_PU_Histograms.cxx
	     @$(CXX) $(CXXFLAGS) $(GLIBS) src/Make_PU_Histograms.cxx -o $(PILEUP)


clean:
	      @rm -f $(OBJS) $(EXE) obj/NtupleProducerDict.cxx obj/NtupleProducerDict.h $(LIB)
