# Makefile

#CC = g++ -g -Wno-writable-strings -Wno-write-strings -fsanitize=address -lTMVA
CC = g++ -g -lTMVA

# gsl -- needed for kinFit, resHTT
GSLINCS = $(shell gsl-config --cflags)
GSLLIBS = $(shell gsl-config --libs)

BOOSTSYS       = $(shell scram tool tag boost BOOST_BASE)
BOOSTINC       = -I$(BOOSTSYS)/include

#NB : first, must compile NtupleProducer, to produce libNtuple.so !
#Added '$(GSLLIBS) -lGenVector' for resHTT
NTUPLEINC    = -I../NtupleProducer/include -I$(CMSSW_RELEASE_BASE)/src $(GSLINCS)
NTUPLELIB    = -L../NtupleProducer -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc630 -lCondFormatsJetMETObjects -lNtuple $(GSLLIBS) -lGenVector

SOVERSION               = $(shell echo $(VERSION) | \
                               sed 's/\(.*\)\.\(.*\)\.\(.*\)/\1:\2:0/')
ROOTCONF                = root-config
ROOTCINT                = $(shell $(ROOTCONF) --prefix)/bin/rootcint
ROOTLIBS                = $(shell $(ROOTCONF) --libs)
ROOTGLIBS               = $(shell $(ROOTCONF) --glibs)
ROOTCFLAGS              = $(shell $(ROOTCONF) --cflags)
INCLUDES                = $(shell $(ROOTCONF) --incdir) $(ROOTCFLAGS)
CINTINCLUDES            = $(shell $(ROOTCONF) --incdir) 
LIBS                    = $(ROOTLIBS)
ROOTFLAGS2 		= `root-config --glibs --cflags`


LDLIBS  = $(ROOTLIBS) $(NTUPLELIB) $(KINFITLIB) 

ADDINC = $(NTUPLEINC) $(KINFITINC) $(BOOSTINC)

SRC = src/tHqMultileponAnalysis.cxx \
 src/Analyzer.cxx \
 src/Lepton.cxx \
 src/BTagCalibrationXStandalone.cxx \
 src/ScaleFactors.cxx \
 src/Pileup.cxx \

INC = include/tHqMultileponAnalysis.h \
 include/Analyzer.h \
 include/Lepton.h \
 include/BTagCalibrationXStandalone.h \
 include/ScaleFactors.h \
 include/Pileup.h \

OBJ = tHqMultileponAnalysis.o \
 Analyzer.o\
 Lepton.o \
 BTagCalibrationXStandalone.o \
 ScaleFactors.o \
 Pileup.o \

EXE = Analyzer
COMPUTE_SUMWEIGHT = test/Get_SumWeightEvents_FlatTrees.exe
SAVE_SUMWEIGHT = test/Store_SumsOfWeights_inMergedNtuple.exe

all:  	     		$(EXE) $(COMPUTE_SUMWEIGHT) $(SAVE_SUMWEIGHT)

$(EXE):     	 	$(OBJ) test/HTT_kinfit/HadTopKinFit.o $@


$(OBJ):     		$(SRC) $(INC)
			@echo ""
			@echo "###################################"
			@echo "-- Creating object file $(OBJ) --"
		        @$(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)
			@echo "-- Done --"
			@echo "###################################"
			@echo ""

$(COMPUTE_SUMWEIGHT):	src/Get_SumWeightEvents_FlatTrees.cxx	
			@echo ""
			@echo "###################################"
			@echo "-- Creating executable $(SUMWEIGHT) --"
			@$(CC) $(ROOTFLAGS2) src/Get_SumWeightEvents_FlatTrees.cxx -o $@
			@echo "-- Done --"
			@echo "###################################"
			@echo ""


$(SAVE_SUMWEIGHT):	src/Store_SumsOfWeights_inMergedNtuple.cxx	
			@echo ""
			@echo "###################################"
			@echo "-- Creating executable $(SAVE_SUMWEIGHT) --"
			@$(CC) $(ROOTFLAGS2) src/Store_SumsOfWeights_inMergedNtuple.cxx -o $@
			@echo "-- Done --"
			@echo "###################################"
			@echo ""


clean: 		
			@rm -f $(OBJ) $(EXE) $(COMPUTE_SUMWEIGHT) $(SAVE_SUMWEIGHT)
