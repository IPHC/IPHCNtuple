# Makefile

#CC = g++ -g -Wno-writable-strings -Wno-write-strings -fsanitize=address -lTMVA
CC = g++ -g -Wno-writable-strings -Wno-write-strings -lTMVA

NTUPLEINC    = -I../NtupleProducer/include -I$(CMSSW_RELEASE_BASE)/src
NTUPLELIB    = -L../NtupleProducer -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc630 -lCondFormatsJetMETObjects \
              -lNtuple

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

ADDINC = $(NTUPLEINC) $(KINFITINC) 

SRC = src/tHqMultileponAnalysis.cxx \
 src/Analyzer.cxx \
 src/Lepton.cxx \
 src/BTagCalibrationXStandalone.cxx \


INC = include/tHqMultileponAnalysis.h \
 include/Analyzer.h \
 include/Lepton.h \
 include/BTagCalibrationXStandalone.h \
 include/Ranges.h \

OBJ = tHqMultileponAnalysis.o \
 Analyzer.o\
 Lepton.o \
 BTagCalibrationXStandalone.o \

EXE = Analyzer

SUMWEIGHT = test/Get_SumWeightEvents_FlatTrees.exe


all:  	     		$(EXE) $(SUMWEIGHT)
#@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        @$(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)
			@echo ""
			@echo "###################################"
			@echo "-- Creating executable $(EXE) --"
			@echo "-- Done --"
			@echo "###################################"
			@echo ""

$(SUMWEIGHT):		
			@echo ""
			@echo "###################################"
			@echo "-- Creating executable $(SUMWEIGHT) --"
			@$(CC) $(ROOTFLAGS2) test/Get_SumWeightEvents_FlatTrees.cxx -o $@
			@echo "-- Done --"
			@echo "###################################"
			@echo ""

clean: 		
			@rm -f $(OBJ) $(EXE) $(SUMWEIGHT)
