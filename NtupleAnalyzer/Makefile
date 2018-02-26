# Makefile

#CC = g++ -g -Wno-writable-strings -Wno-write-strings -fsanitize=address -lTMVA
CC = g++ -g -Wno-writable-strings -Wno-write-strings -lTMVA

NTUPLEINC    = -I../NtupleProducer/include -I$(CMSSW_RELEASE_BASE)/src
NTUPLELIB    = -L../NtupleProducer -L$(CMSSW_RELEASE_BASE)/lib/slc6_amd64_gcc530 -lCondFormatsJetMETObjects \
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

LDLIBS  = $(ROOTLIBS) $(NTUPLELIB) $(KINFITLIB) 

ADDINC = $(NTUPLEINC) $(KINFITINC) 

SRC = src/TTbarHiggsMultileptonAnalysis.cxx \
 src/Analyzer.cxx \
 src/HistoManager.cxx src/Hist.cxx \
 src/Lepton.cxx \
 src/TTbarHiggsBTagEff.cxx \
 src/BTagCalibrationXStandalone.cxx \
 src/TTbarHiggsTFAnalysis.cxx


INC = include/TTbarHiggsMultileptonAnalysis.h \
 include/Analyzer.h \
 include/HistoManager.h include/Hist.h \
 include/Lepton.h \
 include/TTbarHiggsBTagEff.h \
 include/BTagCalibrationXStandalone.h \
 include/Ranges.h \
 include/TTbarHiggsTFAnalysis.h 

OBJ = TTbarHiggsMultileptonAnalysis.o \
 Analyzer.o\
 HistoManager.o Hist.o\
 Lepton.o \
 TTbarHiggsBTagEff.o\
 BTagCalibrationXStandalone.o \
 TTbarHiggsTFAnalysis.o

EXE = Analyzer

all:  	     		$(EXE)
#@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
