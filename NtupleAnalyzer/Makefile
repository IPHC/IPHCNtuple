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

all:  	     		$(EXE)
#@rm -f *.o

$(EXE):     	 	$(OBJ) $@

$(OBJ):     		$(SRC) $(INC)
		        $(CC) -c $(ROOTCFLAGS) $(ADDINC) $(SRC)

clean: 		
			@rm -f $(OBJ) $(EXE)
