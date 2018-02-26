# Makefile for BDT codes

# $< = name of first dependancy
# $@ = name of target ("main")
# -c ask the compiler not to link
# -o specify target
#

#Define variables
CC=g++ #Compiler used
DEBUG = -g
CFLAGS = -c  $(DEBUG) #Compilation options
LFLAGS = $(DEBUG) -L/home/ntonon/root/lib -lTMVA #Modify
# LFLAGS = $(DEBUG) -lTMVA
ROOTFLAGS = `root-config --glibs --cflags`

SRCS = $(wildcard *.cxx) #Source files are all files with .cxx extension
HDRS = $(wildcard *.h) #Header files are all files with .h extension
OBJS = $(SRCS:.cxx=.o) #Object files are all files with .o extension, which have same names as source files

THQ_AN = analysis_main.exe #Name of executable file

#.PHONY : $(wildcard *.o)  #Force to always recompile object


#Instructions
all: $(THQ_AN)

#Obtain executables from object files
$(THQ_AN): analysis_main.o Func_other.o tHq_analysis.o
	@echo "###################################""#"
	@echo "-- Creating executable ./$(THQ_AN) --"
	@$(CC) analysis_main.o Func_other.o tHq_analysis.o -o $@ $(ROOTFLAGS) $(LFLAGS)
	@echo "-- Done --"
	@echo "###################################""#"
	@echo ""


#Obtain objects from source and header files
%.o: %.cxx $(HDRS)
	@echo "-- Compiling --"
	@$(CC) $(ROOTFLAGS) $(CFLAGS) $< -o $@
	@#Example :  g++ `root-config --glibs --cflags` -c -g Func_other.cxx -o Func_other.o


#Erase all objects and executable
clean:
	@-rm -f *.o *.exe
