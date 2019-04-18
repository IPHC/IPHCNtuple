#ifndef MVA
#define MVA

#include <memory>
#include <iostream>
#include <sstream>

#include "TRegexp.h"
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TFile.h"
#include "TThreadSlots.h"
#include "TROOT.h"
#include "Compression.h"

// Input variables
float max_Lep_eta;
float mindr_lep1_jet;
float mindr_lep2_jet;
float met;
float avg_dr_jet;
float MT_met_lep1;
float LepGood_conePt0;
float LepGood_conePt1;
float nJet25_Recl;
float mhtJet25_Recl;

// Input variables (spectators)
float iF_Recl0;
float iF_Recl1;
float iF_Recl2;

// Output variables
float signal_Hj_MVA;
float signal_Hjj_MVA;

// TMVA readers
TMVA::Reader* mva_Hj;
TMVA::Reader* mva_Hjj;

// Function definitions
TMVA::Reader* Book_Hj_MVAReader(  std::string basePath, std::string weightFileName, std::string type);
TMVA::Reader* Book_Hjj_MVAReader( std::string basePath, std::string weightFileName, std::string type);
void Load_MVA();

#endif
