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
float bJet_fromLepTop_csv_var;
float bJet_fromHadTop_csv_var;
float hadTop_pt_var;
float w_fromHadTop_mass_var;
float hadTop_mass_var;
float lep_fromHiggs_pT_var;
float lep_fromTop_pT_var;
float lep_pt_ratio_var;
float dr_lepFromTop_bFromLepTop_var;
float dr_lepFromTop_bFromHadTop_var;
float dr_lepFromHiggs_bFromLepTop_var

// Output variables
float signal_Reco_bTight_MVA;
float signal_Reco_bLoose_MVA;

// TMVA readers
TMVA::Reader* mva_Reco_bTight;
TMVA::Reader* mva_Reco_bLoose;

// Function definitions
TMVA::Reader* Book_Reco_MVAReader(  std::string basePath, std::string weightFileName, std::string type);
void Load_MVA();

#endif
