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
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

//----------------------
// Input variables (spectators)
float iF_Recl0;
float iF_Recl1;
float iF_Recl2;

// Output variables
float signal_2l_TT_MVA;
float signal_2l_TTV_MVA;
float signal_3l_TT_MVA;
float signal_3l_TTV_MVA;

// TMVA readers
TMVA::Reader* mva_2lss_tt;
TMVA::Reader* mva_2lss_ttV;
TMVA::Reader* mva_3l_tt;
TMVA::Reader* mva_3l_ttV;

//----------------------
//-------------- New
Float_t nJet25;
Float_t maxEtaJet25;
Float_t lepCharge;
Float_t nJetEta1 ;
Float_t dEtaFwdJetBJet;
Float_t dEtaFwdJet2BJet;
Float_t dEtaFwdJetClosestLep;
Float_t dPhiHighestPtSSPair;
Float_t minDRll;
Float_t Lep3Pt;

// Function definitions
bool Check_File_Existence(const TString&);
TMVA::Reader* Book_3L_THQ_MVAReader(std::string basePath, std::string weightFileName);
TMVA::Reader* Book_2L_THQ_MVAReader(std::string basePath, std::string weightFileName);
void Load_MVA();



#endif
