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
//-------------- OLD
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
// float signal_2lss_TT_MVA;
// float signal_2lss_TTV_MVA;
float signal_3l_TT_MVA;
float signal_3l_TTV_MVA;

// TMVA readers
// TMVA::Reader* mva_2lss_tt;
// TMVA::Reader* mva_2lss_ttV;
TMVA::Reader* mva_3l_tt;
TMVA::Reader* mva_3l_ttV;

// Function definitions
// TMVA::Reader* Book_2LSS_TT_MVAReader(  std::string basePath, std::string weightFileName, std::string type);
// TMVA::Reader* Book_2LSS_TTV_MVAReader( std::string basePath, std::string weightFileName, std::string type);
// TMVA::Reader* Book_3L_TT_MVAReader(    std::string basePath, std::string weightFileName, std::string type);
// TMVA::Reader* Book_3L_TTV_MVAReader(   std::string basePath, std::string weightFileName, std::string type);

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
TMVA::Reader* Book_3L_THQ_MVAReader(  std::string basePath, std::string weightFileName);
void Load_MVA();



#endif
