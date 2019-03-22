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

#include "tHqMultileponAnalysis.h"

using namespace std;

//----------------------
// Input variables (spectators)
float iF_Recl0;
float iF_Recl1;
float iF_Recl2;

// Output variables
float signal_2lss_TT_MVA;
float signal_2lss_TTV_MVA;
float signal_3l_TT_MVA;
float signal_3l_TTV_MVA;
// float signal_TT_MVA;
// float signal_TTV_MVA;

//--- DEFINE HERE VARIABLES THAT MAY BE NEEDED TO COMPUTE BDT OUTPUT AT NTA LEVEL (i.e. variables used in ttH2017 or tHq2016 BDTs, Hj tagger, ...)
//Because else problem with variable scope

//tHq vars
Float_t nJet25;
Float_t nJetEta1 ;
Float_t maxEtaJet25;
Float_t dEtaFwdJetBJet;
Float_t dEtaFwdJetClosestLep;
Float_t dPhiHighestPtSSPair;
Float_t Lep3Pt;
Float_t minDRll;
Float_t lepCharge;
Float_t dEtaFwdJet2BJet;
Float_t FwdJetPt; //For tHq2017


//ttH2017 vars
Float_t lep1_conePt;
Float_t lep2_conePt;
Float_t lep3_conePt;
Float_t mindr_lep1_jet;
Float_t mindr_lep2_jet;
Float_t mT_lep1;
Float_t mT_lep2;
Float_t max_lep_eta;

//Needed for Hj tagger
Float_t Jet25_lepdrmin; //Min dr between given jet and any lepton
Float_t Jet25_lepdrmax; //Max...
Float_t Jet25_qg; //quark-gluon discriminant for jet
Float_t Jet25_bDiscriminator;
Float_t Jet25_pt;

//Needed for resHTT (resolved hadronic top tagger)
Float_t DeepCSV_bjet_kinFit;
Float_t qg_Wj2;
Float_t pT_bWj1Wj2;
Float_t pT_Wj2;
Float_t m_Wj1Wj2;
Float_t nllKinFit;
Float_t pT_b_o_kinFit_pT_b;

// TMVA readers
TMVA::Reader* mva_2lss_tt;
TMVA::Reader* mva_2lss_ttV;
TMVA::Reader* mva_3l_tt;
TMVA::Reader* mva_3l_ttV;

TMVA::Reader* mva_HjTagger;
TMVA::Reader* mva_resHTT;

// Function definitions
bool Check_File_Existence(const TString&);

TMVA::Reader* Book_TTH_MVAReader(std::string basePath, std::string weightFileName, std::string BDTtype, int nLep);

TMVA::Reader* Book_HJTagger_MVAReader(std::string);

TMVA::Reader* Book_resHTT_MVAReader(std::string);

TMVA::Reader* Book_THQ2016_MVAReader(std::string basePath, std::string weightFileName, int nLep);

TMVA::Reader* Book_THQ2017_MVAReader(std::string basePath, std::string weightFileName, int nLep);

void Load_MVA(std::string analysis_type);

#endif
