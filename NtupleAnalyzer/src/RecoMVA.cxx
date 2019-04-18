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

#include "../include/RecoMVA.h"

//https://root.cern.ch/doc/v606/classTMVA_1_1Reader.html

TMVA::Reader* Book_Reco_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weightMoriond2017/Hj_csv_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("b_from_leptop_bdt.csv", &bJet_fromLepTop_csv_var );
    reader->AddVariable("b_from_hadtop_bdt.csv", &bJet_fromHadTop_csv_var );
    reader->AddVariable("hadTop_tlv_bdt.Pt()", &hadTop_pt_var );
    reader->AddVariable("w_from_hadtop_tlv_bdt.M()", &w_fromHadTop_mass_var );
    reader->AddVariable("hadTop_tlv_bdt.M()", &hadTop_mass_var );
    reader->AddVariable("lep_from_higgs_bdt.obj.pt()", &lep_fromHiggs_pT_var );
    reader->AddVariable("lep_from_leptop_bdt.obj.pt()", &lep_fromTop_pT_var );
    reader->AddVariable("(lep_from_leptop_bdt.obj.pt()-lep_from_higgs_bdt.obj.pt())/(lep_from_leptop_bdt.obj.pt()+lep_from_higgs_bdt.obj.pt())", &lep_pt_ratio_var );
    reader->AddVariable("dr_lepFromTop_bFromLepTop", &dr_lepFromTop_bFromLepTop_var );
    reader->AddVariable("dr_lepFromTop_bFromHadTop", &dr_lepFromTop_bFromHadTop_var );
    reader->AddVariable("dr_lepFromHiggs_bFromLepTop", &dr_lepFromHiggs_bFromLepTop_var );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

void Load_MVA()
{

    //std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
    std::stringstream tmpBuffer;
    std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());

    std::string NtupleAnalyzerMVAPath = std::string("/opt/sbg/scratch1/cms/TTH/");
    mva_Reco_bTight  = Book_Reco_MVAReader( NtupleAnalyzerMVAPath, "weightMoriond2017/weights_factorized_bTight.xml", "test");
    mva_Reco_bLoose  = Book_Reco_MVAReader( NtupleAnalyzerMVAPath, "weightMoriond2017/weights_factorized_bLoose.xml",   "test");

    std::cout.rdbuf(oldStdout);
    //std::cout << "Stdout now restored." << std::endl;

    return;
}

