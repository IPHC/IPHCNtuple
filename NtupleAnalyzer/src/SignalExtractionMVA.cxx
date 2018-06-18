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

#include "../include/SignalExtractionMVA.h"

//https://root.cern.ch/doc/v606/classTMVA_1_1Reader.html

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}

TMVA::Reader* Book_3L_THQ_MVAReader(std::string basePath, std::string weightFileName)
{
    if(!Check_File_Existence(basePath+"/"+weightFileName) )
    {
      cout<<BOLD(FRED("Weight file"<<basePath<<"/"<<weightFileName<<" not found !"))<<endl;
      return 0;
    }

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("nJet25_Recl", &nJet25);
    reader->AddVariable("nJetEta1", &nJetEta1);
    reader->AddVariable("maxEtaJet25", &maxEtaJet25);
    reader->AddVariable("dEtaFwdJetBJet", &dEtaFwdJetBJet);
    reader->AddVariable("dEtaFwdJetClosestLep", &dEtaFwdJetClosestLep);
    reader->AddVariable("dPhiHighestPtSSPair", &dPhiHighestPtSSPair);
    reader->AddVariable("LepGood_conePt[iLepFO_Recl[2]]", &Lep3Pt);
    reader->AddVariable("minDRll", &minDRll);
    reader->AddVariable("LepGood_charge[iLepFO_Recl[0]]+LepGood_charge[iLepFO_Recl[1]]+LepGood_charge[iLepFO_Recl[2]]", &lepCharge);
    reader->AddVariable("dEtaFwdJet2BJet", &dEtaFwdJet2BJet);

    //--- Needed from tHq2016 weight files
    Float_t dummy=0;
    reader->AddSpectator("iLepFO_Recl[0]", &dummy);
    reader->AddSpectator("iLepFO_Recl[1]", &dummy);
    reader->AddSpectator("iLepFO_Recl[2]", &dummy);


    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}


TMVA::Reader* Book_2L_THQ_MVAReader(std::string basePath, std::string weightFileName)
{
    if(!Check_File_Existence(basePath+"/"+weightFileName) )
    {
      cout<<BOLD(FRED("Weight file"<<basePath<<"/"<<weightFileName<<" not found !"))<<endl;
      return 0;
    }

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("nJet25_Recl", &nJet25);
    reader->AddVariable("nJetEta1", &nJetEta1);
    reader->AddVariable("maxEtaJet25", &maxEtaJet25);
    reader->AddVariable("dEtaFwdJetBJet", &dEtaFwdJetBJet);
    reader->AddVariable("dEtaFwdJetClosestLep", &dEtaFwdJetClosestLep);
    reader->AddVariable("dPhiHighestPtSSPair", &dPhiHighestPtSSPair);
    reader->AddVariable("LepGood_conePt[iLepFO_Recl[1]]", &Lep3Pt);
    reader->AddVariable("minDRll", &minDRll);
    reader->AddVariable("LepGood_charge[iLepFO_Recl[0]]+LepGood_charge[iLepFO_Recl[1]]", &lepCharge);
    reader->AddVariable("dEtaFwdJet2BJet", &dEtaFwdJet2BJet);

    //--- Needed from tHq2016 weight files
    // Float_t dummy=0;
    reader->AddSpectator("iLepFO_Recl[0]", (Float_t*) 0);
    reader->AddSpectator("iLepFO_Recl[1]", (Float_t*) 0);
    reader->AddSpectator("iLepFO_Recl[2]", (Float_t*) 0);

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}


void Load_MVA()
{
    //std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
    std::stringstream tmpBuffer;
    // std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());

    std::string NtupleAnalyzerMVAPath = std::string("/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer");
    // mva_3l_tt  = Book_3L_THQ_MVAReader(NtupleAnalyzerMVAPath, "test/weights_2016/thq_vs_tt_3l_BDTG.weights.xml");
	// mva_3l_ttV = Book_3L_THQ_MVAReader(NtupleAnalyzerMVAPath, "test/weights_2016/thq_vs_ttv_3l_BDTG.weights.xml");

	// mva_2lss_tt = Book_2L_THQ_MVAReader(NtupleAnalyzerMVAPath, "test/weights_2016/thq_vs_tt_2lss_BDTG.weights.xml.weights.xml");
	// mva_2lss_ttV = Book_2L_THQ_MVAReader(NtupleAnalyzerMVAPath, "test/weights_2016/thq_vs_ttv_2lss_BDTG.weights.xml.weights.xml");

    // std::cout.rdbuf(oldStdout);
    //std::cout << "Stdout now restored." << std::endl;

    return;
}
