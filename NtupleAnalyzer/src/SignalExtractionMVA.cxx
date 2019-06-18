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

TMVA::Reader* Book_THQ2016_MVAReader(std::string basePath, std::string weightFileName, int nLep)
{
    if(!Check_File_Existence(basePath+"/"+weightFileName) )
    {
      cout<<BOLD(FRED("Weight file"<<basePath<<"/"<<weightFileName<<" not found !"))<<endl;
      return 0;
    }

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    if(nLep == 2) //different def of Lep3Pt & lepCharge
    {
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
    }
    else if(nLep == 3)
    {
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
    }
    else  {cout<<"Error ! Wrong nof leptons (2 or 3)"<<endl; return 0;}

    //--- Needed from tHq2016 weight files
    Float_t dummy=0;
    reader->AddSpectator("iLepFO_Recl[0]", &dummy);
    reader->AddSpectator("iLepFO_Recl[1]", &dummy);
    reader->AddSpectator("iLepFO_Recl[2]", &dummy);


    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

//2017 training, added 'fwdJetPt25' variable
TMVA::Reader* Book_THQ2017_MVAReader(std::string basePath, std::string weightFileName, int nLep)
{
  if(!Check_File_Existence(basePath+"/"+weightFileName) )
  {
    cout<<BOLD(FRED("Weight file"<<basePath<<"/"<<weightFileName<<" not found !"))<<endl;
    return 0;
  }

  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

  if(nLep == 2) //different def of Lep3Pt & lepCharge
  {
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
      reader->AddVariable("fwdJetPt25", &FwdJetPt);
      reader->AddVariable("BDThttTT_eventReco_mvaValue", &resHTT);
      reader->AddVariable("sumJetsPt", &sum_jetPt);
  }
  else if(nLep == 3)
  {
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
      reader->AddVariable("fwdJetPt25", &FwdJetPt);
  }
  else  {cout<<"Error ! Wrong nof leptons (2 or 3)"<<endl; return 0;}

  //--- Needed from tHq2016 weight files
  Float_t dummy=0;
  reader->AddSpectator("iLepFO_Recl[0]", &dummy);
  reader->AddSpectator("iLepFO_Recl[1]", &dummy);
  reader->AddSpectator("iLepFO_Recl[2]", &dummy);

  reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

  return reader;
}

TMVA::Reader* Book_TTH_MVAReader(std::string basePath, std::string weightFileName, std::string BDTtype, int nLep)
{
    if(!Check_File_Existence(basePath+"/"+weightFileName) )
    {
      cout<<BOLD(FRED("Weight file "<<basePath<<"/"<<weightFileName<<" not found !"))<<endl;
      return 0;
    }

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    if(BDTtype == "ttV")
    {
      reader->AddVariable("max(abs(LepGood_eta[iLepFO_Recl[0]]),abs(LepGood_eta[iLepFO_Recl[1]]))", &max_lep_eta);
      reader->AddVariable("nJet25_Recl", &nJet25);
      reader->AddVariable("mindr_lep1_jet", &mindr_lep1_jet);
      reader->AddVariable("mindr_lep2_jet", &mindr_lep2_jet);
      reader->AddVariable("MT_met_lep1", &mT_lep1);
      if(nLep == 2) {reader->AddVariable("LepGood_conePt[iLepFO_Recl[1]]", &lep3_conePt);}
      else {reader->AddVariable("LepGood_conePt[iLepFO_Recl[2]]", &lep3_conePt);}
      reader->AddVariable("LepGood_conePt[iLepFO_Recl[0]]", &lep1_conePt);

    }
    else if(BDTtype == "ttbar")
    {
      reader->AddVariable("max(abs(LepGood_eta[iLepFO_Recl[0]]),abs(LepGood_eta[iLepFO_Recl[1]]))", &max_lep_eta);
      reader->AddVariable("nJet25_Recl", &nJet25);
      reader->AddVariable("mindr_lep1_jet", &mindr_lep1_jet);
      reader->AddVariable("mindr_lep2_jet", &mindr_lep2_jet);
      reader->AddVariable("MT_met_lep1", &mT_lep1);
    }
    else {cout<<"ERROR ! Wrong BDT type name ! (SignalExtraction.cxx)"<<endl;}

    //--- Needed from tHq2016 weight files
    Float_t dummy=0;
    reader->AddSpectator("iLepFO_Recl[0]", &dummy);
    reader->AddSpectator("iLepFO_Recl[1]", &dummy);
    reader->AddSpectator("iLepFO_Recl[2]", &dummy);

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    std::cout<<"--> Booked MVA : "<<basePath<<"/"<<weightFileName<<std::endl;

    return reader;
}


TMVA::Reader* Book_HJTagger_MVAReader(std::string weightPath)
{
  if(!Check_File_Existence(weightPath) )
  {
    cout<<BOLD(FRED("Weight file"<<weightPath<<" not found !"))<<endl;
    return 0;
  }

  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

  reader->AddVariable("Jet25_lepdrmin", &Jet25_lepdrmin);
  reader->AddVariable("max(Jet25_bDiscriminator,0.)", &Jet25_bDiscriminator);
  reader->AddVariable("max(Jet25_qg,0.)", &Jet25_qg);
  reader->AddVariable("Jet25_lepdrmax", &Jet25_lepdrmax);
  reader->AddVariable("Jet25_pt", &Jet25_pt);

  reader->BookMVA("BDTG method", weightPath);

  return reader;
}


TMVA::Reader* Book_resHTT_MVAReader(std::string weightPath)
{
  if(!Check_File_Existence(weightPath) )
  {
    cout<<BOLD(FRED("Weight file"<<weightPath<<" not found !"))<<endl;
    return 0;
  }

  TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

  reader->AddVariable("CSV_b", &DeepCSV_bjet_kinFit);
  reader->AddVariable("qg_Wj2", &qg_Wj2);
  reader->AddVariable("pT_bWj1Wj2", &pT_bWj1Wj2);
  reader->AddVariable("m_Wj1Wj2", &m_Wj1Wj2);
  reader->AddVariable("nllKinFit", &nllKinFit);
  reader->AddVariable("pT_b_o_kinFit_pT_b", &pT_b_o_kinFit_pT_b);
  reader->AddVariable("pT_Wj2", &pT_Wj2);

  reader->BookMVA("BDT method", weightPath);

  return reader;
}



void Load_MVA(TString analysis_type)
{
    // std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
    // std::stringstream tmpBuffer;
    // std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());

    if(analysis_type == "tHq")
    {
      std::string NtupleAnalyzerMVAPath = std::string("/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/BDT_Pallabi");
      mva_3l_tt  = Book_THQ2017_MVAReader(NtupleAnalyzerMVAPath, "thq_vs_tt_3l_BDTG.weights.xml", 3);
      mva_3l_ttV = Book_THQ2017_MVAReader(NtupleAnalyzerMVAPath, "thq_vs_ttv_3l_BDTG.weights.xml", 3);
      mva_2lss_tt  = Book_THQ2017_MVAReader(NtupleAnalyzerMVAPath, "thq_vs_tt_2lss_BDTG.weights.xml", 2);
      mva_2lss_ttV = Book_THQ2017_MVAReader(NtupleAnalyzerMVAPath, "thq_vs_ttv_2lss_BDTG.weights.xml", 2);
    }
    else if(analysis_type == "ttH")
    {
      std::string NtupleAnalyzerMVAPath = std::string("/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/BDT_ttH2017");
      mva_3l_tt  = Book_TTH_MVAReader(NtupleAnalyzerMVAPath, "3l_ttbar_BDTG.weights.xml", "ttbar", 3);
      mva_3l_ttV = Book_TTH_MVAReader(NtupleAnalyzerMVAPath, "3l_ttV_BDTG.weights.xml", "ttV", 3);
      mva_2lss_tt  = Book_TTH_MVAReader(NtupleAnalyzerMVAPath, "2lss_ttbar_BDTG.weights.xml", "ttbar", 2);
      mva_2lss_ttV = Book_TTH_MVAReader(NtupleAnalyzerMVAPath, "2lss_ttV_BDTG.weights.xml", "ttV", 2);
    }
    else if(analysis_type != "FCNC") {cout<<"Error ! Wrong analysis_type value : "<<analysis_type<<endl;}

    // std::string HJTagger_weight_path = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Taggers/Hj_deepcsv_BDTG_2017.weights.xml";
    std::string HJTagger_weight_path = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Taggers/Hj_2017_configA_dcsv_BDTG.weights.xml";
    mva_HjTagger = Book_HJTagger_MVAReader(HJTagger_weight_path);

    std::string resHTT_weight_path = "/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleAnalyzer/test/weights_2017/Taggers/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml";
    mva_resHTT = Book_resHTT_MVAReader(resHTT_weight_path);

    // std::cout.rdbuf(oldStdout);
    // std::cout << "Stdout now restored." << std::endl;

    return;
}
