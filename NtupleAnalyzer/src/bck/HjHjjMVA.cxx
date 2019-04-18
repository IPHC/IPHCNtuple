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

#include "../include/HjHjjMVA.h"

//https://root.cern.ch/doc/v606/classTMVA_1_1Reader.html

TMVA::Reader* Book_Hj_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weightMoriond2017/Hj_csv_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))",  &max_Lep_eta      );
    reader->AddVariable("nJet25_Recl",                                                     &nJet25_Recl      );
    reader->AddVariable("mindr_lep1_jet",                                                  &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                  &mindr_lep2_jet   );
    reader->AddVariable("min(met_pt,400)",                                                 &met              );
    reader->AddVariable("avg_dr_jet",                                                      &avg_dr_jet       );
    reader->AddVariable("MT_met_lep1",                                                     &MT_met_lep1      );

    reader->AddSpectator("iF_Recl[0]",                                                     &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                     &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                     &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

TMVA::Reader* Book_Hjj_MVAReader(std::string basePath, std::string weightFileName, std::string type)
{

    // weightMoriond2017/Hjj_csv_BDTG.weights.xml
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    reader->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_Lep_eta      );
    reader->AddVariable("MT_met_lep1",                                                    &MT_met_lep1      );
    reader->AddVariable("nJet25_Recl",                                                    &nJet25_Recl      );
    reader->AddVariable("mindr_lep1_jet",                                                 &mindr_lep1_jet   );
    reader->AddVariable("mindr_lep2_jet",                                                 &mindr_lep2_jet   );
    reader->AddVariable("LepGood_conePt[iF_Recl[0]]",                                     &LepGood_conePt0  );
    reader->AddVariable("LepGood_conePt[iF_Recl[1]]",                                     &LepGood_conePt1  );

    reader->AddSpectator("iF_Recl[0]",                                                    &iF_Recl0         );
    reader->AddSpectator("iF_Recl[1]",                                                    &iF_Recl1         );
    reader->AddSpectator("iF_Recl[2]",                                                    &iF_Recl2         );

    reader->BookMVA("BDTG method", basePath+"/"+weightFileName);

    return reader;
}

void Load_MVA()
{

    //std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
    std::stringstream tmpBuffer;
    std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());

    std::string NtupleAnalyzerMVAPath = std::string("/opt/sbg/scratch1/cms/TTH/");
    mva_Hj  = Book_Hj_MVAReader(  NtupleAnalyzerMVAPath, "weightMoriond2017/Hj_csv_BDTG.weights.xml", "test");
    mva_Hjj = Book_Hjj_MVAReader( NtupleAnalyzerMVAPath, "weightMoriond2017/Hj_csv_BDTG.weights.xml",   "test");

    std::cout.rdbuf(oldStdout);
    //std::cout << "Stdout now restored." << std::endl;

    return;
}

