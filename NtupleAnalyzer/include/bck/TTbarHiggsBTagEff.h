#ifndef TTbarHiggsBTagEff_H
#define TTbarHiggsBTagEff_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "TriggerObj.h"
#include "Event.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include "Base.h"

#include "HistoManager.h"
#include "Lepton.h"

class TTbarHiggsBTagEff
{

 public:
   TTbarHiggsBTagEff();
   TTbarHiggsBTagEff(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName,
                             //bool isdata, float xsec, float lumi, int nowe, 
			     int nmax);

   ~TTbarHiggsBTagEff();

   void createHistograms();
   void writeHistograms();
  
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   virtual void     Init(TChain *tree);

   virtual void     Loop();
   
   std::vector<Jet>        *vJet        = new std::vector<Jet>();

 
 private:

   HistoManager * theHistoManager;

   TFile * _outputFile;
   TString _outputFileName;

   TString _sampleName;
  
   //bool  _isdata;
   //float _xsec;
   //float _lumi;
   //int   _nowe; // number of weighted events
   int   _nmax; // max number of events to process
   
};

#endif
