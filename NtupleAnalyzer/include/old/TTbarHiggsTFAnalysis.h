#ifndef TTbarHiggsTFAnalysis_H
#define TTbarHiggsTFAnalysis_H

#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include "Event.h"
#include "Muon.h"
#include "Tau.h"
#include "Electron.h"
#include "Jet.h"
#include "Truth.h"
#include <TObject.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>

#include "HistoManager.h"
#include "Lepton.h"

class TTbarHiggsTFAnalysis
{

 public:

   TTbarHiggsTFAnalysis();
   TTbarHiggsTFAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName,
                TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax);

   ~TTbarHiggsTFAnalysis();

   void createHistograms();
   void writeHistograms();

   TChain *fChain;   //!pointer to the analyzed TTree or TChain

   std::vector<Event>	 *vEvent     = new std::vector<Event>();
   std::vector<Electron> *vElectron  = new std::vector<Electron>();
   std::vector<Muon>	 *vMuon      = new std::vector<Muon>();
   std::vector<Tau>	 *vTau       = new std::vector<Tau>();
   std::vector<Jet>	 *vJet       = new std::vector<Jet>();
   std::vector<Truth>	 *vTruth     = new std::vector<Truth>();

   std::vector<Lepton>   vLeptons;
   std::vector<Muon>	 vSelectedMuons;
   std::vector<Electron> vSelectedElectrons;
   std::vector<Tau>	 vSelectedTaus;
   std::vector<Lepton>   vSelectedLeptons;
   std::vector<Muon>	 vFakeMuons;     // inverted MVA cut
   std::vector<Electron> vFakeElectrons; // inverted MVA cut
   std::vector<Lepton>   vFakeLeptons;   // inverted MVA cut
   std::vector<Jet>	 vSelectedJets;

   virtual void     Init(TChain *tree);
   virtual void     Loop();

   void initializeOutputTree();
   void fillOutputTree();

   float GetSignedDPhi(float phi1,float phi2);
   float GetDeltaR(float eta1,float phi1,float eta2,float phi2);

   TTree* tOutput;
   Int_t mc_event;
   Float_t weight;
   Float_t mc_weight;
   Float_t weight_PV; // PU reweighting from PV distribution


 private:

   HistoManager * theHistoManager;

   TFile * outputfile;

   TString _sampleName;
   TString _outputFileName;

   bool _isdata;
   float _xsec;
   float _lumi;
   int   _nowe; // number of weighted events
   int   _nmax; // max number of events to process

   TFile * _file_PVreweighting;
   TH1F* _h_PV;
   
};

#endif
