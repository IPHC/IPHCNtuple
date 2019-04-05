#include "../include/Helper.h"
#include <assert.h>
#include "TH2F.h"
#include "TFile.h"

# define M_PI           3.14159265358979323846

bool SortingJetPt(Jet j1, Jet j2)
{
    if( j1.pt > j2.pt ) return true;
    else return false;
}

bool SortingLeptonPt( Lepton l1, Lepton l2)
{
    if( l1.pt > l2.pt ) return true;
    else return false;
}

bool SortingLeptonConePt( Lepton l1, Lepton l2)
{
    if( l1.conept > l2.conept ) return true;
    else return false;
}

float DeltaRLeptonJet( Lepton l1, Jet j1)
{
    float dEta = l1.eta - j1.eta;
    float dPhi = l1.phi - j1.phi;
    while(dPhi >= M_PI) dPhi -= 2*M_PI;
    while(dPhi < -M_PI) dPhi += 2*M_PI;
    float deltaR = sqrt( dEta*dEta + dPhi*dPhi);
    return deltaR;
}

float DeltaRJets( Jet j1, Jet j2)
{
    float dEta = j1.eta - j2.eta;
    float dPhi = j1.phi - j2.phi;
    while(dPhi >= M_PI) dPhi -= 2*M_PI;
    while(dPhi < -M_PI) dPhi += 2*M_PI;
    float deltaR = sqrt( dEta*dEta + dPhi*dPhi);
    return deltaR;
}

//Taken from : https://github.com/pallabidas/cmgtools-lite/blob/94X_dev_tHq_options/TTHAnalysis/python/plotter/tHq-multilepton/functionsTHQ.cc
//Fwd jet 2017 SFs derived in CR, to get alternate shape corresponding to fwd jet syst (tHq analysis)
float fwdjet_eventWeight_2017_option3_modified(float eta){
/*
Return an event weight based on the data/MC ratio of the maxJetEta25_60
distribution in OS emu events.
All jet pt cut 25 GeV, except if 2.7 < ans(eta) < 3.0, pt > 60.
*/
  eta = fabs(eta);
  if(eta < 0.278) return 0.9826;
  if(eta < 0.556) return 0.9742;
  if(eta < 0.833) return 0.9744;
  if(eta < 1.111) return 0.9685;
  if(eta < 1.389) return 1.0055;
  if(eta < 1.667) return 1.0185;
  if(eta < 1.944) return 0.9973;
  if(eta < 2.222) return 1.0175;
  if(eta < 2.500) return 1.0065;
  if(eta < 2.778) return 1.0715;
  if(eta < 3.056) return 1.2978;
  if(eta < 3.333) return 0.9388;
  if(eta < 3.611) return 1.0956;
  if(eta < 3.889) return 1.0444;
  if(eta < 4.167) return 1.0253;
  if(eta < 4.444) return 0.8298;
  if(eta < 4.722) return 0.5916;
  if(eta < 5.000) return 1.0000;
  return 1.0;
}
