#ifndef Tau_H
#define Tau_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Tau : public Base
{
 public:
   
   Tau();
   virtual ~Tau();
   
   static bool sortPtPredicate(Tau lhs, Tau rhs)   {return (lhs.pt() > rhs.pt());};
   
   int ID()                                        {return _ID;};
   
   void setFakeType(int faketype)                  {_fakeType = faketype;};
   int  fakeType()                                 {return _fakeType;};
   
   void sel();
   void read();
   void init();
   
   float E()               {return _E;};
   float pt()              {return _pt;};
   float ptCor()           {return _ptCor;};
   float ptUnc()           {return _ptUnc;};
   float eta()             {return _eta;};
   float phi()             {return _phi;};
   float m()               {return _m;};
   int   charge()          {return _charge;};
   int   id()              {return _id;};
   
   bool isFakeableTTH()    {return _isFakeableTTH;};
   float lepMVA_TTH()      {return _lepMVA_TTH;};
   bool passTightCharge()  {return _passTightCharge;};
   bool cutEventSel()      {return _cutEventSel;};
   bool noLostHits()       {return _noLostHits;};
   
   float dxy()             {return _dxy;};
   float dz()              {return _dz;};
   
   float decayModeFinding()                    {return _decayModeFinding;};
//   float decayModeFindingOldDMs()              {return _decayModeFindingOldDMs;};
   float decayModeFindingNewDMs()              {return _decayModeFindingNewDMs;};
   
   float byLooseCombinedIsolationDeltaBetaCorr3Hits()              {return _byLooseCombinedIsolationDeltaBetaCorr3Hits;};
   float byMediumCombinedIsolationDeltaBetaCorr3Hits()             {return _byMediumCombinedIsolationDeltaBetaCorr3Hits;};
   float byTightCombinedIsolationDeltaBetaCorr3Hits()              {return _byTightCombinedIsolationDeltaBetaCorr3Hits;};
   float againstElectronVLooseMVA6()                               {return _againstElectronVLooseMVA6;};
   float againstElectronLooseMVA6()                                {return _againstElectronLooseMVA6;};
   float againstElectronMediumMVA6()                               {return _againstElectronMediumMVA6;};
   float againstElectronTightMVA6()                                {return _againstElectronTightMVA6;};
   
   float byVLooseIsolationMVArun2v1DBdR03oldDMwLT()                {return _byVLooseIsolationMVArun2v1DBdR03oldDMwLT;};
   float byLooseIsolationMVArun2v1DBdR03oldDMwLT()                 {return _byLooseIsolationMVArun2v1DBdR03oldDMwLT;};
   float byMediumIsolationMVArun2v1DBdR03oldDMwLT()                {return _byMediumIsolationMVArun2v1DBdR03oldDMwLT;};
   float byTightIsolationMVArun2v1DBdR03oldDMwLT()                 {return _byTightIsolationMVArun2v1DBdR03oldDMwLT;};
   float byVTightIsolationMVArun2v1DBdR03oldDMwLT()                {return _byVTightIsolationMVArun2v1DBdR03oldDMwLT;};
   
   float byCombinedIsolationDeltaBetaCorrRaw3Hits()                {return _byCombinedIsolationDeltaBetaCorrRaw3Hits;};
   float againstMuonLoose3()                                       {return _againstMuonLoose3;};
   float againstMuonTight3()                                       {return _againstMuonTight3;};

   bool hasMCMatch() {return _hasMCMatch;};
   float gen_pt() {return _gen_pt;};
   float gen_eta() {return _gen_eta;};
   float gen_phi() {return _gen_phi;};
   float gen_m() {return _gen_m;};
   float gen_E() {return _gen_E;};
   int gen_status() {return _gen_status;};
   int gen_id() {return _gen_id;};
   int gen_charge() {return _gen_charge;};
   float gen_dr() {return _gen_dr;};
   
 protected:
   
   int   _ID;
   
   int   _fakeType;
   
   float _E;
   float _pt;
   float _ptCor;
   float _ptUnc;
   float _eta;
   float _phi;
   float _m;
   int   _charge;
   int   _id;
   
   float _dxy;
   float _dz;
   
   bool  _isFakeableTTH;
   float _lepMVA_TTH;
   bool _passTightCharge;
   bool _cutEventSel;
   bool _noLostHits;

   int   _decayMode;
   bool  _hasLeadChargedHadrCand;
   float _leadingTrackPt;
   float _leadingTrackCharge;

   float _chargedIsoPtSum;
   float _neutralIsoPtSum;
   float _puCorrPtSum;
   
   float _pfEssential_jet_pt;
   float _pfEssential_jet_eta;
   float _pfEssential_jet_phi;
   float _pfEssential_jet_m;
   float _pfEssential_jetCorr_pt;
   float _pfEssential_jetCorr_eta;
   float _pfEssential_jetCorr_phi;
   float _pfEssential_jetCorr_m;
   float _pfEssential_hasSV;
   float _pfEssential_sv_x;
   float _pfEssential_sv_y;
   float _pfEssential_sv_z;
   float _pfEssential_flightLengthSig;
   float _pfEssential_dxy;
   float _pfEssential_dxy_error;
   float _pfEssential_dxy_Sig;

   float _decayModeFinding;
//   float _decayModeFindingOldDMs;
   float _decayModeFindingNewDMs;
   float _byLooseCombinedIsolationDeltaBetaCorr3Hits;
   float _byMediumCombinedIsolationDeltaBetaCorr3Hits;
   float _byTightCombinedIsolationDeltaBetaCorr3Hits;
   float _againstElectronVLooseMVA6;
   float _againstElectronLooseMVA6;
   float _againstElectronMediumMVA6;
   float _againstElectronTightMVA6;

   float _byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float _byLooseIsolationMVArun2v1DBdR03oldDMwLT;
   float _byMediumIsolationMVArun2v1DBdR03oldDMwLT;
   float _byTightIsolationMVArun2v1DBdR03oldDMwLT;
   float _byVTightIsolationMVArun2v1DBdR03oldDMwLT;

   float _byCombinedIsolationDeltaBetaCorrRaw3Hits;
   float _againstMuonLoose3;
   float _againstMuonTight3;

   bool _hasMCMatch;
   float _gen_pt;
   float _gen_eta;
   float _gen_phi;
   float _gen_m;
   float _gen_E;
   int _gen_status;
   int _gen_id;
   int _gen_charge;
   float _gen_dr;
   
   ClassDef(Tau,1)
};

#endif
