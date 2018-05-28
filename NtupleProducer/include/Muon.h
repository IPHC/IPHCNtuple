#ifndef MUON_H
#define MUON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Muon : public Base
{
 public:
   
   Muon();
   virtual ~Muon();
   
   static bool sortPtPredicate(Muon lhs, Muon rhs)  {return (lhs.pt() > rhs.pt());};
   
   int ID()                                         {return _ID;};
   
   void setFakeType(int faketype)                   {_fakeType = faketype;};
   int fakeType()                                   {return _fakeType;};
   
   void sel();
   void read();
   void init();
   
   float E()                  {return _E;};
   float pt()                 {return _pt;};
   float ptCor()              {return _ptCor;};
   float ptUnc()              {return _ptUnc;};
   float eta()                {return _eta;};
   float phi()                {return _phi;};
   float m()                  {return _m;};
   float conept()             {return _conept;};
   int charge()               {return _charge;};
   int id()                   {return _id;};
   
   bool isLoose()             {return _isLoose;};
   bool isMedium()            {return _isMedium;};
   bool isTight()             {return _isTight;};
   bool isPFMuon()            {return _isPFMuon;};
   
   bool isLooseTTH()          {return _isLooseTTH;};
   bool isFakeableTTH()       {return _isFakeableTTH;};
   bool isTightTTH()          {return _isTightTTH;};
   
   float dxy()                {return _dxy;};
   float dz()                 {return _dz;};
   float iso()                {return _iso;};
   float isoR04()             {return _isoR04;};
   float sip3d()              {return _sip3d;};
   float bestTrackpt()        {return _bestTrack_pt;};
   float bestTrackptError()   {return _bestTrack_ptError;};
   bool  cutEventSel()        {return _cutEventSel;};
   bool  noLostHits()         {return _noLostHits;};
   
   float lepMVA()             {return _lepMVA;};
   
   float lepMVA_miniRelIsoCharged()    {return _lepMVA_miniRelIsoCharged;};
   float lepMVA_miniRelIsoNeutral()    {return _lepMVA_miniRelIsoNeutral;};
   float lepMVA_jetPtRelv2()           {return _lepMVA_jetPtRelv2;};
   float lepMVA_jetPtRatio()    	    {return _lepMVA_jetPtRatio;};
   float lepMVA_jetBTagCSV()    	    {return _lepMVA_jetBTagCSV;};
   float lepMVA_jetBTagDeepCSV()    	    {return _lepMVA_jetBTagDeepCSV;};
   float lepMVA_sip3d()         	    {return _lepMVA_sip3d;};
   float lepMVA_dxy()           	    {return _lepMVA_dxy;};
   float lepMVA_dz()            	    {return _lepMVA_dz;};
   float lepMVA_mvaId()         	    {return _lepMVA_mvaId;};
   float lepMVA_eta()                  {return _lepMVA_eta;};
   float lepMVA_jetNDauChargedMVASel() {return _lepMVA_jetNDauChargedMVASel;};

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
   
   int _ID;
   
   int _fakeType;
   
//   float _E;
//   float _pt;
   float _ptCor;
   float _ptUnc;
//   float _eta;
//   float _phi;
   float _m;
   float _conept;
   int _id;
   
   bool _isLoose;
   bool _isMedium;
   bool _isTight;
   bool _isPFMuon;
   
   bool _isLooseTTH;
   bool _isFakeableTTH;
   
   float _dxy;
   float _dz;
   float _iso;
   float _isoR04;
   float _sip3d;
   float _bestTrack_pt;
   float _bestTrack_ptError;
   bool _cutEventSel;
   bool _noLostHits;

   float _lepMVA_miniRelIsoCharged;
   float _lepMVA_miniRelIsoNeutral;
   float _lepMVA_jetPtRelv2;
   float _lepMVA_jetPtRatio;
   float _lepMVA_jetBTagCSV;
   float _lepMVA_jetBTagDeepCSV;
   float _lepMVA_sip3d;
   float _lepMVA_dxy;
   float _lepMVA_dz;
   float _lepMVA_mvaId;
   float _lepMVA_eta;
   float _lepMVA_jetNDauChargedMVASel;

   float _gen_pt;
   float _gen_eta;
   float _gen_phi;
   float _gen_m;
   float _gen_E;
   int _gen_status;
   int _gen_id;
   int _gen_charge;
   float _gen_dr;
   
   ClassDef(Muon,1)
};

#endif
