#ifndef JET_H
#define JET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Jet : public Base
{
 public:
   
   Jet();
   virtual ~Jet();
   
   static bool sortPtPredicate(Jet lhs, Jet rhs) {return (lhs.pt() > rhs.pt());};
   
   int   ID()          {return _ID;};
   
   void sel(int sync);
   
   float E()           {return _E;};
   float pt()          {return _pt;};
   float eta()         {return _eta;};
   float phi()         {return _phi;};
   float m()           {return _m;};
   float qg()          {return _qg;};
   
   bool  isTight()     {return _isTight;};
   bool  tightJetID()     {return _tightJetID;};
   bool  isLooseTTH()     {return _isLooseTTH;};
   
   int   ntrk()        {return _ntrk;};
   
   float CSVv2()       {return _CSVv2;};
   float cMVAv2()      {return _cMVAv2;};
   float deepCSV()     {return _deepCSV;};
   float deepCSVudsg() {return _deepCSVudsg;};
   float deepCSVb()    {return _deepCSVb;};
   float deepCSVbb()   {return _deepCSVbb;};
   float deepCSVc()    {return _deepCSVc;};
   float deepCSVcc()   {return _deepCSVcc;};
   
   float jet_partonFlavour() {return _jet_partonFlavour  ;};
   float jet_hadronFlavour() {return _jet_hadronFlavour  ;};
   
   float jet_genJet_pt()     {return  _jet_genJet_pt    ;};
   float jet_genJet_E()      {return  _jet_genJet_E   ;};
   /*float jet_genJet_eta()    {return  _jet_genJet_eta ;};
    float jet_genJet_phi()    {return  _jet_genJet_phi ;};
    float jet_genJet_m()      {return  _jet_genJet_m   ;};
    float jet_genJet_E()      {return  _jet_genJet_E   ;};
    float jet_genJet_status() {return  _jet_genJet_status ;};
    float jet_genJet_id()     {return  _jet_genJet_id  ;};*/
   
   float jet_genParton_pt()     {return  _jet_genParton_pt     ;};
   float jet_genParton_E()      {return  _jet_genParton_E      ;};
   /*float jet_genParton_eta()    {return  _jet_genParton_eta    ;};
    float jet_genParton_phi()    {return  _jet_genParton_phi    ;};
    float jet_genParton_m()      {return  _jet_genParton_m      ;};
    float jet_genParton_E()      {return  _jet_genParton_E      ;};
    float jet_genParton_status() {return  _jet_genParton_status ;};*/
   float jet_genParton_id()     {return  _jet_genParton_id     ;};
   
   float JES_uncert()	{return _JES_uncert;}
   float pt_JER()		{return _pt_JER;}
   float pt_JER_down()	{return _pt_JER_down;}
   float pt_JER_up()	{return _pt_JER_up;}

   bool isLooseBTag() {return _isLooseBTag;}
   bool isMediumBTag() {return _isMediumBTag;}
   bool isTightBTag() {return _isTightBTag;}
   
   void  read(bool isdata);
   void  init();
   void  setJESUncertainty(float JES_uncert);
   void  JECUncertainty();
   
 protected:
   
   int _ID;
   
   float _E;
   float _pt;
   float _eta;
   float _phi;
   float _m;
   float _qg;
   
   bool _isTight;
   bool _tightJetID;
   bool _isLooseTTH;
   
   int _ntrk;
   
   float _CSVv2;
   float _cMVAv2;
   float _deepCSV;
   float _deepCSVudsg;
   float _deepCSVb;
   float _deepCSVbb;
   float _deepCSVc;
   float _deepCSVcc;

   bool _isLooseBTag;
   bool _isMediumBTag;
   bool _isTightBTag;
   
   float _jet_partonFlavour  ;
   float _jet_hadronFlavour  ;
   
   float _jet_genJet_pt      ;
   float _jet_genJet_E       ;
   /*float _jet_genJet_eta     ;
    float _jet_genJet_phi     ;
    float _jet_genJet_m       ;
    float _jet_genJet_E       ;
    float _jet_genJet_status  ;
    float _jet_genJet_id      ;*/
   
   float _jet_genParton_pt     ;
   float _jet_genParton_E      ;
   /*float _jet_genParton_eta    ;
    float _jet_genParton_phi    ;
    float _jet_genParton_m      ;
    float _jet_genParton_E      ;
    float _jet_genParton_status ;*/
   float _jet_genParton_id     ;
   
   float _JES_uncert;
   
   float _pt_JER;
   float _pt_JER_down;
   float _pt_JER_up;
   
   ClassDef(Jet,1)
};

#endif
