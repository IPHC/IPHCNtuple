#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(Jet)

Jet::Jet()
{
}

Jet::~Jet()
{
}

void Jet::read(bool isdata)
{
   _ID = idx;
   
   if( CHECK(ntP->jet_E)       )       _E          = ntP->jet_E->at(idx);
   if( CHECK(ntP->jet_pt)      )       _pt         = ntP->jet_pt->at(idx);
   if( CHECK(ntP->jet_eta)     )       _eta        = ntP->jet_eta->at(idx);   
   if( CHECK(ntP->jet_phi)     )       _phi        = ntP->jet_phi->at(idx);   
   if( CHECK(ntP->jet_m)       )       _m          = ntP->jet_m->at(idx);
   if( CHECK(ntP->jet_qgtag)   )       _qg         = ntP->jet_qgtag->at(idx);

   if( CHECK(ntP->jet_tightJetID) )    _tightJetID    = ntP->jet_tightJetID->at(idx);
   
   if ( CHECK(ntP->jet_ntrk)    )          _ntrk           = ntP->jet_ntrk->at(idx);
   if ( CHECK(ntP->jet_CSVv2)   )          _CSVv2          = ntP->jet_CSVv2->at(idx);
   if ( CHECK(ntP->jet_cMVAv2)  )          _cMVAv2         = ntP->jet_cMVAv2->at(idx);
   if ( CHECK(ntP->jet_DeepCSVProbudsg) )  _deepCSVudsg    = ntP->jet_DeepCSVProbudsg->at(idx);
   if ( CHECK(ntP->jet_DeepCSVProbb) )     _deepCSVb       = ntP->jet_DeepCSVProbb->at(idx);
   if ( CHECK(ntP->jet_DeepCSVProbbb) )    _deepCSVbb      = ntP->jet_DeepCSVProbbb->at(idx);
   if ( CHECK(ntP->jet_DeepCSVProbc) )     _deepCSVc       = ntP->jet_DeepCSVProbc->at(idx);
   if ( CHECK(ntP->jet_DeepCSVProbcc) )    _deepCSVcc      = ntP->jet_DeepCSVProbcc->at(idx);
   
   if (!isdata)
     {
        _jet_partonFlavour    = ntP->jet_partonFlavour->at(idx);
        _jet_hadronFlavour    = ntP->jet_hadronFlavour->at(idx);
        _jet_genJet_pt        = ntP->jet_genJet_pt->at(idx); 
        _jet_genJet_E         = ntP->jet_genJet_E->at(idx);
        _jet_genParton_pt     = ntP->jet_genParton_pt->at(idx);
        _jet_genParton_id     = ntP->jet_genParton_id ->at(idx);
        _jet_genParton_E      = ntP->jet_genParton_E->at(idx);      
     }
   
   JECUncertainty();   
}

void Jet::init()
{
   _E       = -100.;
   _pt      = -100.;
   _eta     = -100.;
   _phi     = -100.;
   _m       = -100.;
   _qg      = -100.;
   
   _tightJetID = 0;
   _isLooseTTH = 0;
   
   _ntrk           = -100;
   _CSVv2          = -100.;
   _cMVAv2         = -100.;
   
   _deepCSV        = -100.;
   _deepCSVudsg    = -100.;
   _deepCSVb       = -100.;
   _deepCSVbb      = -100.;
   _deepCSVc       = -100.;
   _deepCSVcc      = -100.;

   _isLooseBTag = 0;
   _isMediumBTag = 0;
   _isTightBTag = 0;
   
   _jet_partonFlavour    = -100.;
   _jet_hadronFlavour    = -100.;
   _jet_genJet_pt        = -100.;
   _jet_genParton_pt     = -100.;
   _jet_genParton_id     = -100.; 
   _jet_genJet_E         = -100.;
   _jet_genParton_E      = -100.; 
   
   _JES_uncert           = 0.;
   _pt_JER               = 0.;
   _pt_JER_down          = 0.;
   _pt_JER_up            = 0.;
}

void Jet::sel(int sync)
{
   float jet_pt_JESup = _pt*(1+_JES_uncert);
   
   bool pass_pt      = ( _pt > 25. || jet_pt_JESup > 25. );
   bool pass_eta     = (fabs(_eta) < 2.4);
   bool pass_tightJetID = (_tightJetID);

   bool pass_lepOverlap = 1;
   bool pass_tauOverlap = 1;
   
   int nMuon = nt->NtMuonFakeable->size();
   int nElectron = nt->NtElectronFakeable->size();
   int nTau = nt->NtTauFakeable->size();
   int nTauTight = nt->NtTauTight->size();
   int nLepTight = nt->NtMuonTight->size()+nt->NtElectronTight->size();
   
   std::vector<Base> *elmuFakeable = new std::vector<Base>;
   for(int ie=0;ie<nElectron;ie++ )
     {
	Base lep = nt->NtElectronFakeable->at(ie);
	lep.iElec = ie;
	elmuFakeable->push_back(lep);
     }
   for(int im=0;im<nMuon;im++ )
     {
	Base lep = nt->NtMuonFakeable->at(im);
	lep.iMuon = im;
	elmuFakeable->push_back(lep);
     }
   std::sort(elmuFakeable->begin(),elmuFakeable->end(),sort_by_pt());
   
   int nLep = elmuFakeable->size();

   bool is_1l2tau = (nLep > 0 && nLepTight <= 1 && nTau >= 2);
   bool is_2lSS = (nLep > 1 && nLepTight <= 2 && nTauTight == 0);
   bool is_2lSS1tau = (nLep > 1 && nLepTight <= 2 && nTauTight > 0);
   bool is_2l2tau = (nLep > 1 && nTau > 1);
   bool is_3l = (nLep > 2 && nTauTight == 0);
   bool is_3l1tau = (nLep > 2 && nTau > 0);
   bool is_4l = (nLep >= 4);
   
   for(int il=0;il<nLep;il++)
     {
	if( is_1l2tau && il >= 1 ) break;
	if( is_2lSS && il >= 2 ) break;
	if( is_2lSS1tau && il >= 2 ) break;
	if( is_2l2tau && il >= 2 ) break;
	if( is_3l && il >= 3 ) break;
	if( is_3l1tau && il >= 3 ) break;
	if( is_4l && il >= 4 ) break;
	
        float dr = GetDeltaR(_eta,_phi,elmuFakeable->at(il)._eta,elmuFakeable->at(il)._phi);
        if( dr < 0.4 ) pass_lepOverlap = 0;
     }      

   for(int it=0;it<nTau;it++)
     {
	if( is_1l2tau && it >= 2 ) break;
	if( is_2lSS ) break;
	if( is_2lSS1tau && it >= 1 ) break;
	if( is_2l2tau && it >= 2 ) break;
	if( is_3l ) break;
	if( is_3l1tau && it >= 1 ) break;
	if( is_4l ) break;

        float dr = GetDeltaR(_eta,_phi,nt->NtTauFakeable->at(it).eta(),nt->NtTauFakeable->at(it).phi());
        if( dr < 0.4 ) pass_tauOverlap = 0;
     }      
   
   delete elmuFakeable;
   
   pass_lepOverlap = (pass_lepOverlap && sync != 1) || (sync == 1);
   pass_tauOverlap = (pass_tauOverlap && sync != 1) || (sync == 1);
   
   _isLooseTTH = ( pass_pt &&
		   pass_eta &&
		   pass_tightJetID &&
		   pass_lepOverlap &&
		   pass_tauOverlap
		 );
   
   _isLooseBTag  =  (_deepCSVb+_deepCSVbb > 0.1522);
   _isMediumBTag =  (_deepCSVb+_deepCSVbb > 0.4941);
   _isTightBTag  =  (_deepCSVb+_deepCSVbb > 0.8001);
}

void Jet::setJESUncertainty(float JES_uncert)
{
   _JES_uncert = JES_uncert;
}

void Jet::JECUncertainty()
{  
   // JER taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
   // 
   double _cJER[5];
   double _cJER_down[5];
   double _cJER_up[5];
   
   _cJER[0] = 1.052; // 0.0-0.5
   _cJER[1] = 1.057; // 0.5-1.1
   _cJER[2] = 1.096; // 1.1-1.7
   _cJER[3] = 1.134; // 1.7-2.3
   _cJER[4] = 1.288; // 2.3-5.0
   
   _cJER_down[0] = 0.990;
   _cJER_down[1] = 1.001;
   _cJER_down[2] = 1.032;
   _cJER_down[3] = 1.042;
   _cJER_down[4] = 1.089;

   _cJER_up[0] = 1.115;
   _cJER_up[1] = 1.114;
   _cJER_up[2] = 1.161;
   _cJER_up[3] = 1.228;
   _cJER_up[4] = 1.488;
   
   int etaIdx = -1;
   if( fabs(_eta) >= 0.  && fabs(_eta) < 0.5 ) etaIdx = 0;
   if( fabs(_eta) >= 0.5 && fabs(_eta) < 1.1 ) etaIdx = 1;
   if( fabs(_eta) >= 1.1 && fabs(_eta) < 1.7 ) etaIdx = 2;
   if( fabs(_eta) >= 1.7 && fabs(_eta) < 2.3 ) etaIdx = 3;
   if( fabs(_eta) >= 2.3 && fabs(_eta) < 5.0 ) etaIdx = 4;	

   double pt_uncorr = (_pt - _jet_genJet_pt) / _cJER[etaIdx] + _jet_genJet_pt;
  
   _pt_JER	= _jet_genJet_pt + _cJER[etaIdx]*(pt_uncorr-_jet_genJet_pt);
   _pt_JER_down = _jet_genJet_pt + _cJER_down[etaIdx]*(pt_uncorr-_jet_genJet_pt);
   _pt_JER_up   = _jet_genJet_pt + _cJER_up[etaIdx]*(pt_uncorr-_jet_genJet_pt);		    
}
