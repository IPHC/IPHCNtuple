#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(JetExt)

JetExt::JetExt()
{
}

JetExt::~JetExt()
{
}

void JetExt::read(bool isdata)
{
   ID = idx;

   E          = ntP->jet_E->at(idx);
   pt         = ntP->jet_pt->at(idx);
   eta        = ntP->jet_eta->at(idx);
   phi        = ntP->jet_phi->at(idx);
   m          = ntP->jet_m->at(idx);
   qgtag      = ntP->jet_qgtag->at(idx);
   pileupJetId = ntP->jet_pileupJetId->at(idx);

   tightJetID    = ntP->jet_tightJetID->at(idx);

   ntrk           = ntP->jet_ntrk->at(idx);
   CSVv2          = ntP->jet_CSVv2->at(idx);
   cMVAv2         = ntP->jet_cMVAv2->at(idx);
   deepCSVudsg    = ntP->jet_DeepCSVProbudsg->at(idx);
   deepCSVb       = ntP->jet_DeepCSVProbb->at(idx);
   deepCSVbb      = ntP->jet_DeepCSVProbbb->at(idx);
   deepCSVc       = ntP->jet_DeepCSVProbc->at(idx);
   deepCSVcc      = ntP->jet_DeepCSVProbcc->at(idx);

   DeepCSVbtag	  = ntP->jet_DeepCSVProbb->at(idx) + ntP->jet_DeepCSVProbbb->at(idx);

   if( !isdata )
     {
        jet_partonFlavour    = ntP->jet_partonFlavour->at(idx);
        jet_hadronFlavour    = ntP->jet_hadronFlavour->at(idx);
        jet_genJet_pt        = ntP->jet_genJet_pt->at(idx);
        jet_genJet_E         = ntP->jet_genJet_E->at(idx);
        jet_genParton_pt     = ntP->jet_genParton_pt->at(idx);
        jet_genParton_id     = ntP->jet_genParton_id ->at(idx);
        jet_genParton_E      = ntP->jet_genParton_E->at(idx);
     }

   // JECUncertainty();
}

void JetExt::init()
{
   E       = -100.;
   pt      = -100.;
   eta     = -100.;
   phi     = -100.;
   m       = -100.;
   qgtag      = -100.;
   pileupJetId = -100.;

   tightJetID = 0;
   isLooseTTH = 0;

   ntrk           = -100;
   CSVv2          = -100.;
   cMVAv2         = -100.;

   deepCSV        = -100.;
   deepCSVudsg    = -100.;
   deepCSVb       = -100.;
   deepCSVbb      = -100.;
   deepCSVc       = -100.;
   deepCSVcc      = -100.;

   DeepCSVbtag    = -100;

   isLooseBTag = 0;
   isMediumBTag = 0;
   isTightBTag = 0;

   jet_partonFlavour    = -100.;
   jet_hadronFlavour    = -100.;
   jet_genJet_pt        = -100.;
   jet_genParton_pt     = -100.;
   jet_genParton_id     = -100.;
   jet_genJet_E         = -100.;
   jet_genParton_E      = -100.;

   JES_uncert           = 0.;
   pt_JER_down          = 0.;
   pt_JER_up            = 0.;
}

void JetExt::sel(int sync)
{
   float jet_pt_JESup = pt*(1+JES_uncert);

   bool pass_pt = (pt > 25.);

   bool pass_eta = false;

   if(sync == 0)
   {
	 pass_pt = pass_pt || (jet_pt_JESup > 25.);
	 pass_eta = (fabs(eta) < 5.0); //tHq2017 //Default cut now (always store tHq events)
   }
   else //For synchro with ttH, never care about fwd jets
   {
   	pass_eta = (fabs(eta) < 2.4); //ttH2017
   }

   bool pass_tightJetID = (tightJetID);

   bool pass_lepOverlap = 1;
   bool pass_tauOverlap = 1;

   int nMuon = nt->NtMuonFakeableExt->size();
   int nElectron = nt->NtElectronFakeableExt->size();
   int nTau = nt->NtTauFakeableExt->size();
   int nTauLoose = nt->NtTauLooseExt->size();
   int nTauMedium = nt->NtTauMediumExt->size();
   int nLepTight = nt->NtMuonTightExt->size()+nt->NtElectronTightExt->size();

   std::vector<Base> *elmuFakeable = new std::vector<Base>;
   for(int ie=0;ie<nElectron;ie++ )
     {
	Base lep = nt->NtElectronFakeableExt->at(ie);
	lep.iElec = ie;
	elmuFakeable->push_back(lep);
     }
   for(int im=0;im<nMuon;im++ )
     {
	Base lep = nt->NtMuonFakeableExt->at(im);
	lep.iMuon = im;
	elmuFakeable->push_back(lep);
     }
   std::sort(elmuFakeable->begin(),elmuFakeable->end(),sort_by_conept());

   int nLep = elmuFakeable->size();

   /*
   bool is_1l2tau = (nLep > 0 && nLepTight <= 1 && nTau >= 2);
   bool is_2lSS = (nLep > 1 && nLepTight <= 2 && nTauLoose == 0);
   bool is_2lSS1tau = (nLep > 1 && nLepTight <= 2 && nTauLoose > 0 && nTauMedium <= 1);
   bool is_2l2tau = (nLep > 1 && nLepTight <= 2 && nTau > 1);
   bool is_3l = (nLep > 2 && nLepTight <= 3 && nTauLoose == 0);
   bool is_3l1tau = (nLep > 2 && nLepTight <= 3 && nTau > 0);
   bool is_4l = (nLep > 3 && nLepTight > 3);
   */

   for(int il=0;il<nLep;il++)
     {
     	//Sergio (Oviedo/CERN) confirmed that they actually clean jets w.r.t *all* FO leptons/taus in *all* categs
     	/*
	if( is_1l2tau && il >= 1 ) break;
	if( is_2lSS && il >= 2 ) break;
	if( is_2lSS1tau && il >= 2 ) break;
	if( is_2l2tau && il >= 2 ) break;
	if( is_3l && il >= 3 ) break;
	if( is_3l1tau && il >= 3 ) break;
	if( is_4l && il >= 4 ) break;
	*/

        float dr = GetDeltaR(eta,phi,elmuFakeable->at(il).eta,elmuFakeable->at(il).phi);
        if( dr < 0.4 ) pass_lepOverlap = 0;
     }

   for(int it=0;it<nTau;it++)
     {
     	//Sergio (Oviedo/CERN) confirmed that they actually clean jets w.r.t *all* FO leptons/taus in *all* categs
     	/*
	if( is_1l2tau && it >= 2 ) break;
	if( is_2lSS ) break;
	if( is_2lSS1tau && it >= 1 ) break;
	if( is_2l2tau && it >= 2 ) break;
	if( is_3l ) break;
	if( is_3l1tau && it >= 1 ) break;
	if( is_4l ) break;
	*/


        float dr = GetDeltaR(eta,phi,nt->NtTauFakeableExt->at(it).eta,nt->NtTauFakeableExt->at(it).phi);
        if( dr < 0.4 ) pass_tauOverlap = 0;
     }

   delete elmuFakeable;

//   pass_lepOverlap = (pass_lepOverlap && sync != 1) || (sync == 1);
//   pass_tauOverlap = (pass_tauOverlap && sync != 1) || (sync == 1);

   isLooseTTH = ( pass_pt &&
		  pass_eta &&
		  pass_tightJetID &&
		  pass_lepOverlap &&
		  pass_tauOverlap
		);

   isLooseBTag  =  (deepCSVb+deepCSVbb > 0.1522);
   isMediumBTag =  (deepCSVb+deepCSVbb > 0.4941);
   isTightBTag  =  (deepCSVb+deepCSVbb > 0.8001);

   for(int d=0;d<evdebug->size();d++)
     {
	double evId = ntP->ev_id;
	if( evId == evdebug->at(d) )
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << evId << std::endl;
	     std::cout << "  jet #" << ID << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  eta = " << eta << std::endl;
	     std::cout << "  phi = " << phi << std::endl;
	     std::cout << "  jet_pt_JESup = " << jet_pt_JESup << std::endl;
	     std::cout << "  isLooseTTH = " << isLooseTTH << std::endl;
	     std::cout << "  pass_pt = " << pass_pt << std::endl;
	     std::cout << "  pass_eta = " << pass_eta << std::endl;
	     std::cout << "  tightJetID = " << tightJetID << std::endl;
	     std::cout << "  pass_lepOverlap = " << pass_lepOverlap << std::endl;
	     std::cout << "  pass_tauOverlap = " << pass_tauOverlap << std::endl;
	  }
     }
}

void JetExt::setJESUncertainty(float unc)
{
   JES_uncert = unc;
}


void JetExt::apply_JER_smearing(bool isdata, float JER_res, float JER_sf, float JER_sf_up, float JER_sf_down)
{
    if(isdata) {JES_uncert = 0; pt_JER_up = pt; pt_JER_down = pt;}

    //Check if jet is gen matched
    float dpt = 0;
    float dR = 0;
    float dpt_min = 99999;
    int idx_matched_genjet = -1;
	float resol = JER_res*pt;

    for(int ij=0; ij<ntP->genJet_n; ij++)
    {
        dpt = fabs(pt - ntP->genJet_pt->at(ij));
        dR = GetDeltaR(eta, phi, ntP->genJet_eta->at(ij), ntP->genJet_phi->at(ij));

        if(dR < 0.2 && dpt < (3*fabs(resol)) )
        {
            if( dpt <= dpt_min )
            {
                idx_matched_genjet = ij;
                dpt_min = dpt;
            }
        }
    }

    float JER_corr      = 1.;
    float JER_corr_up   = 1.;
    float JER_corr_down = 1.;

    if(idx_matched_genjet >= 0) //if jet matched
    {
        float genpt = ntP->genJet_pt->at(idx_matched_genjet);

        if(genpt >= 0.)
        {
            JER_corr = 1 + (JER_sf - 1) * (pt - genpt) / pt;
            JER_corr_up = 1 + (JER_sf_up - 1) * (pt - genpt) / pt;
            JER_corr_down = 1 + (JER_sf_down - 1) * (pt - genpt) / pt;
        }
    }
    else
    {
        TRandom3 rnd;
        float smear = rnd.Gaus(0., JER_res);

        float JER_corr = 1 + smear * sqrt(std::max((float) 0., JER_sf*JER_sf - 1) );
        float JER_corr_up = 1 + smear * sqrt(std::max((float) 0., JER_sf_up*JER_sf_up - 1) );
        float JER_corr_down = 1 + smear * sqrt(std::max((float) 0., JER_sf_down*JER_sf_down - 1) );
    }

    if(JER_corr<0) {JER_corr = 1;;}
    if(JER_corr_up<0) {JER_corr_up = 1.;}
    if(JER_corr_down<0) {JER_corr_down = 1.;}

    //Store up/down variations
    pt_JER_up = pt * JER_corr_up;
    pt_JER_down = pt * JER_corr_down;

    //Correct jet pt by JER correcting factor
    pt*= JER_corr;

    return;
}


//Obsolete !
/*
void JetExt::JECUncertainty()
{
   // JER taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
   //
   double cJER[5];
   double cJER_down[5];
   double cJER_up[5];

   cJER[0] = 1.052; // 0.0-0.5
   cJER[1] = 1.057; // 0.5-1.1
   cJER[2] = 1.096; // 1.1-1.7
   cJER[3] = 1.134; // 1.7-2.3
   cJER[4] = 1.288; // 2.3-5.0

   cJER_down[0] = 0.990;
   cJER_down[1] = 1.001;
   cJER_down[2] = 1.032;
   cJER_down[3] = 1.042;
   cJER_down[4] = 1.089;

   cJER_up[0] = 1.115;
   cJER_up[1] = 1.114;
   cJER_up[2] = 1.161;
   cJER_up[3] = 1.228;
   cJER_up[4] = 1.488;

   int etaIdx = -1;
   if( fabs(eta) >= 0.  && fabs(eta) < 0.5 ) etaIdx = 0;
   if( fabs(eta) >= 0.5 && fabs(eta) < 1.1 ) etaIdx = 1;
   if( fabs(eta) >= 1.1 && fabs(eta) < 1.7 ) etaIdx = 2;
   if( fabs(eta) >= 1.7 && fabs(eta) < 2.3 ) etaIdx = 3;
   if( fabs(eta) >= 2.3 && fabs(eta) < 5.0 ) etaIdx = 4;

   double pt_uncorr = (pt - jet_genJet_pt) / cJER[etaIdx] + jet_genJet_pt;

   pt_JER	= jet_genJet_pt + cJER[etaIdx]*(pt_uncorr-jet_genJet_pt);
   pt_JER_down = jet_genJet_pt + cJER_down[etaIdx]*(pt_uncorr-jet_genJet_pt);
   pt_JER_up   = jet_genJet_pt + cJER_up[etaIdx]*(pt_uncorr-jet_genJet_pt);
}
*/
