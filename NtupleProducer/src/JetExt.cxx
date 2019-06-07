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
   QGdiscr    = ntP->jet_qgtag->at(idx);
   pileupJetId = ntP->jet_pileupJetId->at(idx);

   tightJetID    = ntP->jet_tightJetID->at(idx);
   looseJetID    = ntP->jet_looseJetID->at(idx);

   ntrk           = ntP->jet_ntrk->at(idx);
   CSVv2          = ntP->jet_CSVv2->at(idx);
   cMVAv2         = ntP->jet_cMVAv2->at(idx);

   deepCSVudsg    = ntP->jet_DeepCSVProbudsg->at(idx);
   deepCSVb       = ntP->jet_DeepCSVProbb->at(idx);
   deepCSVbb      = ntP->jet_DeepCSVProbbb->at(idx);
   deepCSVc       = ntP->jet_DeepCSVProbc->at(idx);
   deepCSVcc      = ntP->jet_DeepCSVProbcc->at(idx);

   DeepCSVbtag	  = deepCSVb + deepCSVbb;

   deepFlavourb       = ntP->jet_DeepFlavourProbb->at(idx);
   deepFlavourbb      = ntP->jet_DeepFlavourProbbb->at(idx);
   deepFlavourlepb    = ntP->jet_DeepFlavourProblepb->at(idx);
   deepFlavourc       = ntP->jet_DeepFlavourProbc->at(idx);
   deepFlavouruds     = ntP->jet_DeepFlavourProbuds->at(idx);
   deepFlavourg       = ntP->jet_DeepFlavourProbg->at(idx);

   DeepFlavourbtag    = deepFlavourb + deepFlavourbb + deepFlavourlepb;
   
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
}

void JetExt::init()
{
   E       = -100.;
   pt      = -100.;
   eta     = -100.;
   phi     = -100.;
   m       = -100.;
   QGdiscr      = -100.;
   pileupJetId = -100.;

   tightJetID = 0;
   looseJetID = 0;
   isLooseTTH = 0;
   isLooseFwdTTH = 0;
   isSoftLooseTTH = 0;

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

   deepFlavourb      = -100.;
   deepFlavourbb     = -100.;
   deepFlavourlepb   = -100.;
   deepFlavourc      = -100.;
   deepFlavouruds    = -100.;
   deepFlavourg      = -100.;
   
   DeepFlavourbtag    = -100;
   
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
   
   JER_corr = 1;
   JER_corr_up = 1;
   JER_corr_down = 1;

   pt_JES_up           = 0.;
   pt_JES_down         = 0.;
   pt_JER_up           = 0.;
   pt_JER_down         = 0.;
   
   E_JES_up           = 0.;
   E_JES_down         = 0.;
   E_JER_up           = 0.;
   E_JER_down         = 0.;
}

void JetExt::sel(int sync, bool DEBUG, int year)
{
   bool pass_pt = (pt > 25.);

   bool pass_eta = false;
   bool pass_eta_fwd = (fabs(eta) > 2.4 && fabs(eta) < 5.0);
   bool pass_pt_fwd = (fabs(eta) <= 2.7 || fabs(eta) >= 3.0 || (fabs(eta) > 2.7 && fabs(eta) < 3.0 && pt > 60.));
   
   if(sync == 0)
   {
	 pass_pt = pass_pt || pt_JES_up > 25. || pt_JER_up > 25.; //Also keep jets with variations >25 GeV, for systematics
	 pass_eta = (fabs(eta) < 5.0); //tHq2017 //Default cut (to keep events of tHq analysis) -- Decide in NTAnalyzer if want to use fwd jets or not
	 //pass_eta = (fabs(eta) < 2.4); //ttH cut -- not used here for now
   }
   else //For synchro with ttH, never care about fwd jets
   {
	pass_pt = (pt > 25.);
   	pass_eta = (fabs(eta) < 2.4); //ttH cut
   }

   bool pass_tightJetID = (tightJetID);
   bool pass_looseJetID = (looseJetID);
   bool pass_JetID = ((year == 2016 && looseJetID) || (year != 2016 && tightJetID));

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
        if( dr < 0.4 && ID == elmuFakeable->at(il).matchedJetId ) pass_lepOverlap = 0;
     }

   for(int it=0;it<nTau;it++)
     {
     	//CERN actually clean jets w.r.t *all* FO leptons/taus in *all* categs
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
        if( dr < 0.4 && ID == nt->NtTauFakeableExt->at(it).matchedJetId ) pass_tauOverlap = 0;
     }

   delete elmuFakeable;

//   pass_lepOverlap = (pass_lepOverlap && sync != 1) || (sync == 1);
//   pass_tauOverlap = (pass_tauOverlap && sync != 1) || (sync == 1);

   isLooseTTH = ( pass_pt &&
		  pass_eta &&
		  pass_JetID &&
		  pass_lepOverlap &&
		  pass_tauOverlap
		);

   isLooseFwdTTH = ( pass_pt &&
		     pass_pt_fwd &&
		     pass_eta_fwd &&
		     pass_JetID &&
		     pass_lepOverlap &&
		     pass_tauOverlap
		);
   
   //NEW -- similar as "isLooseTTH" sel, but for jets with pT between 15 and 25 => Save nof soft jets, for FCNC input variable
   isSoftLooseTTH = ( pt > 15 && pt < 25 &&
		  pass_eta &&
		  pass_JetID &&
		  pass_lepOverlap &&
		  pass_tauOverlap
		);

   if( year == 2016 )
     {	
	isLooseBTag  =  (DeepFlavourbtag > 0.0614);
	isMediumBTag =  (DeepFlavourbtag > 0.3093);
	isTightBTag  =  (DeepFlavourbtag > 0.7221);
     }
   else if( year == 2017 )
     {
	isLooseBTag  =  (DeepFlavourbtag > 0.0521);
	isMediumBTag =  (DeepFlavourbtag > 0.3033);
	isTightBTag  =  (DeepFlavourbtag > 0.7489);
     }   
   else
     {
	isLooseBTag  =  (DeepFlavourbtag > 0.0494);
	isMediumBTag =  (DeepFlavourbtag > 0.2770);
	isTightBTag  =  (DeepFlavourbtag > 0.7264);
     }   

	  if(DEBUG)
	  {
	     std::cout << "------------------------------" << std::endl;
	     std::cout << "Event #" << std::setprecision(12)<< ntP->ev_id << std::endl;
	     std::cout << "  jet #" << ID << std::endl;
	     std::cout << "  pt = " << pt << std::endl;
	     std::cout << "  eta = " << eta << std::endl;
	     std::cout << "  phi = " << phi << std::endl;
	     std::cout << "  pt_JES_up = " << pt_JES_up << std::endl;
	     std::cout << "  isLooseTTH = " << isLooseTTH << std::endl;
	     std::cout << "  isLooseFwdTTH = " << isLooseFwdTTH << std::endl;
	     std::cout << "  pass_pt = " << pass_pt << std::endl;
	     std::cout << "  pass_eta = " << pass_eta << std::endl;
	     std::cout << "  DeepFlavourbtag = " << DeepFlavourbtag << std::endl;
	     std::cout << "  isLooseBTag = " << isLooseBTag << std::endl;
	     std::cout << "  isMediumBTag = " << isMediumBTag << std::endl;
	     std::cout << "  JetID = " << pass_JetID << std::endl;
	     std::cout << "  (( jet_neutralEmEnergyFraction = " <<
	     ntP->jet_neutralEmEnergyFraction->at(idx) << std::endl;
	     std::cout << "  pass_lepOverlap = " << pass_lepOverlap << std::endl;
	     std::cout << "  pass_tauOverlap = " << pass_tauOverlap << std::endl;
	  }
     
     return;
}

//[0]: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
//[1] : https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetSmearer.py
//[2] : https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
//NB : found some inconsistencies in docs/codes : 
//1) In [0], it is said to do for non-matched jets : JER_corr = 1 + gaus(0, JER_res) * srt(JER_sf*JER_sf-1). But in [1,2], they do : JER_corr = gaus(0, JER_res*srt(JER_sf*JER_sf-1)) !
//=> Looks like should follow stochastic smearing formula from [1,2] instead (looks more reasonnable)
void JetExt::apply_JER_smearing(bool isdata, float JER_res, float JER_sf, float JER_sf_up, float JER_sf_down)
{
    bool debug = false;
    
    if(debug) {cout<<endl<<endl<<"Debug... :"<<endl;}

    if(isdata) {pt_JER_up = pt; pt_JER_down = pt; return;}
    
    if(JER_sf < 0) {JER_sf = 0;}
    if(JER_sf_up < 0) {JER_sf_up = 0;}
    if(JER_sf_down < 0) {JER_sf_down = 0;}

    //Check if jet is gen matched
    float dpt = 0;
    float dR = 0;
    float dR_min = 99999;
    int idx_matched_genjet = -1;
    float resol = JER_res*pt;
    
    if(debug) {cout<<"pt = "<<pt<<endl;}
    if(debug) {cout<<"eta = "<<eta<<endl;}
    if(debug) {cout<<"JER_sf = "<<JER_sf<<endl;}
    if(debug) {cout<<"JER_sf_up = "<<JER_sf_up<<endl;}
    if(debug) {cout<<"JER_sf_down = "<<JER_sf_down<<endl;}
    if(debug) {cout<<"JER_res = "<<JER_res<<endl;}
    if(debug) {cout<<"resol = "<<resol<<endl;}

    if(debug) {cout<<"------"<<endl;}
    for(int ij=0; ij<ntP->genJet_n; ij++)
    {
        dpt = fabs(pt - ntP->genJet_pt->at(ij));
        dR = GetDeltaR(eta, phi, ntP->genJet_eta->at(ij), ntP->genJet_phi->at(ij));
	
	if(debug) {cout<<"dpt = "<<dpt<<endl;}
	if(debug) {cout<<"dR = "<<dR<<endl;}

        if(dR < (0.4/2) && dpt < (3 * fabs(resol)) ) //What should be used ? cf func comments
        //if(dR < 0.4 && dpt < (3 * fabs(resol)) ) //What should be used ? cf func comments
        {
            if(dR < dR_min)
            {
                idx_matched_genjet = ij;
                dR_min = dR;
            }
        }
    }
    if(debug) {cout<<"------"<<endl;}

    JER_corr      = 1.;
    JER_corr_up   = 1.;
    JER_corr_down = 1.;

    if(idx_matched_genjet >= 0) //if jet matched
    {   
        if(debug) {cout<<"jet matched OK"<<endl;}

        float genpt = ntP->genJet_pt->at(idx_matched_genjet);

        if(genpt >= 0.)
        {
            JER_corr = 1 + (JER_sf - 1) * (pt - genpt) / pt;
            JER_corr_up = 1 + (JER_sf_up - 1) * (pt - genpt) / pt;
            JER_corr_down = 1 + (JER_sf_down - 1) * (pt - genpt) / pt;
        }
    }
    /*
    else if(JER_sf > 1)
    {   
    	//Case 1: we have a "good" gen jet matched to the reco jet //NB : as advised on Twiki ! But gives too large corrections, must be wrong !
    
        TRandom3 rnd(0); //Must set arg to 0 to get unique seed !
        float smear = rnd.Gaus(0., JER_res);
	
        JER_corr = 1 + smear * sqrt(std::max((float) 0., JER_sf*JER_sf - 1) );
        JER_corr_up = 1 + smear * sqrt(std::max((float) 0., JER_sf_up*JER_sf_up - 1) );
        JER_corr_down = 1 + smear * sqrt(std::max((float) 0., JER_sf_down*JER_sf_down - 1) );

	if(debug) {cout<<"Jet not matched XX ! smear with GAUS(0, "<<sqrt(JER_res)<<")="<<smear<<" * "<<sqrt(std::max((float) 0., JER_sf*JER_sf - 1))<<" ==> "<<JER_corr<<endl;}
    }
    */
    else if(JER_sf > 1)
    {   
    	//Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation //NB : as found in other codes ! Reasonnable corrections

        TRandom3 rnd(0);
	//rnd.SetSeed((uint32_t)(jet.userInt("deterministicSeed")));
        float smear = rnd.Gaus(0., JER_res * sqrt(std::max((float) 0., JER_sf*JER_sf - 1) ));
        float smear_up = rnd.Gaus(0., JER_res * sqrt(std::max((float) 0., JER_sf_up*JER_sf_up - 1) ));
        float smear_down = rnd.Gaus(0., JER_res * sqrt(std::max((float) 0., JER_sf_down*JER_sf_down - 1) ));
	
	if(debug) {cout<<"Jet not matched XX ! smear with GAUS(0, "<<JER_res * sqrt(std::max((float) 0., JER_sf*JER_sf - 1) )<<") ==> "<<smear<<endl;}

        JER_corr = 1 + smear;
        JER_corr_up = 1 + smear_up;
        JER_corr_down = 1 + smear_down;
    }
    else 
    {
    	//Case 3: we cannot smear this jet, as we don't have a generator level jet and the resolution in data is better than the resolution in the simulation,
        //so we would need to randomly "unsmear" the jet, which is impossible
	
	if(debug) {cout<<"Could not smear jet (JER_sf < 1)"<<endl;}
    }
    
    if(debug) {cout<<"JER_corr = "<<JER_corr<<endl;}
    
    if (pt * JER_corr < 0.01) 
    {
     	// Negative or too small smearFactor. We would change direction of the jet and this is not what we want.
        // Recompute the smearing factor in order to have jet.energy() == MIN_JET_ENERGY
        double newSmearFactor = 0.01 / pt;
        if (debug) 
	{
        	std::cout << "The smearing factor (" << JER_corr << ") is either negative or too small. Fixing it to " << newSmearFactor << " to avoid change of direction." << std::endl;
        }
        JER_corr = newSmearFactor;
	JER_corr_up = newSmearFactor;
	JER_corr_down = newSmearFactor;
    }
    
    if(JER_corr<0) {JER_corr = 1.;}
    if(JER_corr_up<0) {JER_corr_up = 1.;}
    if(JER_corr_down<0) {JER_corr_down = 1.;}

    //Store up/down variations
    pt_JER_up = pt * JER_corr_up;
    pt_JER_down = pt * JER_corr_down;
    E_JER_up = E * JER_corr_up;
    E_JER_down = E * JER_corr_down;

    //Correct jet pt/E by JER correcting factor
    pt*= JER_corr;
    E*= JER_corr;

    return;
}

void JetExt::setJESUncertainty(bool isdata, float unc)
{
    bool debug = false;
     
    if(isdata) {pt_JES_up = pt; pt_JES_down = pt;}

    pt_JES_up = pt * (1+fabs(unc));
    pt_JES_down = pt * (1-fabs(unc));
    
    E_JES_up = E * (1+fabs(unc));
    E_JES_down = E * (1-fabs(unc));

    if(debug) //NB : checked that JES uncert. is the same in FlatTreeProducer and NTProd, when not applying JER
    {
    	cout<<"* pt "<<pt<<" / eta "<<eta<<" => FTP : JES = "<<ntP->jet_Unc->at(idx)<<" / NTP : JES = "<<unc<<endl;
    }
    
    return;
}
