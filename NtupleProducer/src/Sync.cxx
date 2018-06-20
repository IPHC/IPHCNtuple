#include "include/Sync.h"
#include "include/NtupleProducer.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

Sync::Sync(std::string fname,int syncFlag)
{
   fname_out = fname;
   sync = syncFlag;
}

Sync::~Sync()
{  
   if( sync != 0 )
     {	
	m_file->Write();
	m_file->Close();
     }   
}

void Sync::Init()
{
   if( sync != 0 ) m_file = new TFile(fname_out.c_str(),"RECREATE");
   if( sync == 1 ) m_tree = new TTree("syncTree","Sync Ntuple");
   if( sync == 2 )
     {
	m_tree_1l2tau_SR = new TTree("syncTree_1l2tau_SR","Sync Ntuple");
	m_tree_1l2tau_Fake = new TTree("syncTree_1l2tau_Fake","Sync Ntuple");
	m_tree_2lSS_SR = new TTree("syncTree_2lSS_SR","Sync Ntuple");
	m_tree_2lSS_Fake = new TTree("syncTree_2lSS_Fake","Sync Ntuple");
	m_tree_2lSS_Flip = new TTree("syncTree_2lSS_Flip","Sync Ntuple");
	m_tree_2lSS1tau_SR = new TTree("syncTree_2lSS1tau_SR","Sync Ntuple");
	m_tree_2lSS1tau_Fake = new TTree("syncTree_2lSS1tau_Fake","Sync Ntuple");
	m_tree_2lSS1tau_Flip = new TTree("syncTree_2lSS1tau_Flip","Sync Ntuple");
	m_tree_2l2tau_SR = new TTree("syncTree_2l2tau_SR","Sync Ntuple");
	m_tree_2l2tau_Fake = new TTree("syncTree_2l2tau_Fake","Sync Ntuple");
	m_tree_3l_SR = new TTree("syncTree_3l_SR","Sync Ntuple");
	m_tree_3l_Fake = new TTree("syncTree_3l_Fake","Sync Ntuple");
	m_tree_3l1tau_SR = new TTree("syncTree_3l1tau_SR","Sync Ntuple");
	m_tree_3l1tau_Fake = new TTree("syncTree_3l1tau_Fake","Sync Ntuple");
	m_tree_4l_SR = new TTree("syncTree_4l_SR","Sync Ntuple");
	m_tree_4l_Fake = new TTree("syncTree_4l_Fake","Sync Ntuple");
	m_tree_ttWctrl_SR = new TTree("syncTree_ttWctrl_SR","Sync Ntuple");
	m_tree_ttWctrl_Fake = new TTree("syncTree_ttWctrl_Fake","Sync Ntuple");
	m_tree_ttWctrl_Flip = new TTree("syncTree_ttWctrl_Flip","Sync Ntuple");
	m_tree_ttZctrl_SR = new TTree("syncTree_ttZctrl_SR","Sync Ntuple");
	m_tree_ttZctrl_Fake = new TTree("syncTree_ttZctrl_Fake","Sync Ntuple");
	
	m_hist_overlap = new TH2F("overlap","overlap",21,0.,21.,21,0.,21.);
     }   
}

void Sync::setBranchAddress()
{
   if( sync == 1 ) createBranch(m_tree);
   if( sync == 2 )
     {
	createBranch(m_tree_1l2tau_SR);
	createBranch(m_tree_1l2tau_Fake);
	createBranch(m_tree_2lSS_SR);
	createBranch(m_tree_2lSS_Fake);
	createBranch(m_tree_2lSS_Flip);
	createBranch(m_tree_2lSS1tau_SR);
	createBranch(m_tree_2lSS1tau_Fake);
	createBranch(m_tree_2lSS1tau_Flip);
	createBranch(m_tree_2l2tau_SR);
	createBranch(m_tree_2l2tau_Fake);
	createBranch(m_tree_3l_SR);
	createBranch(m_tree_3l_Fake);
	createBranch(m_tree_3l1tau_SR);
	createBranch(m_tree_3l1tau_Fake);
	createBranch(m_tree_4l_SR);
	createBranch(m_tree_4l_Fake);
	createBranch(m_tree_ttWctrl_SR);
	createBranch(m_tree_ttWctrl_Fake);
	createBranch(m_tree_ttWctrl_Flip);
	createBranch(m_tree_ttZctrl_SR);
	createBranch(m_tree_ttZctrl_Fake);
     }
}

void Sync::createBranch(TTree *tr)
{   
   tr->Branch("nEvent",&nEvent,"nEvent/I");
   tr->Branch("ls",&ls,"ls/I");
   tr->Branch("run",&run,"run/I");
   tr->Branch("n_presel_mu",&n_presel_mu,"n_presel_mu/I");
   tr->Branch("n_fakeablesel_mu",&n_fakeablesel_mu,"n_fakeablesel_mu/I");
   tr->Branch("n_mvasel_mu",&n_mvasel_mu,"n_mvasel_mu/I");
   tr->Branch("n_presel_ele",&n_presel_ele,"n_presel_ele/I");
   tr->Branch("n_fakeablesel_ele",&n_fakeablesel_ele,"n_fakeablesel_ele/I");
   tr->Branch("n_mvasel_ele",&n_mvasel_ele,"n_mvasel_ele/I");
   tr->Branch("n_presel_tau",&n_presel_tau,"n_presel_tau/I");
   tr->Branch("n_presel_jet",&n_presel_jet,"n_presel_jet/I");
   
   tr->Branch("mu1_pt",&mu1_pt,"mu1_pt/F");
   tr->Branch("mu1_conept",&mu1_conept,"mu1_conept/F");
   tr->Branch("mu1_eta",&mu1_eta,"mu1_eta/F");
   tr->Branch("mu1_phi",&mu1_phi,"mu1_phi/F");
   tr->Branch("mu1_E",&mu1_E,"mu1_E/F");
   tr->Branch("mu1_charge",&mu1_charge,"mu1_charge/I");
   tr->Branch("mu1_miniRelIso",&mu1_miniRelIso,"mu1_miniRelIso/F");
   tr->Branch("mu1_miniIsoCharged",&mu1_miniIsoCharged,"mu1_miniIsoCharged/F");
   tr->Branch("mu1_miniIsoNeutral",&mu1_miniIsoNeutral,"mu1_miniIsoNeutral/F");
   tr->Branch("mu1_PFRelIso04",&mu1_PFRelIso04,"mu1_PFRelIso04/F");
   tr->Branch("mu1_jetNDauChargedMVASel",&mu1_jetNDauChargedMVASel,"mu1_jetNDauChargedMVASel/F");
   tr->Branch("mu1_jetPtRel",&mu1_jetPtRel,"mu1_jetPtRel/F");
   tr->Branch("mu1_jetPtRatio",&mu1_jetPtRatio,"mu1_jetPtRatio/F");
   tr->Branch("mu1_jetCSV",&mu1_jetCSV,"mu1_jetCSV/F");
   tr->Branch("mu1_sip3D",&mu1_sip3D,"mu1_sip3D/F");
   tr->Branch("mu1_dxy",&mu1_dxy,"mu1_dxy/F");
   tr->Branch("mu1_dxyAbs",&mu1_dxyAbs,"mu1_dxyAbs/F");
   tr->Branch("mu1_dz",&mu1_dz,"mu1_dz/F");
   tr->Branch("mu1_segmentCompatibility",&mu1_segmentCompatibility,"mu1_segmentCompatibility/F");
   tr->Branch("mu1_leptonMVA",&mu1_leptonMVA,"mu1_leptonMVA/F");
   tr->Branch("mu1_mediumID",&mu1_mediumID,"mu1_mediumID/F");
   tr->Branch("mu1_dpt_div_pt",&mu1_dpt_div_pt,"mu1_dpt_div_pt/F");
   tr->Branch("mu1_isfakeablesel",&mu1_isfakeablesel,"mu1_isfakeablesel/I");
   tr->Branch("mu1_ismvasel",&mu1_ismvasel,"mu1_ismvasel/I");

   tr->Branch("mu2_pt",&mu2_pt,"mu2_pt/F");
   tr->Branch("mu2_conept",&mu2_conept,"mu2_conept/F");
   tr->Branch("mu2_eta",&mu2_eta,"mu2_eta/F");
   tr->Branch("mu2_phi",&mu2_phi,"mu2_phi/F");
   tr->Branch("mu2_E",&mu2_E,"mu2_E/F");
   tr->Branch("mu2_charge",&mu2_charge,"mu2_charge/I");
   tr->Branch("mu2_miniRelIso",&mu2_miniRelIso,"mu2_miniRelIso/F");
   tr->Branch("mu2_miniIsoCharged",&mu2_miniIsoCharged,"mu2_miniIsoCharged/F");
   tr->Branch("mu2_miniIsoNeutral",&mu2_miniIsoNeutral,"mu2_miniIsoNeutral/F");
   tr->Branch("mu2_PFRelIso04",&mu2_PFRelIso04,"mu2_PFRelIso04/F");
   tr->Branch("mu2_jetNDauChargedMVASel",&mu2_jetNDauChargedMVASel,"mu2_jetNDauChargedMVASel/F");
   tr->Branch("mu2_jetPtRel",&mu2_jetPtRel,"mu2_jetPtRel/F");
   tr->Branch("mu2_jetPtRatio",&mu2_jetPtRatio,"mu2_jetPtRatio/F");
   tr->Branch("mu2_jetCSV",&mu2_jetCSV,"mu2_jetCSV/F");
   tr->Branch("mu2_sip3D",&mu2_sip3D,"mu2_sip3D/F");
   tr->Branch("mu2_dxy",&mu2_dxy,"mu2_dxy/F");
   tr->Branch("mu2_dxyAbs",&mu2_dxyAbs,"mu2_dxyAbs/F");
   tr->Branch("mu2_dz",&mu2_dz,"mu2_dz/F");
   tr->Branch("mu2_segmentCompatibility",&mu2_segmentCompatibility,"mu2_segmentCompatibility/F");
   tr->Branch("mu2_leptonMVA",&mu2_leptonMVA,"mu2_leptonMVA/F");
   tr->Branch("mu2_mediumID",&mu2_mediumID,"mu2_mediumID/F");
   tr->Branch("mu2_dpt_div_pt",&mu2_dpt_div_pt,"mu2_dpt_div_pt/F");
   tr->Branch("mu2_isfakeablesel",&mu2_isfakeablesel,"mu2_isfakeablesel/I");
   tr->Branch("mu2_ismvasel",&mu2_ismvasel,"mu2_ismvasel/I");
   
   tr->Branch("ele1_pt",&ele1_pt,"ele1_pt/F");
   tr->Branch("ele1_conept",&ele1_conept,"ele1_conept/F");
   tr->Branch("ele1_eta",&ele1_eta,"ele1_eta/F");
   tr->Branch("ele1_phi",&ele1_phi,"ele1_phi/F");
   tr->Branch("ele1_E",&ele1_E,"ele1_E/F");
   tr->Branch("ele1_charge",&ele1_charge,"ele1_charge/I");
   tr->Branch("ele1_miniRelIso",&ele1_miniRelIso,"ele1_miniRelIso/F");
   tr->Branch("ele1_miniIsoCharged",&ele1_miniIsoCharged,"ele1_miniIsoCharged/F");
   tr->Branch("ele1_miniIsoNeutral",&ele1_miniIsoNeutral,"ele1_miniIsoNeutral/F");
   tr->Branch("ele1_PFRelIso04",&ele1_PFRelIso04,"ele1_PFRelIso04/F");
   tr->Branch("ele1_jetNDauChargedMVASel",&ele1_jetNDauChargedMVASel,"ele1_jetNDauChargedMVASel/F");
   tr->Branch("ele1_jetPtRel",&ele1_jetPtRel,"ele1_jetPtRel/F");
   tr->Branch("ele1_jetPtRatio",&ele1_jetPtRatio,"ele1_jetPtRatio/F");
   tr->Branch("ele1_jetCSV",&ele1_jetCSV,"ele1_jetCSV/F");
   tr->Branch("ele1_sip3D",&ele1_sip3D,"ele1_sip3D/F");
   tr->Branch("ele1_dxy",&ele1_dxy,"ele1_dxy/F");
   tr->Branch("ele1_dxyAbs",&ele1_dxyAbs,"ele1_dxyAbs/F");
   tr->Branch("ele1_dz",&ele1_dz,"ele1_dz/F");
   tr->Branch("ele1_ntMVAeleID",&ele1_ntMVAeleID,"ele1_ntMVAeleID/F");
   tr->Branch("ele1_leptonMVA",&ele1_leptonMVA,"ele1_leptonMVA/F");
   tr->Branch("ele1_isChargeConsistent",&ele1_isChargeConsistent,"ele1_isChargeConsistent/I");
   tr->Branch("ele1_passesConversionVeto",&ele1_passesConversionVeto,"ele1_passesConversionVeto/F");
   tr->Branch("ele1_nMissingHits",&ele1_nMissingHits,"ele1_nMissingHits/F");
   tr->Branch("ele1_sigmaEtaEta",&ele1_sigmaEtaEta,"ele1_sigmaEtaEta/F");
   tr->Branch("ele1_HoE",&ele1_HoE,"ele1_HoE/F");
   tr->Branch("ele1_deltaEta",&ele1_deltaEta,"ele1_deltaEta/F");
   tr->Branch("ele1_deltaPhi",&ele1_deltaPhi,"ele1_deltaPhi/F");
   tr->Branch("ele1_OoEminusOoP",&ele1_OoEminusOoP,"ele1_OoEminusOoP/F");
   tr->Branch("ele1_isfakeablesel",&ele1_isfakeablesel,"ele1_isfakeablesel/I");
   tr->Branch("ele1_ismvasel",&ele1_ismvasel,"ele1_ismvasel/I");

   tr->Branch("ele2_pt",&ele2_pt,"ele2_pt/F");
   tr->Branch("ele2_conept",&ele2_conept,"ele2_conept/F");
   tr->Branch("ele2_eta",&ele2_eta,"ele2_eta/F");
   tr->Branch("ele2_phi",&ele2_phi,"ele2_phi/F");
   tr->Branch("ele2_E",&ele2_E,"ele2_E/F");
   tr->Branch("ele2_charge",&ele2_charge,"ele2_charge/I");
   tr->Branch("ele2_miniRelIso",&ele2_miniRelIso,"ele2_miniRelIso/F");
   tr->Branch("ele2_miniIsoCharged",&ele2_miniIsoCharged,"ele2_miniIsoCharged/F");
   tr->Branch("ele2_miniIsoNeutral",&ele2_miniIsoNeutral,"ele2_miniIsoNeutral/F");
   tr->Branch("ele2_PFRelIso04",&ele2_PFRelIso04,"ele2_PFRelIso04/F");
   tr->Branch("ele2_jetNDauChargedMVASel",&ele2_jetNDauChargedMVASel,"ele2_jetNDauChargedMVASel/F");
   tr->Branch("ele2_jetPtRel",&ele2_jetPtRel,"ele2_jetPtRel/F");
   tr->Branch("ele2_jetPtRatio",&ele2_jetPtRatio,"ele2_jetPtRatio/F");
   tr->Branch("ele2_jetCSV",&ele2_jetCSV,"ele2_jetCSV/F");
   tr->Branch("ele2_sip3D",&ele2_sip3D,"ele2_sip3D/F");
   tr->Branch("ele2_dxy",&ele2_dxy,"ele2_dxy/F");
   tr->Branch("ele2_dxyAbs",&ele2_dxyAbs,"ele2_dxyAbs/F");
   tr->Branch("ele2_dz",&ele2_dz,"ele2_dz/F");
   tr->Branch("ele2_ntMVAeleID",&ele2_ntMVAeleID,"ele2_ntMVAeleID/F");
   tr->Branch("ele2_leptonMVA",&ele2_leptonMVA,"ele2_leptonMVA/F");
   tr->Branch("ele2_isChargeConsistent",&ele2_isChargeConsistent,"ele2_isChargeConsistent/I");
   tr->Branch("ele2_passesConversionVeto",&ele2_passesConversionVeto,"ele2_passesConversionVeto/F");
   tr->Branch("ele2_nMissingHits",&ele2_nMissingHits,"ele2_nMissingHits/F");
   tr->Branch("ele2_sigmaEtaEta",&ele2_sigmaEtaEta,"ele2_sigmaEtaEta/F");
   tr->Branch("ele2_HoE",&ele2_HoE,"ele2_HoE/F");
   tr->Branch("ele2_deltaEta",&ele2_deltaEta,"ele2_deltaEta/F");
   tr->Branch("ele2_deltaPhi",&ele2_deltaPhi,"ele2_deltaPhi/F");
   tr->Branch("ele2_OoEminusOoP",&ele2_OoEminusOoP,"ele2_OoEminusOoP/F");
   tr->Branch("ele2_isfakeablesel",&ele2_isfakeablesel,"ele2_isfakeablesel/I");
   tr->Branch("ele2_ismvasel",&ele2_ismvasel,"ele2_ismvasel/I");

   tr->Branch("tau1_pt",&tau1_pt,"tau1_pt/F");
   tr->Branch("tau1_eta",&tau1_eta,"tau1_eta/F");
   tr->Branch("tau1_phi",&tau1_phi,"tau1_phi/F");
   tr->Branch("tau1_E",&tau1_E,"tau1_E/F");
   tr->Branch("tau1_charge",&tau1_charge,"tau1_charge/I");
   tr->Branch("tau1_dxy",&tau1_dxy,"tau1_dxy/F");
   tr->Branch("tau1_dz",&tau1_dz,"tau1_dz/F");
   tr->Branch("tau1_decayModeFindingOldDMs",&tau1_decayModeFindingOldDMs,"tau1_decayModeFindingOldDMs/F");
   tr->Branch("tau1_decayModeFindingNewDMs",&tau1_decayModeFindingNewDMs,"tau1_decayModeFindingNewDMs/F");
   tr->Branch("tau1_byCombinedIsolationDeltaBetaCorr3Hits",&tau1_byCombinedIsolationDeltaBetaCorr3Hits,"tau1_byCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits",&tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits,"tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits",&tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits,"tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3Hits",&tau1_byTightCombinedIsolationDeltaBetaCorr3Hits,"tau1_byTightCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_rawMVArun2v1DBdR03oldDMwLT",&tau1_rawMVArun2v1DBdR03oldDMwLT,"tau1_rawMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau1_againstMuonLoose3",&tau1_againstMuonLoose3,"tau1_againstMuonLoose3/F");
   tr->Branch("tau1_againstMuonTight3",&tau1_againstMuonTight3,"tau1_againstMuonTight3/F");
   tr->Branch("tau1_againstElectronVLooseMVA6",&tau1_againstElectronVLooseMVA6,"tau1_againstElectronVLooseMVA6/F");
   tr->Branch("tau1_againstElectronLooseMVA6",&tau1_againstElectronLooseMVA6,"tau1_againstElectronLooseMVA6/F");
   tr->Branch("tau1_againstElectronMediumMVA6",&tau1_againstElectronMediumMVA6,"tau1_againstElectronMediumMVA6/F");
   tr->Branch("tau1_againstElectronTightMVA6",&tau1_againstElectronTightMVA6,"tau1_againstElectronTightMVA6/F");

   tr->Branch("tau2_pt",&tau2_pt,"tau2_pt/F");
   tr->Branch("tau2_eta",&tau2_eta,"tau2_eta/F");
   tr->Branch("tau2_phi",&tau2_phi,"tau2_phi/F");
   tr->Branch("tau2_E",&tau2_E,"tau2_E/F");
   tr->Branch("tau2_charge",&tau2_charge,"tau2_charge/I");
   tr->Branch("tau2_dxy",&tau2_dxy,"tau2_dxy/F");
   tr->Branch("tau2_dz",&tau2_dz,"tau2_dz/F");
   tr->Branch("tau2_decayModeFindingOldDMs",&tau2_decayModeFindingOldDMs,"tau2_decayModeFindingOldDMs/F");
   tr->Branch("tau2_decayModeFindingNewDMs",&tau2_decayModeFindingNewDMs,"tau2_decayModeFindingNewDMs/F");
   tr->Branch("tau2_byCombinedIsolationDeltaBetaCorr3Hits",&tau2_byCombinedIsolationDeltaBetaCorr3Hits,"tau2_byCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits",&tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits,"tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits",&tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits,"tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau2_byTightCombinedIsolationDeltaBetaCorr3Hits",&tau2_byTightCombinedIsolationDeltaBetaCorr3Hits,"tau2_byTightCombinedIsolationDeltaBetaCorr3Hits/F");
   tr->Branch("tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   tr->Branch("tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_rawMVArun2v1DBdR03oldDMwLT",&tau2_rawMVArun2v1DBdR03oldDMwLT,"tau2_rawMVArun2v1DBdR03oldDMwLT/F");
   tr->Branch("tau2_againstMuonLoose3",&tau2_againstMuonLoose3,"tau2_againstMuonLoose3/F");
   tr->Branch("tau2_againstMuonTight3",&tau2_againstMuonTight3,"tau2_againstMuonTight3/F");
   tr->Branch("tau2_againstElectronVLooseMVA6",&tau2_againstElectronVLooseMVA6,"tau2_againstElectronVLooseMVA6/F");
   tr->Branch("tau2_againstElectronLooseMVA6",&tau2_againstElectronLooseMVA6,"tau2_againstElectronLooseMVA6/F");
   tr->Branch("tau2_againstElectronMediumMVA6",&tau2_againstElectronMediumMVA6,"tau2_againstElectronMediumMVA6/F");
   tr->Branch("tau2_againstElectronTightMVA6",&tau2_againstElectronTightMVA6,"tau2_againstElectronTightMVA6/F");
   
   tr->Branch("jet1_pt",&jet1_pt,"jet1_pt/F");
   tr->Branch("jet1_eta",&jet1_eta,"jet1_eta/F");
   tr->Branch("jet1_phi",&jet1_phi,"jet1_phi/F");
   tr->Branch("jet1_E",&jet1_E,"jet1_E/F");
   tr->Branch("jet1_CSV",&jet1_CSV,"jet1_CSV/F");

   tr->Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
   tr->Branch("jet2_eta",&jet2_eta,"jet2_eta/F");
   tr->Branch("jet2_phi",&jet2_phi,"jet2_phi/F");
   tr->Branch("jet2_E",&jet2_E,"jet2_E/F");
   tr->Branch("jet2_CSV",&jet2_CSV,"jet2_CSV/F");

   tr->Branch("jet3_pt",&jet3_pt,"jet3_pt/F");
   tr->Branch("jet3_eta",&jet3_eta,"jet3_eta/F");
   tr->Branch("jet3_phi",&jet3_phi,"jet3_phi/F");
   tr->Branch("jet3_E",&jet3_E,"jet3_E/F");
   tr->Branch("jet3_CSV",&jet3_CSV,"jet3_CSV/F");

   tr->Branch("jet4_pt",&jet4_pt,"jet4_pt/F");
   tr->Branch("jet4_eta",&jet4_eta,"jet4_eta/F");
   tr->Branch("jet4_phi",&jet4_phi,"jet4_phi/F");
   tr->Branch("jet4_E",&jet4_E,"jet4_E/F");
   tr->Branch("jet4_CSV",&jet4_CSV,"jet4_CSV/F");
   
   tr->Branch("PFMET",&PFMET,"PFMET/F");
   tr->Branch("PFMETphi",&PFMETphi,"PFMETphi/F");
   tr->Branch("MHT",&MHT,"MHT/F");
   tr->Branch("metLD",&metLD,"metLD/F");
   tr->Branch("isGenMatched",&isGenMatched,"isGenMatched/I");
   
   tr->Branch("lep1_conePt",&lep1_conePt,"lep1_conePt/F");
   tr->Branch("lep2_conePt",&lep2_conePt,"lep2_conePt/F");
   tr->Branch("lep3_conePt",&lep3_conePt,"lep3_conePt/F");
   
   tr->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
   tr->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
   tr->Branch("mindr_lep3_jet",&mindr_lep3_jet,"mindr_lep3_jet/F");

   tr->Branch("mindr_tau1_jet",&mindr_tau1_jet,"mindr_tau1_jet/F");
   tr->Branch("mindr_tau2_jet",&mindr_tau2_jet,"mindr_tau2_jet/F");
   tr->Branch("mindr_tau3_jet",&mindr_tau3_jet,"mindr_tau3_jet/F");
   
   tr->Branch("avg_dr_jet",&avg_dr_jet,"avg_dr_jet/F");
   tr->Branch("avr_dr_lep_tau",&avr_dr_lep_tau,"avr_dr_lep_tau/F");
   tr->Branch("max_dr_lep_tau",&max_dr_lep_tau,"max_dr_lep_tau/F");
   tr->Branch("mindr_tau_jet",&mindr_tau_jet,"mindr_tau_jet/F");
   tr->Branch("min_dr_lep_tau",&min_dr_lep_tau,"min_dr_lep_tau/F");
   tr->Branch("min_dr_lep_jet",&min_dr_lep_jet,"min_dr_lep_jet/F");
   tr->Branch("dr_leps",&dr_leps,"dr_leps/F");
   tr->Branch("dr_taus",&dr_taus,"dr_taus/F");
   tr->Branch("dR_lep_tau_ss",&dR_lep_tau_ss,"dR_lep_tau_ss/F");
   tr->Branch("dr_lep1_tau",&dr_lep1_tau,"dr_lep1_tau/F");
   tr->Branch("dr_lep2_tau",&dr_lep2_tau,"dr_lep2_tau/F");
   tr->Branch("max_lep_eta",&max_lep_eta,"max_lep_eta/F");
   tr->Branch("mT_lep1",&mT_lep1,"mT_lep1/F");
   tr->Branch("mT_lep2",&mT_lep2,"mT_lep2/F");
   tr->Branch("mTauTauVis",&mTauTauVis,"mTauTauVis/F");
   tr->Branch("mTauTauVis1",&mTauTauVis1,"mTauTauVis1/F");
   tr->Branch("mTauTauVis2",&mTauTauVis2,"mTauTauVis2/F");
   tr->Branch("mbb",&mbb,"mbb/F");
   tr->Branch("mbb_loose",&mbb_loose,"mbb_loose/F");
   tr->Branch("cosThetaS_hadTau",&cosThetaS_hadTau,"cosThetaS_hadTau/F");
   tr->Branch("HTT",&HTT,"HTT/F");
   tr->Branch("HadTop_pt",&HadTop_pt,"HadTop_pt/F");
   tr->Branch("Hj_tagger",&Hj_tagger,"Hj_tagger/F");
   tr->Branch("nBJetLoose",&nBJetLoose,"nBJetLoose/I");
   
   tr->Branch("mvaOutput_plainKin_ttV",&mvaOutput_plainKin_ttV,"mvaOutput_plainKin_ttV/F");
   tr->Branch("mvaOutput_plainKin_ttbar",&mvaOutput_plainKin_ttbar,"mvaOutput_plainKin_ttbar/F");
   tr->Branch("mvaOutput_1l_2tau_HTT_SUM_VT",&mvaOutput_1l_2tau_HTT_SUM_VT,"mvaOutput_1l_2tau_HTT_SUM_VT/F");
   tr->Branch("mvaOutput_2l_2tau_plainKin_1B_VT",&mvaOutput_2l_2tau_plainKin_1B_VT,"mvaOutput_2l_2tau_plainKin_1B_VT/F");
   tr->Branch("mvaOutput_2l_2tau_plainKin_SUM_VT",&mvaOutput_2l_2tau_plainKin_SUM_VT,"mvaOutput_2l_2tau_plainKin_SUM_VT/F");
   tr->Branch("mvaOutput_2lss_ttV",&mvaOutput_2lss_ttV,"mvaOutput_2lss_ttV/F");
   tr->Branch("mvaOutput_2lss_ttbar",&mvaOutput_2lss_ttbar,"mvaOutput_2lss_ttbar/F");
   tr->Branch("mvaOutput_2lss_1tau_plainKin_ttbar",&mvaOutput_2lss_1tau_plainKin_ttbar,"mvaOutput_2lss_1tau_plainKin_ttbar/F");
   tr->Branch("mvaOutput_2lss_1tau_plainKin_ttV",&mvaOutput_2lss_1tau_plainKin_ttV,"mvaOutput_2lss_1tau_plainKin_ttV/F");
   tr->Branch("mvaOutput_2lss_1tau_plainKin_1B_M",&mvaOutput_2lss_1tau_plainKin_1B_M,"mvaOutput_2lss_1tau_plainKin_1B_M/F");
   tr->Branch("mvaOutput_2lss_1tau_plainKin_SUM_M",&mvaOutput_2lss_1tau_plainKin_SUM_M,"mvaOutput_2lss_1tau_plainKin_SUM_M/F");
   tr->Branch("mvaOutput_2lss_1tau_HTT_SUM_M",&mvaOutput_2lss_1tau_HTT_SUM_M,"mvaOutput_2lss_1tau_HTT_SUM_M/F");
   tr->Branch("mvaOutput_2lss_1tau_HTTMEM_SUM_M",&mvaOutput_2lss_1tau_HTTMEM_SUM_M,"mvaOutput_2lss_1tau_HTTMEM_SUM_M/F");
   tr->Branch("mvaOutput_3l_ttV",&mvaOutput_3l_ttV,"mvaOutput_3l_ttV/F");
   tr->Branch("mvaOutput_3l_ttbar",&mvaOutput_3l_ttbar,"mvaOutput_3l_ttbar/F");
   tr->Branch("mvaOutput_3l_1tau_plainKin_SUM_M",&mvaOutput_3l_1tau_plainKin_SUM_M,"mvaOutput_3l_1tau_plainKin_SUM_M/F");
   tr->Branch("mvaOutput_3l_1tau_plainKin_1B_M",&mvaOutput_3l_1tau_plainKin_1B_M,"mvaOutput_3l_1tau_plainKin_1B_M/F");
   
   tr->Branch("FR_weight",&FR_weight,"FR_weight/F");
   tr->Branch("triggerSF_weight",&triggerSF_weight,"triggerSF_weight/F");
   tr->Branch("leptonSF_weight",&leptonSF_weight,"leptonSF_weight/F");
   tr->Branch("tauSF_weight",&tauSF_weight,"tauSF_weight/F");
   tr->Branch("bTagSF_weight",&bTagSF_weight,"bTagSF_weight/F");
   tr->Branch("PU_weight",&PU_weight,"PU_weight/F");
   tr->Branch("MC_weight",&MC_weight,"MC_weight/F");
   
   tr->Branch("Integral_ttH",&Integral_ttH,"Integral_ttH/F");
   tr->Branch("Integral_ttZ",&Integral_ttZ,"Integral_ttZ/F");
   tr->Branch("Integral_ttZ_Zll",&Integral_ttZ_Zll,"Integral_ttZ_Zll/F");
   tr->Branch("Integral_ttbar",&Integral_ttbar,"Integral_ttbar/F");
   tr->Branch("integration_type",&integration_type,"integration_type/F");
   tr->Branch("memOutput_LR",&memOutput_LR,"memOutput_LR/F");
}

void Sync::initVar()
{
   nEvent = -9999;
   ls = -9999;
   run = -9999;
   n_presel_mu = -9999;
   n_fakeablesel_mu = -9999;
   n_mvasel_mu = -9999;
   n_presel_ele = -9999;
   n_fakeablesel_ele = -9999;
   n_mvasel_ele = -9999;
   n_presel_tau = -9999;
   n_presel_jet = -9999;
   
   mu1_pt = -9999;
   mu1_conept = -9999;
   mu1_eta = -9999;
   mu1_phi = -9999;
   mu1_E = -9999;
   mu1_charge = -9999;
   mu1_miniRelIso = -9999;
   mu1_miniIsoCharged = -9999;
   mu1_miniIsoNeutral = -9999;
   mu1_PFRelIso04 = -9999;
   mu1_jetNDauChargedMVASel = -9999;
   mu1_jetPtRel = -9999;
   mu1_jetPtRatio = -9999;
   mu1_jetCSV = -9999;
   mu1_sip3D = -9999;
   mu1_dxy = -9999;
   mu1_dxyAbs = -9999;
   mu1_dz = -9999;
   mu1_segmentCompatibility = -9999;
   mu1_leptonMVA = -9999;
   mu1_mediumID = -9999;
   mu1_dpt_div_pt = -9999;
   mu1_isfakeablesel = -9999;
   mu1_ismvasel = -9999;

   mu2_pt = -9999;
   mu2_conept = -9999;
   mu2_eta = -9999;
   mu2_phi = -9999;
   mu2_E = -9999;
   mu2_charge = -9999;
   mu2_miniRelIso = -9999;
   mu2_miniIsoCharged = -9999;
   mu2_miniIsoNeutral = -9999;
   mu2_PFRelIso04 = -9999;
   mu2_jetNDauChargedMVASel = -9999;
   mu2_jetPtRel = -9999;
   mu2_jetPtRatio = -9999;
   mu2_jetCSV = -9999;
   mu2_sip3D = -9999;
   mu2_dxy = -9999;
   mu2_dxyAbs = -9999;
   mu2_dz = -9999;
   mu2_segmentCompatibility = -9999;
   mu2_leptonMVA = -9999;
   mu2_mediumID = -9999;
   mu2_dpt_div_pt = -9999;
   mu2_isfakeablesel = -9999;
   mu2_ismvasel = -9999;
   
   ele1_pt = -9999;
   ele1_conept = -9999;
   ele1_eta = -9999;
   ele1_phi = -9999;
   ele1_E = -9999;
   ele1_charge = -9999;
   ele1_miniRelIso = -9999;
   ele1_miniIsoCharged = -9999;
   ele1_miniIsoNeutral = -9999;
   ele1_PFRelIso04 = -9999;
   ele1_jetNDauChargedMVASel = -9999;
   ele1_jetPtRel = -9999;
   ele1_jetPtRatio = -9999;
   ele1_jetCSV = -9999;
   ele1_sip3D = -9999;
   ele1_dxy = -9999;
   ele1_dxyAbs = -9999;
   ele1_dz = -9999;
   ele1_ntMVAeleID = -9999;
   ele1_leptonMVA = -9999;
   ele1_isChargeConsistent = -9999;
   ele1_passesConversionVeto = -9999;
   ele1_nMissingHits = -9999;
   ele1_sigmaEtaEta = -9999;
   ele1_HoE = -9999;
   ele1_deltaEta = -9999;
   ele1_deltaPhi = -9999;
   ele1_OoEminusOoP = -9999;
   ele1_isfakeablesel = -9999;
   ele1_ismvasel = -9999;

   ele2_pt = -9999;
   ele2_conept = -9999;
   ele2_eta = -9999;
   ele2_phi = -9999;
   ele2_E = -9999;
   ele2_charge = -9999;
   ele2_miniRelIso = -9999;
   ele2_miniIsoCharged = -9999;
   ele2_miniIsoNeutral = -9999;
   ele2_PFRelIso04 = -9999;
   ele2_jetNDauChargedMVASel = -9999;
   ele2_jetPtRel = -9999;
   ele2_jetPtRatio = -9999;
   ele2_jetCSV = -9999;
   ele2_sip3D = -9999;
   ele2_dxy = -9999;
   ele2_dxyAbs = -9999;
   ele2_dz = -9999;
   ele2_ntMVAeleID = -9999;
   ele2_leptonMVA = -9999;
   ele2_isChargeConsistent = -9999;
   ele2_passesConversionVeto = -9999;
   ele2_nMissingHits = -9999;
   ele2_sigmaEtaEta = -9999;
   ele2_HoE = -9999;
   ele2_deltaEta = -9999;
   ele2_deltaPhi = -9999;
   ele2_OoEminusOoP = -9999;
   ele2_isfakeablesel = -9999;
   ele2_ismvasel = -9999;

   tau1_pt = -9999;
   tau1_eta = -9999;
   tau1_phi = -9999;
   tau1_E = -9999;
   tau1_charge = -9999;
   tau1_dxy = -9999;
   tau1_dz = -9999;
   tau1_decayModeFindingOldDMs = -9999;
   tau1_decayModeFindingNewDMs = -9999;
   tau1_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_rawMVArun2v1DBdR03oldDMwLT = -9999;
   tau1_againstMuonLoose3 = -9999;
   tau1_againstMuonTight3 = -9999;
   tau1_againstElectronVLooseMVA6 = -9999;
   tau1_againstElectronLooseMVA6 = -9999;
   tau1_againstElectronMediumMVA6 = -9999;
   tau1_againstElectronTightMVA6 = -9999;
   
   tau2_pt = -9999;
   tau2_eta = -9999;
   tau2_phi = -9999;
   tau2_E = -9999;
   tau2_charge = -9999;
   tau2_dxy = -9999;
   tau2_dz = -9999;
   tau2_decayModeFindingOldDMs = -9999;
   tau2_decayModeFindingNewDMs = -9999;
   tau2_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau2_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
   tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
   tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_rawMVArun2v1DBdR03oldDMwLT = -9999;
   tau2_againstMuonLoose3 = -9999;
   tau2_againstMuonTight3 = -9999;
   tau2_againstElectronVLooseMVA6 = -9999;
   tau2_againstElectronLooseMVA6 = -9999;
   tau2_againstElectronMediumMVA6 = -9999;
   tau2_againstElectronTightMVA6 = -9999;
   
   jet1_pt = -9999;
   jet1_eta = -9999;
   jet1_phi = -9999;
   jet1_E = -9999;
   jet1_CSV = -9999;

   jet2_pt = -9999;
   jet2_eta = -9999;
   jet2_phi = -9999;
   jet2_E = -9999;
   jet2_CSV = -9999;

   jet3_pt = -9999;
   jet3_eta = -9999;
   jet3_phi = -9999;
   jet3_E = -9999;
   jet3_CSV = -9999;

   jet4_pt = -9999;
   jet4_eta = -9999;
   jet4_phi = -9999;
   jet4_E = -9999;
   jet4_CSV = -9999;
   
   PFMET = -9999;
   PFMETphi = -9999;
   MHT = -9999;
   metLD = -9999;
   isGenMatched = 1;
   
   lep1_conePt = -9999;
   lep2_conePt = -9999;
   lep3_conePt = -9999;
   
   mindr_lep1_jet = -9999;
   mindr_lep2_jet = -9999;
   mindr_lep3_jet = -9999;
   
   mindr_tau1_jet = -9999;
   mindr_tau2_jet = -9999;
   mindr_tau3_jet = -9999;
   
   avg_dr_jet = -9999;
   avr_dr_lep_tau = -9999;
   max_dr_lep_tau = -9999;
   mindr_tau_jet = -9999;
   min_dr_lep_tau = -9999;
   min_dr_lep_jet = -9999;
   dr_leps = -9999;
   dr_taus = -9999;
   dR_lep_tau_ss = -9999;
   dr_lep1_tau = -9999;
   dr_lep2_tau = -9999;
   max_lep_eta = -9999;
   mT_lep1 = -9999;
   mT_lep2 = -9999;
   mTauTauVis = -9999;
   mTauTauVis1 = -9999;
   mTauTauVis2 = -9999;
   mbb = -9999;
   mbb_loose = -9999;
   cosThetaS_hadTau = -9999;
   HTT = -9999;
   HadTop_pt = -9999;
   Hj_tagger = -9999;
   nBJetLoose = -9999;
   
   mvaOutput_plainKin_ttV = -9999;
   mvaOutput_plainKin_ttbar = -9999;
   mvaOutput_1l_2tau_HTT_SUM_VT = -9999;
   mvaOutput_2l_2tau_plainKin_1B_VT = -9999;
   mvaOutput_2l_2tau_plainKin_SUM_VT = -9999;
   mvaOutput_2lss_ttV = -9999;
   mvaOutput_2lss_ttbar = -9999;
   mvaOutput_2lss_1tau_plainKin_ttbar = -9999;
   mvaOutput_2lss_1tau_plainKin_ttV = -9999;
   mvaOutput_2lss_1tau_plainKin_1B_M = -9999;
   mvaOutput_2lss_1tau_plainKin_SUM_M = -9999;
   mvaOutput_2lss_1tau_HTT_SUM_M = -9999;
   mvaOutput_2lss_1tau_HTTMEM_SUM_M = -9999;
   mvaOutput_3l_ttV = -9999;
   mvaOutput_3l_ttbar = -9999;
   mvaOutput_3l_1tau_plainKin_SUM_M = -9999;
   mvaOutput_3l_1tau_plainKin_1B_M = -9999;
   
   FR_weight = -9999;
   triggerSF_weight = -9999;
   leptonSF_weight = -9999;
   tauSF_weight = -9999;
   bTagSF_weight = -9999;
   PU_weight = -9999;
   MC_weight = -9999;
   
   Integral_ttH = -9999;
   Integral_ttZ = -9999;
   Integral_ttZ_Zll = -9999;
   Integral_ttbar = -9999;
   integration_type = -9999;
   memOutput_LR = -9999;
}

void Sync::get(Ntuple *nt,int n_presel_el,int n_presel_mu,int n_presel_tau,int n_presel_jet)
{
   nEvent = ntP->ev_id;
   ls = ntP->ev_lumi;
   run = ntP->ev_run;
   n_presel_mu = n_presel_mu;
   n_fakeablesel_mu = -9999;
   n_mvasel_mu = -9999;
   n_presel_ele = n_presel_el;
   n_fakeablesel_ele = -9999;
   n_mvasel_ele = -9999;
   n_presel_tau = n_presel_tau;
   n_presel_jet = n_presel_jet;
   
   int nMuon = nt->NtMuonLooseExt->size();
   
   if( nMuon > 0 )
     {	
	MuonExt mu = nt->NtMuonLooseExt->at(0);
	
	mu1_pt                     = mu.pt;
	mu1_conept                 = mu.conept;
	mu1_eta                    = mu.eta;
	mu1_phi                    = mu.phi;
	mu1_E                      = mu.E;
	mu1_charge                 = mu.charge;
	mu1_miniRelIso             = mu.iso;
	mu1_miniIsoCharged         = mu.lepMVA_miniRelIsoCharged*mu1_pt;
	mu1_miniIsoNeutral         = mu.lepMVA_miniRelIsoNeutral*mu1_pt;
	mu1_PFRelIso04             = mu.isoR04;
	mu1_jetNDauChargedMVASel   = mu.lepMVA_jetNDauChargedMVASel;
	mu1_jetPtRel               = mu.lepMVA_jetPtRelv2;
	mu1_jetPtRatio             = mu.lepMVA_jetPtRatio;
//	mu1_jetCSV                 = mu.lepMVA_jetBTagCSV;
	mu1_jetCSV                 = mu.lepMVA_jetBTagDeepCSV;
	mu1_sip3D                  = fabs(mu.sip3d);
	mu1_dxy                    = mu.dxy;
	mu1_dxyAbs                 = fabs(mu.dxy);
	mu1_dz                     = mu.dz;
	mu1_segmentCompatibility   = mu.lepMVA_mvaId;
	mu1_leptonMVA              = mu.lepMVA;
	mu1_mediumID               = mu.isMedium;
	mu1_dpt_div_pt             = mu.bestTrack_ptError/mu.bestTrack_pt;
	mu1_isfakeablesel          = mu.isFakeableTTH;
	mu1_ismvasel               = mu.isTightTTH;
     }

   if( nMuon > 1 )
     {	
	MuonExt mu = nt->NtMuonLooseExt->at(1);
	
	mu2_pt                     = mu.pt;
	mu2_conept                 = mu.conept;
	mu2_eta                    = mu.eta;
	mu2_phi                    = mu.phi;
	mu2_E                      = mu.E;
	mu2_charge                 = mu.charge;
	mu2_miniRelIso             = mu.iso;
	mu2_miniIsoCharged         = mu.lepMVA_miniRelIsoCharged*mu1_pt;
	mu2_miniIsoNeutral         = mu.lepMVA_miniRelIsoNeutral*mu1_pt;
	mu2_PFRelIso04             = mu.isoR04;
	mu2_jetNDauChargedMVASel   = mu.lepMVA_jetNDauChargedMVASel;
	mu2_jetPtRel               = mu.lepMVA_jetPtRelv2;
	mu2_jetPtRatio             = mu.lepMVA_jetPtRatio;
//	mu2_jetCSV                 = mu.lepMVA_jetBTagCSV;
	mu2_jetCSV                 = mu.lepMVA_jetBTagDeepCSV;
	mu2_sip3D                  = fabs(mu.sip3d);
	mu2_dxy                    = mu.dxy;
	mu2_dxyAbs                 = fabs(mu.dxy);
	mu2_dz                     = mu.dz;
	mu2_segmentCompatibility   = mu.lepMVA_mvaId;
	mu2_leptonMVA              = mu.lepMVA;
	mu2_mediumID               = mu.isMedium;
	mu2_dpt_div_pt             = mu.bestTrack_ptError/mu.bestTrack_pt;
	mu2_isfakeablesel          = mu.isFakeableTTH;
	mu2_ismvasel               = mu.isTightTTH;
     }

   int nElectron = nt->NtElectronLooseExt->size();
   
   if( nElectron > 0 )
     {
	ElectronExt ele = nt->NtElectronLooseExt->at(0);
	
	ele1_pt                     = ele.pt;
	ele1_conept                 = ele.conept;
	ele1_eta                    = ele.eta;
	ele1_phi                    = ele.phi;
	ele1_E                      = ele.E;
	ele1_charge                 = ele.charge;
	ele1_miniRelIso             = ele.miniIso;
	ele1_miniIsoCharged         = ele.lepMVA_miniRelIsoCharged*ele1_pt;
	ele1_miniIsoNeutral         = ele.lepMVA_miniRelIsoNeutral*ele1_pt;
	ele1_PFRelIso04             = ele.isoR04;
	ele1_jetNDauChargedMVASel   = ele.lepMVA_jetNDauChargedMVASel;
	ele1_jetPtRel               = ele.lepMVA_jetPtRelv2;
	ele1_jetPtRatio             = ele.lepMVA_jetPtRatio;
//	ele1_jetCSV                 = ele.lepMVA_jetBTagCSV;
	ele1_jetCSV                 = ele.lepMVA_jetBTagDeepCSV;
	ele1_sip3D                  = fabs(ele.sip3d);
	ele1_dxy                    = ele.dxy;
	ele1_dxyAbs                 = fabs(ele.dxy);
	ele1_dz                     = ele.dz;
	ele1_ntMVAeleID             = ele.lepMVA_mvaId;
	ele1_leptonMVA              = ele.lepMVA;
	ele1_isChargeConsistent     = ele.isGsfCtfScPixChargeConsistent;
	ele1_passesConversionVeto   = ele.passCV;
	ele1_nMissingHits           = ele.nlosthits;
	ele1_sigmaEtaEta            = ele.sigmaIetaIeta;
	ele1_HoE                    = ele.hadronicOverEm;
	ele1_deltaEta               = ele.deltaEtaSuperClusterTrackAtVtx;
	ele1_deltaPhi               = ele.deltaPhiSuperClusterTrackAtVtx;
	ele1_OoEminusOoP            = ele.ooEmooP;
	ele1_isfakeablesel          = ele.isFakeableTTH;
	ele1_ismvasel               = ele.isTightTTH;
     }

   if( nElectron > 1 )
     {	   
	ElectronExt ele = nt->NtElectronLooseExt->at(1);
	
	ele2_pt                     = ele.pt;
	ele2_conept                 = ele.conept;
	ele2_eta                    = ele.eta;
	ele2_phi                    = ele.phi;
	ele2_E                      = ele.E;
	ele2_charge                 = ele.charge;
	ele2_miniRelIso             = ele.miniIso;
	ele2_miniIsoCharged         = ele.lepMVA_miniRelIsoCharged*ele1_pt;
	ele2_miniIsoNeutral         = ele.lepMVA_miniRelIsoNeutral*ele1_pt;
	ele2_PFRelIso04             = ele.isoR04;
	ele2_jetNDauChargedMVASel   = ele.lepMVA_jetNDauChargedMVASel;
	ele2_jetPtRel               = ele.lepMVA_jetPtRelv2;
	ele2_jetPtRatio             = ele.lepMVA_jetPtRatio;
//	ele2_jetCSV                 = ele.lepMVA_jetBTagCSV;
	ele2_jetCSV                 = ele.lepMVA_jetBTagDeepCSV;
	ele2_sip3D                  = fabs(ele.sip3d);
	ele2_dxy                    = ele.dxy;
	ele2_dxyAbs                 = fabs(ele.dxy);
	ele2_dz                     = ele.dz;
	ele2_ntMVAeleID             = ele.lepMVA_mvaId;
	ele2_leptonMVA              = ele.lepMVA;
	ele2_isChargeConsistent     = ele.isGsfCtfScPixChargeConsistent;
	ele2_passesConversionVeto   = ele.passCV;
	ele2_nMissingHits           = ele.nlosthits;
	ele2_sigmaEtaEta            = ele.sigmaIetaIeta;
	ele2_HoE                    = ele.hadronicOverEm;
	ele2_deltaEta               = ele.deltaEtaSuperClusterTrackAtVtx;
	ele2_deltaPhi               = ele.deltaPhiSuperClusterTrackAtVtx;
	ele2_OoEminusOoP            = ele.ooEmooP;
	ele2_isfakeablesel          = ele.isFakeableTTH;
	ele2_ismvasel               = ele.isTightTTH;
     }

   int nTau = nt->NtTauFakeableExt->size();

   if( nTau > 0 )
     {
	TauExt tau = nt->NtTauFakeableExt->at(0);

	tau1_pt                                               = tau.pt;
	tau1_eta                                              = tau.eta;
	tau1_phi                                              = tau.phi;
	tau1_E                                                = tau.E;
	tau1_charge                                           = tau.charge;
	tau1_dxy                                              = tau.dxy;
	tau1_dz                                               = tau.dz;
	tau1_decayModeFindingOldDMs                           = tau.decayModeFinding;
	tau1_decayModeFindingNewDMs                           = tau.decayModeFindingNewDMs;
	tau1_byCombinedIsolationDeltaBetaCorr3Hits            = -9999;
	tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits       = tau.byLooseCombinedIsolationDeltaBetaCorr3Hits;
	tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits      = tau.byMediumCombinedIsolationDeltaBetaCorr3Hits;
	tau1_byTightCombinedIsolationDeltaBetaCorr3Hits       = tau.byTightCombinedIsolationDeltaBetaCorr3Hits;
	tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03   = -9999;
	tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03  = -9999;
	tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03   = -9999;
	tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT         = tau.byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
	tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT          = tau.byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT         = tau.byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT          = tau.byTightIsolationMVArun2v1DBdR03oldDMwLT;
	tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT         = tau.byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	tau1_rawMVArun2v1DBdR03oldDMwLT                       = -9999;
	tau1_againstMuonLoose3                                = tau.againstMuonLoose3;
	tau1_againstMuonTight3                                = tau.againstMuonTight3;
	tau1_againstElectronVLooseMVA6                        = tau.againstElectronVLooseMVA6;
	tau1_againstElectronLooseMVA6                         = tau.againstElectronLooseMVA6;
	tau1_againstElectronMediumMVA6                        = tau.againstElectronMediumMVA6;
	tau1_againstElectronTightMVA6                         = tau.againstElectronTightMVA6;
     }

   if( nTau > 1 )
     {	      
	TauExt tau = nt->NtTauFakeableExt->at(1);

	tau2_pt                                               = tau.pt;
	tau2_eta                                              = tau.eta;
	tau2_phi                                              = tau.phi;
	tau2_E                                                = tau.E;
	tau2_charge                                           = tau.charge;
	tau2_dxy                                              = tau.dxy;
	tau2_dz                                               = tau.dz;
	tau2_decayModeFindingOldDMs                           = tau.decayModeFinding;
	tau2_decayModeFindingNewDMs                           = tau.decayModeFindingNewDMs;
	tau2_byCombinedIsolationDeltaBetaCorr3Hits            = -9999;
	tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits       = tau.byLooseCombinedIsolationDeltaBetaCorr3Hits;
	tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits      = tau.byMediumCombinedIsolationDeltaBetaCorr3Hits;
	tau2_byTightCombinedIsolationDeltaBetaCorr3Hits       = tau.byTightCombinedIsolationDeltaBetaCorr3Hits;
	tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03   = -9999;
	tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03  = -9999;
	tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03   = -9999;
	tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT         = tau.byVLooseIsolationMVArun2v1DBdR03oldDMwLT;
	tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT          = tau.byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT         = tau.byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT          = tau.byTightIsolationMVArun2v1DBdR03oldDMwLT;
	tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT         = tau.byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	tau2_rawMVArun2v1DBdR03oldDMwLT                       = -9999;
	tau2_againstMuonLoose3                                = tau.againstMuonLoose3;
	tau2_againstMuonTight3                                = tau.againstMuonTight3;
	tau2_againstElectronVLooseMVA6                        = tau.againstElectronVLooseMVA6;
	tau2_againstElectronLooseMVA6                         = tau.againstElectronLooseMVA6;
	tau2_againstElectronMediumMVA6                        = tau.againstElectronMediumMVA6;
	tau2_againstElectronTightMVA6                         = tau.againstElectronTightMVA6;
     }

   int nJet = nt->NtJetLooseExt->size();
   
   if( nJet > 0 )
     {
	JetExt jet = nt->NtJetLooseExt->at(0);
	
	jet1_pt    = jet.pt;
	jet1_eta   = jet.eta;
	jet1_phi   = jet.phi;
	jet1_E     = jet.E;
	jet1_CSV   = jet.deepCSVb+jet.deepCSVbb;
     }   

   if( nJet > 1 )
     {
	JetExt jet = nt->NtJetLooseExt->at(1);

	jet2_pt    = jet.pt;
	jet2_eta   = jet.eta;
	jet2_phi   = jet.phi;
	jet2_E     = jet.E;
	jet2_CSV   = jet.deepCSVb+jet.deepCSVbb;
     }   

   if( nJet > 2 )
     {
	JetExt jet = nt->NtJetLooseExt->at(2);
	
	jet3_pt    = jet.pt;
	jet3_eta   = jet.eta;
	jet3_phi   = jet.phi;
	jet3_E     = jet.E;
	jet3_CSV   = jet.deepCSVb+jet.deepCSVbb;
     }   

   if( nJet > 3 )
     {
	JetExt jet = nt->NtJetLooseExt->at(3);
	
	jet4_pt    = jet.pt;
	jet4_eta   = jet.eta;
	jet4_phi   = jet.phi;
	jet4_E     = jet.E;
	jet4_CSV   = jet.deepCSVb+jet.deepCSVbb;
     }
   
   TLorentzVector jet;
   float jet_px = 0;
   float jet_py = 0;
   for(int ij=0;ij<nJet;ij++)
     {	     
	jet.SetPtEtaPhiE(nt->NtJetLooseExt->at(ij).pt,
			 nt->NtJetLooseExt->at(ij).eta,
			 nt->NtJetLooseExt->at(ij).phi,
			 nt->NtJetLooseExt->at(ij).E);
	jet_px += jet.Px();
	jet_py += jet.Py();
     }
   for(int ij=0;ij<nt->NtElectronFakeableExt->size();ij++)
     {	     
	jet.SetPtEtaPhiE(nt->NtElectronFakeableExt->at(ij).pt,
			 nt->NtElectronFakeableExt->at(ij).eta,
			 nt->NtElectronFakeableExt->at(ij).phi,
			 nt->NtElectronFakeableExt->at(ij).E);
	jet_px += jet.Px();
	jet_py += jet.Py();
	
	if( !nt->NtElectronFakeableExt->at(ij).hasMCMatch )
	  isGenMatched = 0;
     }
   for(int ij=0;ij<nt->NtMuonFakeableExt->size();ij++)
     {
	jet.SetPtEtaPhiE(nt->NtMuonFakeableExt->at(ij).pt,
			 nt->NtMuonFakeableExt->at(ij).eta,
			 nt->NtMuonFakeableExt->at(ij).phi,
			 nt->NtMuonFakeableExt->at(ij).E);
	jet_px += jet.Px();
	jet_py += jet.Py();

	if( !nt->NtMuonFakeableExt->at(ij).hasMCMatch )
	  isGenMatched = 0;
     }
   for(int ij=0;ij<nTau;ij++)
     {	     
	jet.SetPtEtaPhiE(nt->NtTauFakeableExt->at(ij).pt,
			 nt->NtTauFakeableExt->at(ij).eta,
			 nt->NtTauFakeableExt->at(ij).phi,
			 nt->NtTauFakeableExt->at(ij).E);
	jet_px += jet.Px();
	jet_py += jet.Py();
	
	if( !nt->NtTauFakeableExt->at(ij).hasMCMatch )
	  isGenMatched = 0;
     }
   MHT = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );
   metLD = ntP->met_pt*0.00397 + MHT*0.00265;
   PFMET = ntP->met_pt;
   PFMETphi = ntP->met_phi;
}

bool Sync::fill(Ntuple *nt,EventExt *ev)
{
   bool pass = 1;
   
   if( sync == 1 ) m_tree->Fill();
   else
     {
	bool trig_e = ev->trig_e;
	bool trig_etau = ev->trig_etau;
	bool trig_m = ev->trig_m;
	bool trig_mtau = ev->trig_mtau;
	
	bool trig_ee = ev->trig_ee;
	bool trig_em = ev->trig_em;
	bool trig_mm = ev->trig_mm;
	
	bool trig_eee = ev->trig_eee;
	bool trig_eem = ev->trig_eem;
	bool trig_emm = ev->trig_emm;
	bool trig_mmm = ev->trig_mmm;

	std::vector<Base> *elmuLoose = new std::vector<Base>;
	for(int ie=0;ie<nt->NtElectronLooseExt->size();ie++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtElectronLooseExt->at(ie)));
	     lep->iElec = ie;
	     elmuLoose->push_back(*lep);
	  }
	for(int im=0;im<nt->NtMuonLooseExt->size();im++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtMuonLooseExt->at(im)));
	     lep->iMuon = im;
	     elmuLoose->push_back(*lep);
	  }
	std::sort(elmuLoose->begin(),elmuLoose->end(),sort_by_conept());

	std::vector<Base> *elmuTight = new std::vector<Base>;
	for(int ie=0;ie<nt->NtElectronTightExt->size();ie++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtElectronTightExt->at(ie)));
	     lep->iElec = ie;
	     elmuTight->push_back(*lep);
	  }
	for(int im=0;im<nt->NtMuonTightExt->size();im++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtMuonTightExt->at(im)));
	     lep->iMuon = im;
	     elmuTight->push_back(*lep);
	  }
	std::sort(elmuTight->begin(),elmuTight->end(),sort_by_conept());

	std::vector<Base> *elmuFakeable = new std::vector<Base>;
	int nElectronFakeable = nt->NtElectronFakeableExt->size();
	for(int ie=0;ie<nElectronFakeable;ie++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtElectronFakeableExt->at(ie)));
	     lep->iElec = ie;
	     elmuFakeable->push_back(*lep);
	  }
	int nMuonFakeable = nt->NtMuonFakeableExt->size();
	for(int im=0;im<nMuonFakeable;im++ )
	  {
	     Base* lep = dynamic_cast<Base*>(&(nt->NtMuonFakeableExt->at(im)));
	     lep->iMuon = im;
	     elmuFakeable->push_back(*lep);
	  }
	std::sort(elmuFakeable->begin(),elmuFakeable->end(),sort_by_conept());
	
	int nLepLoose = elmuLoose->size();
	int nLepFakeable = elmuFakeable->size();
	int nLepTight = elmuTight->size();
	int nTauFakeable = nt->NtTauFakeableExt->size();
	int nTauLoose = nt->NtTauLooseExt->size();
	int nTauMedium = nt->NtTauMediumExt->size();
	
	int nJetLoose = nt->NtJetLooseExt->size();
	int nJetLooseBL = 0;
	int nJetLooseBM = 0;
	for(int ij=0;ij<nJetLoose;ij++)
	  {
	     if( nt->NtJetLooseExt->at(ij).isMediumBTag ) nJetLooseBM++;
	     if( nt->NtJetLooseExt->at(ij).isLooseBTag ) nJetLooseBL++;
	  }

	int nSFOS = 0;
	float mll_min = 10E+10;
	float mll_z_min = 10E+10;
	float mllll = -1;
	for(int il=0;il<nLepLoose;il++)
	  {
	     TLorentzVector *l1 = new TLorentzVector();
	     l1->SetPtEtaPhiE(elmuLoose->at(il).pt,
			      elmuLoose->at(il).eta,
			      elmuLoose->at(il).phi,
			      elmuLoose->at(il).E);

	     for(int ill=il+1;ill<nLepLoose;ill++)
	       {
		  TLorentzVector *l2 = new TLorentzVector();
		  l2->SetPtEtaPhiE(elmuLoose->at(ill).pt,
				   elmuLoose->at(ill).eta,
				   elmuLoose->at(ill).phi,
				   elmuLoose->at(ill).E);

		  float mll = (*l1+*l2).M();
		  float mll_z = fabs((*l1+*l2).M()-91.2);

		  bool SFOS = (((elmuLoose->at(il).iElec >= 0 && elmuLoose->at(ill).iElec >= 0) ||
				(elmuLoose->at(il).iMuon >= 0 && elmuLoose->at(ill).iMuon >= 0)) &&
			       elmuLoose->at(il).charge*elmuLoose->at(ill).charge < 0);
		  if( SFOS ) nSFOS++;
		  
		  if( mll < mll_min ) mll_min = mll;
		  if( mll_z < mll_z_min && SFOS ) mll_z_min = mll_z;
		  
		  delete l2;

		  for(int illl=ill+1;illl<nLepLoose;illl++)
		    {
		       TLorentzVector *l3 = new TLorentzVector();
		       l3->SetPtEtaPhiE(elmuLoose->at(illl).pt,
					elmuLoose->at(illl).eta,
					elmuLoose->at(illl).phi,
					elmuLoose->at(illl).E);
		       
		       for(int illll=illl+1;illll<nLepLoose;illll++)
			 {
			    TLorentzVector *l4 = new TLorentzVector();
			    l4->SetPtEtaPhiE(elmuLoose->at(illll).pt,
					     elmuLoose->at(illll).eta,
					     elmuLoose->at(illll).phi,
					     elmuLoose->at(illll).E);

			    bool SFOS34 = (((elmuLoose->at(illl).iElec >= 0 && elmuLoose->at(illll).iElec >= 0) ||
					    (elmuLoose->at(illl).iMuon >= 0 && elmuLoose->at(illll).iMuon >= 0)) &&
					   elmuLoose->at(illl).charge*elmuLoose->at(illll).charge < 0);
			    
			    if( !SFOS34 ) continue;
			    
			    float mllll = (*l1+*l2+*l3+*l4).M();

			    delete l4;
			 }
		       
		       delete l3;
		    }		  
	       }	    
	     
	     delete l1;
	  }		
	bool pass_mll = (mll_min > 12.);
	bool pass_mll_z = (mll_z_min > 10.);

	float mee_z_min = 10E+10;
	for(int ie=0;ie<nt->NtElectronLooseExt->size();ie++)
	  {
	     TLorentzVector *e1 = new TLorentzVector();
	     e1->SetPtEtaPhiE(nt->NtElectronLooseExt->at(ie).pt,
			      nt->NtElectronLooseExt->at(ie).eta,
			      nt->NtElectronLooseExt->at(ie).phi,
			      nt->NtElectronLooseExt->at(ie).E);

	     for(int iee=ie+1;iee<nt->NtElectronLooseExt->size();iee++)
	       {
		  TLorentzVector *e2 = new TLorentzVector();
		  e2->SetPtEtaPhiE(nt->NtElectronLooseExt->at(iee).pt,
				   nt->NtElectronLooseExt->at(iee).eta,
				   nt->NtElectronLooseExt->at(iee).phi,
				   nt->NtElectronLooseExt->at(iee).E);		  

		  float mee = fabs((*e1+*e2).M()-91.2);
		  
		  if( mee < mee_z_min ) mee_z_min = mee;
		  
		  delete e2;
	       }	    
	     
	     delete e1;
	  }		
	bool pass_mee = (mee_z_min > 10.);

	bool pass_1l2tau_SR = 0;
	bool pass_1l2tau_Fake = 0;
	  {	     
	     bool pass_nlep = (nLepFakeable > 0 && nLepTight <= 1 && nTauFakeable >= 2);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && (trig_e || trig_etau)) ||
		    (nMuonFakeable >= 1 && (trig_m || trig_mtau));
		  bool pass_fakeable_pt = (elmuFakeable->at(0).iMuon >= 0 && elmuFakeable->at(0).conept > 25.) ||
		    (elmuFakeable->at(0).iElec >= 0 && elmuFakeable->at(0).conept > 30.);
		  bool pass_fakeable_eta = (fabs(elmuFakeable->at(0).eta) < 2.1);
		  bool pass_tau_pt = (nt->NtTauFakeableExt->at(0).pt > 30. && nt->NtTauFakeableExt->at(1).pt > 20.);
		  bool pass_njet = (nJetLoose >= 3 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH);
		  bool pass_tau_id = (nt->NtTauFakeableExt->at(0).isMediumTTH &&
				      nt->NtTauFakeableExt->at(1).isMediumTTH);
		  bool pass_tau_charge = (nt->NtTauFakeableExt->at(0).charge*nt->NtTauFakeableExt->at(1).charge < 0);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch &&
				     nt->NtTauFakeableExt->at(0).hasMCMatch && nt->NtTauFakeableExt->at(1).hasMCMatch);
		  
		  bool pass_fake_os = (nt->NtTauFakeableExt->at(0).charge*nt->NtTauFakeableExt->at(1).charge < 0);
		  
		  pass_1l2tau_SR = (pass_trig && pass_fakeable_pt && pass_fakeable_eta && pass_tau_pt && pass_njet && pass_mll &&
				    pass_tight && pass_tau_id && pass_tau_charge && pass_truth && pass_fake_os);

		  pass_1l2tau_Fake = (pass_trig && pass_fakeable_pt && pass_fakeable_eta && pass_tau_pt && pass_njet && pass_mll &&
				      (!pass_tight || !pass_tau_id) && pass_fake_os);
		  
		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 1l2tau" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_fakeable_eta = " << pass_fakeable_eta << std::endl;
			    std::cout << "  pass_tau_pt = " << pass_tau_pt << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_tau_id = " << pass_tau_id << std::endl;
			    std::cout << "  && pass_tau_charge = " << pass_tau_charge << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_1l2tau_SR = " << pass_1l2tau_SR << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  && pass_fake_os = " << pass_fake_os << std::endl;
			    std::cout << "  = pass_1l2tau_Fake = " << pass_1l2tau_Fake << std::endl;
			 }		  
		    }		  
	       }
	  }

	bool pass_2lSS_SR = 0;
	bool pass_2lSS_Fake = 0;
	bool pass_2lSS_Flip = 0;	
	bool pass_ttWctrl_SR = 0;
	bool pass_ttWctrl_Fake = 0;
	bool pass_ttWctrl_Flip = 0;
	  {	     
	     bool pass_nlep = (nLepFakeable > 1 && nLepTight <= 2);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && trig_e) ||
		    (nMuonFakeable >= 1 && trig_m) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 25. && elmuFakeable->at(1).conept > 15.);
		  bool pass_tight_charge = (elmuFakeable->at(0).tightCharge && elmuFakeable->at(1).tightCharge);
		  bool pass_tau_veto = (nTauLoose == 0);
		  bool pass_njet = (nJetLoose >= 4 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  bool pass_njet3 = (nJetLoose == 3 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  bool pass_metLD = (metLD > 0.2);
		  if( !(elmuFakeable->at(0).iElec >= 0) || !(elmuFakeable->at(1).iElec >= 0) ) pass_metLD = 1;
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH);
		  bool pass_ss = (elmuFakeable->at(0).charge*elmuFakeable->at(1).charge > 0);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch);

		  pass_2lSS_SR = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet &&
				  pass_tight && pass_ss && pass_truth);

		  pass_ttWctrl_SR = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet3 &&
				     pass_tight && pass_ss && pass_truth);

		  pass_2lSS_Fake = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet &&
				    pass_ss && !pass_tight);

		  pass_ttWctrl_Fake = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet3 &&
				       pass_ss && !pass_tight);

		  bool pass_has_elec = (elmuFakeable->at(0).iElec >= 0 || elmuFakeable->at(1).iElec >= 0);
		  
		  pass_2lSS_Flip = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet &&
				    pass_has_elec && !pass_ss && pass_tight && pass_truth);

		  pass_ttWctrl_Flip = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_veto && pass_mee && pass_metLD && pass_njet3 &&
				       pass_has_elec && !pass_ss && pass_tight && pass_truth);

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 2lSS" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_tight_charge = " << pass_tight_charge << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_tau_veto = " << pass_tau_veto << std::endl;
			    std::cout << "  pass_mee = " << pass_mee << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_2lSS_SR = " << pass_2lSS_SR << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  = pass_2lSS_Fake = " << pass_2lSS_Fake << std::endl;
			    std::cout << "  && pass_has_elec = " << pass_has_elec << std::endl;
			    std::cout << "  && !pass_ss = " << !pass_ss << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_2lSS_Flip = " << pass_2lSS_Flip << std::endl;
			    std::cout << "Category ttWctrl" << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_tight_charge = " << pass_tight_charge << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_tau_veto = " << pass_tau_veto << std::endl;
			    std::cout << "  pass_mee = " << pass_mee << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_njet3 = " << pass_njet3 << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_ttWctrl_SR = " << pass_ttWctrl_SR << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  = pass_ttWctrl_Fake = " << pass_ttWctrl_Fake << std::endl;
			    std::cout << "  && pass_has_elec = " << pass_has_elec << std::endl;
			    std::cout << "  && !pass_ss = " << !pass_ss << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_ttWctrl_Flip = " << pass_ttWctrl_Flip << std::endl;
			 }		  
		    }		  
	       }
	  }

	bool pass_2lSS1tau_SR = 0;
	bool pass_2lSS1tau_Fake = 0;
	bool pass_2lSS1tau_Flip = 0;
	  {
	     bool pass_nlep = (nLepFakeable > 1 && nLepTight <= 2 && nTauFakeable >= 1);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && (trig_e || trig_etau)) ||
		    (nMuonFakeable >= 1 && (trig_m || trig_mtau)) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 25.);
		  pass_fakeable_pt = pass_fakeable_pt && ((elmuLoose->at(1).iMuon >= 0) ? (elmuFakeable->at(1).conept > 10.) : (elmuFakeable->at(1).conept > 15.));
		  bool pass_tight_charge = (elmuFakeable->at(0).tightCharge && elmuFakeable->at(1).tightCharge);
		  bool pass_tau = (nTauLoose >= 1);
		  bool pass_njet = (nJetLoose >= 3 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  bool pass_metLD = (metLD > 0.2);
		  if( !(elmuFakeable->at(0).iElec >= 0) || !(elmuFakeable->at(1).iElec >= 0) ) pass_metLD = 1;
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH);
		  bool pass_ss = (elmuFakeable->at(0).charge*elmuFakeable->at(1).charge > 0);
		  bool pass_tau_tight = (nTauMedium <= 1);
		  bool pass_os_tau = (elmuFakeable->at(0).charge*nt->NtTauFakeableExt->at(0).charge < 0) &&
		    (elmuFakeable->at(1).charge*nt->NtTauFakeableExt->at(0).charge < 0);
		  bool pass_ss_tau = (elmuFakeable->at(0).iElec >= 0 && elmuFakeable->at(0).charge*nt->NtTauFakeableExt->at(0).charge > 0) ||
		    (elmuFakeable->at(1).iElec >= 0 && elmuFakeable->at(1).charge*nt->NtTauFakeableExt->at(0).charge > 0);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch);

		  pass_2lSS1tau_SR = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_tight && pass_mee && pass_metLD && pass_njet &&
				      pass_tight && pass_ss && pass_tau && pass_os_tau && pass_truth);

		  pass_2lSS1tau_Fake = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_tight && pass_mee && pass_metLD && pass_njet &&
					pass_ss && !pass_tight && pass_tau && pass_os_tau);
		  
		  bool pass_has_elec = (elmuFakeable->at(0).iElec >= 0 || elmuFakeable->at(1).iElec >= 0);
		  
		  pass_2lSS1tau_Flip = (pass_trig && pass_fakeable_pt && pass_tight_charge && pass_mll && pass_tau_tight && pass_mee && pass_metLD && pass_njet &&
					pass_has_elec && !pass_ss && pass_tight && pass_tau && pass_ss_tau && pass_truth);

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 2lSS1tau" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_tight_charge = " << pass_tight_charge << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_tau_tight = " << pass_tau_tight << std::endl;
			    std::cout << "  pass_mee = " << pass_mee << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && pass_tau = " << pass_tau << std::endl;
			    std::cout << "  && pass_os_tau = " << pass_os_tau << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_2lSS1tau_SR = " << pass_2lSS1tau_SR << std::endl;
			    std::cout << "  && pass_ss = " << pass_ss << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  && pass_tau = " << pass_tau << std::endl;
			    std::cout << "  && pass_os_tau = " << pass_os_tau << std::endl;
			    std::cout << "  = pass_2lSS1tau_Fake = " << pass_2lSS1tau_Fake << std::endl;
			    std::cout << "  && pass_has_elec = " << pass_has_elec << std::endl;
			    std::cout << "  && !pass_ss = " << !pass_ss << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_tau = " << pass_tau << std::endl;
			    std::cout << "  && pass_ss_tau = " << pass_ss_tau << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_2lSS1tau_Flip = " << pass_2lSS1tau_Flip << std::endl;
			 }		  
		    }		  
	       }
	  }

	bool pass_2l2tau_SR = 0;
	bool pass_2l2tau_Fake = 0;
	  {
	     bool pass_nlep = (nLepFakeable >= 2 && nTauFakeable >= 2);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && (trig_e || trig_etau)) ||
		    (nMuonFakeable >= 1 && (trig_m || trig_mtau)) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 25. && elmuFakeable->at(1).conept > 10.);
		  if( elmuFakeable->at(1).iElec >= 0 ) pass_fakeable_pt = pass_fakeable_pt && (elmuFakeable->at(1).conept > 15.);
		  bool pass_metLD = (((metLD > 0.2 && nSFOS == 0) || (metLD > 0.3 && nSFOS > 0)) || nJetLoose >= 4);
		  int lep_charge_sum = elmuFakeable->at(0).charge+elmuFakeable->at(1).charge+nt->NtTauFakeableExt->at(0).charge+nt->NtTauFakeableExt->at(1).charge;
		  bool pass_charge = (lep_charge_sum == 0);
		  bool pass_njet = (nJetLoose >= 2 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH);
		  bool pass_tau_tight = (nt->NtTauFakeableExt->at(0).isMediumTTH && nt->NtTauFakeableExt->at(1).isMediumTTH);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch && nt->NtTauFakeableExt->at(0).hasMCMatch && nt->NtTauFakeableExt->at(1).hasMCMatch);

		  pass_2l2tau_SR = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				    pass_tight && pass_tau_tight && pass_truth);

		  pass_2l2tau_Fake = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				      (!pass_tight || !pass_tau_tight));

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 2l2tau" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_mll_z = " << pass_mll_z << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_charge = " << pass_charge << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_tau_tight = " << pass_tau_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_2l2tau_SR = " << pass_2l2tau_SR << std::endl;
			    std::cout << "  && (!pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  || pass_tau_tight) = " << pass_tau_tight << std::endl;
			    std::cout << "  = pass_2l2tau_Fake = " << pass_2l2tau_Fake << std::endl;
			 }		  
		    }
	       }	     
	  }
	
	bool pass_3l_SR = 0;
	bool pass_3l_Fake = 0;
	bool pass_ttZctrl_SR = 0;
	bool pass_ttZctrl_Fake = 0;
	  {
	     bool pass_nlep = (nLepFakeable > 2);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && trig_e) ||
		    (nMuonFakeable >= 1 && trig_m) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em) ||
		    (nElectronFakeable >= 3 && trig_eee) ||
		    (nElectronFakeable >= 2 && nMuonFakeable >= 1 && trig_eem) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 2 && trig_emm) ||
		    (nMuonFakeable >= 3 && trig_mmm);
		  bool pass_nlep = (nLepTight <= 3);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 25. && elmuFakeable->at(1).conept > 15. && elmuFakeable->at(2).conept > 10.);
		  bool pass_tau_veto = (nTauLoose == 0);
		  bool pass_metLD = (((metLD > 0.2 && nSFOS == 0) || (metLD > 0.3 && nSFOS > 0)) || nJetLoose >= 4);
		  int lep_charge_sum = elmuFakeable->at(0).charge+elmuFakeable->at(1).charge+elmuFakeable->at(2).charge;
		  bool pass_charge = (lep_charge_sum == -1 || lep_charge_sum == 1);
		  bool pass_njet = (nJetLoose >= 2 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  bool pass_mllll = (mllll > 140 || mllll < 0);
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH && elmuFakeable->at(2).isTightTTH);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch && elmuFakeable->at(2).hasMCMatch);

		  pass_3l_SR = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_tau_veto && pass_metLD && pass_charge && pass_njet && pass_mllll &&
				pass_tight && pass_truth);

		  pass_ttZctrl_SR = (pass_trig && pass_fakeable_pt && pass_mll && !pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				     pass_tight && pass_truth);
		  
		  pass_3l_Fake = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_tau_veto && pass_metLD && pass_charge && pass_njet && pass_mllll &&
				  !pass_tight);

		  pass_ttZctrl_Fake = (pass_trig && pass_fakeable_pt && pass_mll && !pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				       !pass_tight);

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 3l" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_mll_z = " << pass_mll_z << std::endl;
			    std::cout << "  pass_tau_veto = " << pass_tau_veto << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_charge = " << pass_charge << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  pass_mllll = " << pass_mllll << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_3l_SR = " << pass_3l_SR << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  = pass_3l_Fake = " << pass_3l_Fake << std::endl;
			    std::cout << "Category ttZctrl" << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  !pass_mll_z = " << !pass_mll_z << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_charge = " << pass_charge << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_ttZctrl_SR = " << pass_ttZctrl_SR << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  = pass_ttZctrl_Fake = " << pass_ttZctrl_Fake << std::endl;
			 }		  
		    }		  
	       }	     
	  }	

	bool pass_3l1tau_SR = 0;
	bool pass_3l1tau_Fake = 0;
	  {
	     bool pass_nlep = (nLepFakeable > 2 && nTauFakeable >= 1);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && (trig_e || trig_etau)) ||
		    (nMuonFakeable >= 1 && (trig_m || trig_mtau)) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em) ||
		    (nElectronFakeable >= 3 && trig_eee) ||
		    (nElectronFakeable >= 2 && nMuonFakeable >= 1 && trig_eem) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 2 && trig_emm) ||
		    (nMuonFakeable >= 3 && trig_mmm);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 20. && elmuFakeable->at(1).conept > 10. && elmuFakeable->at(2).conept > 10.);
		  bool pass_metLD = (((metLD > 0.2 && nSFOS == 0) || (metLD > 0.3 && nSFOS > 0)) || nJetLoose >= 4);
		  int lep_charge_sum = elmuFakeable->at(0).charge+elmuFakeable->at(1).charge+elmuFakeable->at(2).charge+nt->NtTauFakeableExt->at(0).charge;
		  bool pass_charge = (lep_charge_sum == 0);
		  bool pass_njet = (nJetLoose >= 2 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  
		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH && elmuFakeable->at(2).isTightTTH);
		  bool pass_tau_tight = (nTauLoose >= 1);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch && elmuFakeable->at(2).hasMCMatch);

		  pass_3l1tau_SR = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				    pass_tight && pass_tau_tight && pass_truth);

		  pass_3l1tau_Fake = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet &&
				      !pass_tight && pass_tau_tight);

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 3l1tau" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_mll_z = " << pass_mll_z << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_charge = " << pass_charge << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_tau_tight = " << pass_tau_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_3l1tau_SR = " << pass_3l1tau_SR << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  && pass_tau_tight = " << pass_tau_tight << std::endl;
			    std::cout << "  = pass_3l1tau_Fake = " << pass_3l1tau_Fake << std::endl;
			 }		  
		    }
	       }
	  }

	bool pass_4l_SR = 0;
	bool pass_4l_Fake = 0;
	  {
	     bool pass_nlep = (nLepFakeable > 3);
	     if( pass_nlep )
	       {
		  bool pass_trig = (nElectronFakeable >= 1 && trig_e) ||
		    (nMuonFakeable >= 1 && trig_m) ||
		    (nElectronFakeable >= 2 && trig_ee) ||
		    (nMuonFakeable >= 2 && trig_mm) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 1 && trig_em) ||
		    (nElectronFakeable >= 3 && trig_eee) ||
		    (nElectronFakeable >= 2 && nMuonFakeable >= 1 && trig_eem) ||
		    (nElectronFakeable >= 1 && nMuonFakeable >= 2 && trig_emm) ||
		    (nMuonFakeable >= 3 && trig_mmm);
		  bool pass_fakeable_pt = (elmuFakeable->at(0).conept > 25. && elmuFakeable->at(1).conept > 15. && elmuFakeable->at(2).conept > 15. && elmuFakeable->at(3).conept > 10.);
		  bool pass_metLD = (((metLD > 0.2 && nSFOS == 0) || (metLD > 0.3 && nSFOS > 0)) || nJetLoose >= 4);
		  int lep_charge_sum = elmuFakeable->at(0).charge+elmuFakeable->at(1).charge+elmuFakeable->at(2).charge+elmuFakeable->at(3).charge;
		  bool pass_charge = (lep_charge_sum == 0);
		  bool pass_njet = (nJetLoose >= 2 && (nJetLooseBL >= 2 || nJetLooseBM >= 1));
		  bool pass_mllll = (mllll > 140 || mllll < 0);

		  bool pass_tight = (elmuFakeable->at(0).isTightTTH && elmuFakeable->at(1).isTightTTH && elmuFakeable->at(2).isTightTTH && elmuFakeable->at(3).isTightTTH);
		  bool pass_truth = (elmuFakeable->at(0).hasMCMatch && elmuFakeable->at(1).hasMCMatch && elmuFakeable->at(2).hasMCMatch && elmuFakeable->at(3).hasMCMatch);

		  pass_4l_SR = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet && pass_mllll &&
				pass_tight && pass_truth);
		  
		  pass_4l_Fake = (pass_trig && pass_fakeable_pt && pass_mll && pass_mll_z && pass_metLD && pass_charge && pass_njet && pass_mllll &&
				  !pass_tight);

		  for(int d=0;d<evdebug->size();d++)
		    {		       
		       int evId = ntP->ev_id;
		       if( evId == evdebug->at(d) )
			 {
			    std::cout << "------------------------------" << std::endl;
			    std::cout << "Event #" << evId << std::endl;
			    std::cout << "Category 4l" << std::endl;
			    std::cout << "  pass_trig = " << pass_trig << std::endl;
			    std::cout << "  pass_fakeable_pt = " << pass_fakeable_pt << std::endl;
			    std::cout << "  pass_mll = " << pass_mll << std::endl;
			    std::cout << "  pass_mll_z = " << pass_mll_z << std::endl;
			    std::cout << "  pass_metLD = " << pass_metLD << std::endl;
			    std::cout << "  pass_charge = " << pass_charge << std::endl;
			    std::cout << "  pass_njet = " << pass_njet << std::endl;
			    std::cout << "  pass_mllll = " << pass_mllll << std::endl;
			    std::cout << "  && pass_tight = " << pass_tight << std::endl;
			    std::cout << "  && pass_truth = " << pass_truth << std::endl;
			    std::cout << "  = pass_4l_SR = " << pass_4l_SR << std::endl;
			    std::cout << "  && !pass_tight = " << !pass_tight << std::endl;
			    std::cout << "  = pass_4l_Fake = " << pass_4l_Fake << std::endl;
			 }		  
		    }
	       }	     
	  }

	ev->is_1l2tau_SR = pass_1l2tau_SR;
	ev->is_1l2tau_Fake = pass_1l2tau_Fake;
	ev->is_2lSS_SR = pass_2lSS_SR;
	ev->is_2lSS_SR = pass_2lSS_SR;
	ev->is_2lSS_Fake = pass_2lSS_Fake;
	ev->is_2lSS_Flip = pass_2lSS_Flip;
	ev->is_2lSS1tau_SR = pass_2lSS1tau_SR;
	ev->is_2lSS1tau_Fake = pass_2lSS1tau_Fake;
	ev->is_2lSS1tau_Flip = pass_2lSS1tau_Flip;
	ev->is_2l2tau_SR = pass_2l2tau_SR;
	ev->is_2l2tau_Fake = pass_2l2tau_Fake;
	ev->is_3l_SR = pass_3l_SR;
	ev->is_3l_Fake = pass_3l_Fake;
	ev->is_3l1tau_SR = pass_3l1tau_SR;
	ev->is_3l1tau_Fake = pass_3l1tau_Fake;
	ev->is_4l_SR = pass_4l_SR;
	ev->is_4l_Fake = pass_4l_Fake;
	ev->is_ttWctrl_SR = pass_ttWctrl_SR;
	ev->is_ttWctrl_Fake = pass_ttWctrl_Fake;
	ev->is_ttWctrl_Flip = pass_ttWctrl_Flip;
	ev->is_ttZctrl_SR = pass_ttZctrl_SR;
	ev->is_ttZctrl_Fake = pass_ttZctrl_Fake;

	if( sync == 2 )
	  {	     
	     if( pass_1l2tau_SR ) m_tree_1l2tau_SR->Fill();
	     if( pass_1l2tau_Fake ) m_tree_1l2tau_Fake->Fill();
	     if( pass_2lSS_SR ) m_tree_2lSS_SR->Fill();
	     if( pass_2lSS_Fake ) m_tree_2lSS_Fake->Fill();
	     if( pass_2lSS_Flip ) m_tree_2lSS_Flip->Fill();
	     if( pass_2lSS1tau_SR ) m_tree_2lSS1tau_SR->Fill();
	     if( pass_2lSS1tau_Fake ) m_tree_2lSS1tau_Fake->Fill();
	     if( pass_2lSS1tau_Flip ) m_tree_2lSS1tau_Flip->Fill();
	     if( pass_2l2tau_SR ) m_tree_2l2tau_SR->Fill();
	     if( pass_2l2tau_Fake ) m_tree_2l2tau_Fake->Fill();
	     if( pass_3l_SR ) m_tree_3l_SR->Fill();
	     if( pass_3l_Fake ) m_tree_3l_Fake->Fill();
	     if( pass_3l1tau_SR ) m_tree_3l1tau_SR->Fill();
	     if( pass_3l1tau_Fake ) m_tree_3l1tau_Fake->Fill();
	     if( pass_4l_SR ) m_tree_4l_SR->Fill();
	     if( pass_4l_Fake ) m_tree_4l_Fake->Fill();
	     if( pass_ttWctrl_SR ) m_tree_ttWctrl_SR->Fill();
	     if( pass_ttWctrl_Fake ) m_tree_ttWctrl_Fake->Fill();
	     if( pass_ttWctrl_Flip ) m_tree_ttWctrl_Flip->Fill();
	     if( pass_ttZctrl_SR ) m_tree_ttZctrl_SR->Fill();
	     if( pass_ttZctrl_Fake ) m_tree_ttZctrl_Fake->Fill();
	  }	
	
	delete elmuLoose;
	delete elmuTight;
	delete elmuFakeable;
	
	int pass_cat[21];
	std::string name_cat[21];
	
	pass_cat[0] = pass_1l2tau_SR;        name_cat[0] = "1l2tau_SR";
	pass_cat[1] = pass_1l2tau_Fake;      name_cat[1] = "1l2tau_Fake";
	pass_cat[2] = pass_2lSS_SR;          name_cat[2] = "2lSS_SR";
	pass_cat[3] = pass_2lSS_Fake;        name_cat[3] = "2lSS_Fake";
	pass_cat[4] = pass_2lSS_Flip;        name_cat[4] = "2lSS_Flip";
	pass_cat[5] = pass_2lSS1tau_SR;      name_cat[5] = "2lSS1tau_SR";
	pass_cat[6] = pass_2lSS1tau_Fake;    name_cat[6] = "2lSS1tau_Fake";
	pass_cat[7] = pass_2lSS1tau_Flip;    name_cat[7] = "2lSS1tau_Flip";
	pass_cat[8] = pass_2l2tau_SR;        name_cat[8] = "2l2tau_SR";
	pass_cat[9] = pass_2l2tau_Fake;      name_cat[9] = "2l2tau_Fake";
	pass_cat[10] = pass_3l_SR;           name_cat[10] = "3l_SR";
	pass_cat[11] = pass_3l_Fake;         name_cat[11] = "3l_Fake";
	pass_cat[12] = pass_3l1tau_SR;       name_cat[12] = "3l1tau_SR";	 
	pass_cat[13] = pass_3l1tau_Fake;     name_cat[13] = "3l1tau_Fake";
	pass_cat[14] = pass_4l_SR;           name_cat[14] = "4l_SR";
	pass_cat[15] = pass_4l_Fake;         name_cat[15] = "4l_Fake";
	pass_cat[16] = pass_ttWctrl_SR;      name_cat[16] = "ttWctrl_SR";
	pass_cat[17] = pass_ttWctrl_Fake;    name_cat[17] = "ttWctrl_Fake";
	pass_cat[18] = pass_ttWctrl_Flip;    name_cat[18] = "ttWctrl_Flip";
	pass_cat[19] = pass_ttZctrl_SR;      name_cat[19] = "ttZctrl_SR";
	pass_cat[20] = pass_ttZctrl_Fake;    name_cat[20] = "ttZctrl_Fake";

	int pass_event = 0;
	for(int ic=0;ic<21;ic++)
	  {	     
	     pass_event += pass_cat[ic];
	     if( sync != 0 )
	       {	    
		  m_hist_overlap->GetXaxis()->SetBinLabel(ic+1,name_cat[ic].c_str());
		  m_hist_overlap->GetYaxis()->SetBinLabel(ic+1,name_cat[ic].c_str());
	       }		  
	  }	
	
	if( pass_event == 0 ) pass = 0;
	else
	  {
	     if( sync != 0 )
	       {	    
		  for(int ic=0;ic<21;ic++)
		    {
		       if( pass_cat[ic] == 1 )
			 {
			    for(int icc=0;icc<21;icc++)
			      {		       
				 if( pass_cat[icc] == 1 )
				   {			    
				      m_hist_overlap->Fill(ic,icc);
				   }			    
			      }		       
			 }		  
		    }
	       }		  
	  }	
     }   
   
   return pass;
}
