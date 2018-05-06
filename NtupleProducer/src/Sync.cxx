#include "include/Sync.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>

Sync::Sync(std::string fname_out)
{
   _fname_out = fname_out;
}

Sync::~Sync()
{   
   m_file->Write();
   m_file->Close();
}

void Sync::Init(int sync)
{
   m_file = new TFile(_fname_out.c_str(),"RECREATE");
   if( sync == 1 ) m_tree = new TTree("syncTree","Sync Ntuple");
   else
     {
	m_tree_1l2tau_SR = new TTree("syncTree_1l2tau_SR","Sync Ntuple");
	m_tree_1l2tau_Fake = new TTree("syncTree_1l2tau_Fake","Sync Ntuple");
     }   
}

void Sync::setBranchAddress(int sync)
{
   if( sync == 1 ) createBranch(m_tree);
   else
     {
	createBranch(m_tree_1l2tau_SR);
	createBranch(m_tree_1l2tau_Fake);
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
   isGenMatched = -9999;
   
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
   nEvent = nt->NtEvent->at(0).id();
   ls = nt->NtEvent->at(0).lumi();
   run = nt->NtEvent->at(0).run();
   n_presel_mu = n_presel_mu;
   n_fakeablesel_mu = -9999;
   n_mvasel_mu = -9999;
   n_presel_ele = n_presel_el;
   n_fakeablesel_ele = -9999;
   n_mvasel_ele = -9999;
   n_presel_tau = n_presel_tau;
   n_presel_jet = n_presel_jet;
   
   int nMuon = nt->NtMuonLoose->size();
   
   if( nMuon > 0 )
     {	
	mu1_pt = nt->NtMuonLoose->at(0).pt();
	mu1_conept = nt->NtMuonLoose->at(0).conept();
	mu1_eta = nt->NtMuonLoose->at(0).eta();
	mu1_phi = nt->NtMuonLoose->at(0).phi();
	mu1_E = nt->NtMuonLoose->at(0).E();
	mu1_charge = nt->NtMuonLoose->at(0).charge();
	mu1_miniRelIso = nt->NtMuonLoose->at(0).iso();
	mu1_miniIsoCharged = nt->NtMuonLoose->at(0).lepMVA_miniRelIsoCharged()*mu1_pt;
	mu1_miniIsoNeutral = nt->NtMuonLoose->at(0).lepMVA_miniRelIsoNeutral()*mu1_pt;
	mu1_PFRelIso04 = nt->NtMuonLoose->at(0).isoR04();
	mu1_jetNDauChargedMVASel = nt->NtMuonLoose->at(0).lepMVA_jetNDauChargedMVASel();
	mu1_jetPtRel = nt->NtMuonLoose->at(0).lepMVA_jetPtRelv2();
	mu1_jetPtRatio = nt->NtMuonLoose->at(0).lepMVA_jetPtRatio();
	mu1_jetCSV = nt->NtMuonLoose->at(0).lepMVA_jetBTagCSV();
	mu1_sip3D = fabs(nt->NtMuonLoose->at(0).sip3d());
	mu1_dxy = nt->NtMuonLoose->at(0).dxy();
	mu1_dxyAbs = fabs(nt->NtMuonLoose->at(0).dxy());
	mu1_dz = nt->NtMuonLoose->at(0).dz();
	mu1_segmentCompatibility = nt->NtMuonLoose->at(0).lepMVA_mvaId();
	mu1_leptonMVA = nt->NtMuonLoose->at(0).lepMVA();
	mu1_mediumID = nt->NtMuonLoose->at(0).isMedium();
	mu1_dpt_div_pt = nt->NtMuonLoose->at(0).bestTrackptError()/nt->NtMuonLoose->at(0).bestTrackpt();
	mu1_isfakeablesel = nt->NtMuonLoose->at(0).isFakeableTTH();
	mu1_ismvasel = nt->NtMuonLoose->at(0).isTightTTH();
     }   

   if( nMuon > 1 )
     {	
	mu2_pt = nt->NtMuonLoose->at(1).pt();
	mu2_conept = nt->NtMuonLoose->at(1).conept();
	mu2_eta = nt->NtMuonLoose->at(1).eta();
	mu2_phi = nt->NtMuonLoose->at(1).phi();
	mu2_E = nt->NtMuonLoose->at(1).E();
	mu2_charge = nt->NtMuonLoose->at(1).charge();
	mu2_miniRelIso = nt->NtMuonLoose->at(1).iso();
	mu2_miniIsoCharged = nt->NtMuonLoose->at(1).lepMVA_miniRelIsoCharged()*mu2_pt;
	mu2_miniIsoNeutral = nt->NtMuonLoose->at(1).lepMVA_miniRelIsoNeutral()*mu2_pt;
	mu2_PFRelIso04 = nt->NtMuonLoose->at(1).isoR04();
	mu2_jetNDauChargedMVASel = nt->NtMuonLoose->at(1).lepMVA_jetNDauChargedMVASel();
	mu2_jetPtRel = nt->NtMuonLoose->at(1).lepMVA_jetPtRelv2();
	mu2_jetPtRatio = nt->NtMuonLoose->at(1).lepMVA_jetPtRatio();
	mu2_jetCSV = nt->NtMuonLoose->at(1).lepMVA_jetBTagCSV();
	mu2_sip3D = fabs(nt->NtMuonLoose->at(1).sip3d());
	mu2_dxy = nt->NtMuonLoose->at(1).dxy();
	mu2_dxyAbs = fabs(nt->NtMuonLoose->at(1).dxy());
	mu2_dz = nt->NtMuonLoose->at(1).dz();
	mu2_segmentCompatibility = nt->NtMuonLoose->at(1).lepMVA_mvaId();
	mu2_leptonMVA = nt->NtMuonLoose->at(1).lepMVA();
	mu2_mediumID = nt->NtMuonLoose->at(1).isMedium();
	mu2_dpt_div_pt = nt->NtMuonLoose->at(1).bestTrackptError()/nt->NtMuonLoose->at(1).bestTrackpt();
	mu2_isfakeablesel = nt->NtMuonLoose->at(1).isFakeableTTH();
	mu2_ismvasel = nt->NtMuonLoose->at(1).isTightTTH();
     }

   int nElectron = nt->NtElectronLoose->size();
   
   if( nElectron > 0 )
     {	   
	ele1_pt = nt->NtElectronLoose->at(0).pt();
	ele1_conept = nt->NtElectronLoose->at(0).conept();
	ele1_eta = nt->NtElectronLoose->at(0).eta();
	ele1_phi = nt->NtElectronLoose->at(0).phi();
	ele1_E = nt->NtElectronLoose->at(0).E();
	ele1_charge = nt->NtElectronLoose->at(0).charge();
	ele1_miniRelIso = nt->NtElectronLoose->at(0).miniIso();
	ele1_miniIsoCharged = nt->NtElectronLoose->at(0).lepMVA_miniRelIsoCharged()*ele1_pt;
	ele1_miniIsoNeutral = nt->NtElectronLoose->at(0).lepMVA_miniRelIsoNeutral()*ele1_pt;
	ele1_PFRelIso04 = nt->NtElectronLoose->at(0).isoR04();
	ele1_jetNDauChargedMVASel = nt->NtElectronLoose->at(0).lepMVA_jetNDauChargedMVASel();
	ele1_jetPtRel = nt->NtElectronLoose->at(0).lepMVA_jetPtRelv2();;
	ele1_jetPtRatio = nt->NtElectronLoose->at(0).lepMVA_jetPtRatio();
	ele1_jetCSV = nt->NtElectronLoose->at(0).lepMVA_jetBTagCSV();
	ele1_sip3D = fabs(nt->NtElectronLoose->at(0).sip3d());
	ele1_dxy = nt->NtElectronLoose->at(0).dxy();
	ele1_dxyAbs = fabs(nt->NtElectronLoose->at(0).dxy());
	ele1_dz = nt->NtElectronLoose->at(0).dz();
	ele1_ntMVAeleID = nt->NtElectronLoose->at(0).lepMVA_mvaId();
	ele1_leptonMVA = nt->NtElectronLoose->at(0).lepMVA();
	ele1_isChargeConsistent = nt->NtElectronLoose->at(0).isGsfScPixChargeConsistent();
	ele1_passesConversionVeto = nt->NtElectronLoose->at(0).passCV();
	ele1_nMissingHits = nt->NtElectronLoose->at(0).nlosthits();
	ele1_sigmaEtaEta = nt->NtElectronLoose->at(0).sigmaIetaIeta();
	ele1_HoE = nt->NtElectronLoose->at(0).hadronicOverEm();
	ele1_deltaEta = nt->NtElectronLoose->at(0).deltaEtaSuperClusterTrackAtVtx();
	ele1_deltaPhi = nt->NtElectronLoose->at(0).deltaPhiSuperClusterTrackAtVtx();
	ele1_OoEminusOoP = nt->NtElectronLoose->at(0).ooEmooP();
	ele1_isfakeablesel = nt->NtElectronLoose->at(0).isFakeableTTH();
	ele1_ismvasel = nt->NtElectronLoose->at(0).isTightTTH();
     }   

   if( nElectron > 1 )
     {	   
	ele2_pt = nt->NtElectronLoose->at(1).pt();
	ele2_conept = nt->NtElectronLoose->at(1).conept();
	ele2_eta = nt->NtElectronLoose->at(1).eta();
	ele2_phi = nt->NtElectronLoose->at(1).phi();
	ele2_E = nt->NtElectronLoose->at(1).E();
	ele2_charge = nt->NtElectronLoose->at(1).charge();
	ele2_miniRelIso = nt->NtElectronLoose->at(1).miniIso();
	ele2_miniIsoCharged = nt->NtElectronLoose->at(1).lepMVA_miniRelIsoCharged()*ele2_pt;
	ele2_miniIsoNeutral = nt->NtElectronLoose->at(1).lepMVA_miniRelIsoNeutral()*ele2_pt;
	ele2_PFRelIso04 = nt->NtElectronLoose->at(1).isoR04();
	ele2_jetNDauChargedMVASel = nt->NtElectronLoose->at(1).lepMVA_jetNDauChargedMVASel();
	ele2_jetPtRel = nt->NtElectronLoose->at(1).lepMVA_jetPtRelv2();;
	ele2_jetPtRatio = nt->NtElectronLoose->at(1).lepMVA_jetPtRatio();
	ele2_jetCSV = nt->NtElectronLoose->at(1).lepMVA_jetBTagCSV();
	ele2_sip3D = fabs(nt->NtElectronLoose->at(1).sip3d());
	ele2_dxy = nt->NtElectronLoose->at(1).dxy();
	ele2_dxyAbs = fabs(nt->NtElectronLoose->at(1).dxy());
	ele2_dz = nt->NtElectronLoose->at(1).dz();
	ele2_ntMVAeleID = nt->NtElectronLoose->at(1).lepMVA_mvaId();
	ele2_leptonMVA = nt->NtElectronLoose->at(1).lepMVA();
	ele2_isChargeConsistent = nt->NtElectronLoose->at(1).isGsfScPixChargeConsistent();
	ele2_passesConversionVeto = nt->NtElectronLoose->at(1).passCV();
	ele2_nMissingHits = nt->NtElectronLoose->at(1).nlosthits();
	ele2_sigmaEtaEta = nt->NtElectronLoose->at(1).sigmaIetaIeta();
	ele2_HoE = nt->NtElectronLoose->at(1).hadronicOverEm();
	ele2_deltaEta = nt->NtElectronLoose->at(1).deltaEtaSuperClusterTrackAtVtx();
	ele2_deltaPhi = nt->NtElectronLoose->at(1).deltaPhiSuperClusterTrackAtVtx();
	ele2_OoEminusOoP = nt->NtElectronLoose->at(1).ooEmooP();
	ele2_isfakeablesel = nt->NtElectronLoose->at(1).isFakeableTTH();
	ele2_ismvasel = nt->NtElectronLoose->at(1).isTightTTH();
     }   

   int nTau = nt->NtTauFakeable->size();

   if( nTau > 0 )
     {	      
	tau1_pt = nt->NtTauFakeable->at(0).pt();
	tau1_eta = nt->NtTauFakeable->at(0).eta();
	tau1_phi = nt->NtTauFakeable->at(0).phi();
	tau1_E = nt->NtTauFakeable->at(0).E();
	tau1_charge = nt->NtTauFakeable->at(0).charge();
	tau1_dxy = nt->NtTauFakeable->at(0).dxy();
	tau1_dz = nt->NtTauFakeable->at(0).dz();
	tau1_decayModeFindingOldDMs = nt->NtTauFakeable->at(0).decayModeFinding();
	tau1_decayModeFindingNewDMs = nt->NtTauFakeable->at(0).decayModeFindingNewDMs();
	tau1_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(0).byLooseCombinedIsolationDeltaBetaCorr3Hits();
	tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(0).byMediumCombinedIsolationDeltaBetaCorr3Hits();
	tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(0).byTightCombinedIsolationDeltaBetaCorr3Hits();
	tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(0).byVLooseIsolationMVArun2v1DBdR03oldDMwLT();
	tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(0).byLooseIsolationMVArun2v1DBdR03oldDMwLT();
	tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(0).byMediumIsolationMVArun2v1DBdR03oldDMwLT();
	tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(0).byTightIsolationMVArun2v1DBdR03oldDMwLT();
	tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(0).byVTightIsolationMVArun2v1DBdR03oldDMwLT();
	tau1_rawMVArun2v1DBdR03oldDMwLT = -9999;
	tau1_againstMuonLoose3 = nt->NtTauFakeable->at(0).againstMuonLoose3();
	tau1_againstMuonTight3 = nt->NtTauFakeable->at(0).againstMuonTight3();
	tau1_againstElectronVLooseMVA6 = nt->NtTauFakeable->at(0).againstElectronVLooseMVA6();
	tau1_againstElectronLooseMVA6 = nt->NtTauFakeable->at(0).againstElectronLooseMVA6();
	tau1_againstElectronMediumMVA6 = nt->NtTauFakeable->at(0).againstElectronMediumMVA6();
	tau1_againstElectronTightMVA6 = nt->NtTauFakeable->at(0).againstElectronTightMVA6();
     }   

   if( nTau > 1 )
     {	      
	tau2_pt = nt->NtTauFakeable->at(1).pt();
	tau2_eta = nt->NtTauFakeable->at(1).eta();
	tau2_phi = nt->NtTauFakeable->at(1).phi();
	tau2_E = nt->NtTauFakeable->at(1).E();
	tau2_charge = nt->NtTauFakeable->at(1).charge();
	tau2_dxy = nt->NtTauFakeable->at(1).dxy();
	tau2_dz = nt->NtTauFakeable->at(1).dz();
	tau2_decayModeFindingOldDMs = nt->NtTauFakeable->at(1).decayModeFinding();
	tau2_decayModeFindingNewDMs = nt->NtTauFakeable->at(1).decayModeFindingNewDMs();
	tau2_byCombinedIsolationDeltaBetaCorr3Hits = -9999;
	tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(1).byLooseCombinedIsolationDeltaBetaCorr3Hits();
	tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(1).byMediumCombinedIsolationDeltaBetaCorr3Hits();
	tau2_byTightCombinedIsolationDeltaBetaCorr3Hits = nt->NtTauFakeable->at(1).byTightCombinedIsolationDeltaBetaCorr3Hits();
	tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
	tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(1).byVLooseIsolationMVArun2v1DBdR03oldDMwLT();
	tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(1).byLooseIsolationMVArun2v1DBdR03oldDMwLT();
	tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(1).byMediumIsolationMVArun2v1DBdR03oldDMwLT();
	tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(1).byTightIsolationMVArun2v1DBdR03oldDMwLT();
	tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT = nt->NtTauFakeable->at(1).byVTightIsolationMVArun2v1DBdR03oldDMwLT();
	tau2_rawMVArun2v1DBdR03oldDMwLT = -9999;
	tau2_againstMuonLoose3 = nt->NtTauFakeable->at(1).againstMuonLoose3();
	tau2_againstMuonTight3 = nt->NtTauFakeable->at(1).againstMuonTight3();
	tau2_againstElectronVLooseMVA6 = nt->NtTauFakeable->at(1).againstElectronVLooseMVA6();
	tau2_againstElectronLooseMVA6 = nt->NtTauFakeable->at(1).againstElectronLooseMVA6();
	tau2_againstElectronMediumMVA6 = nt->NtTauFakeable->at(1).againstElectronMediumMVA6();
	tau2_againstElectronTightMVA6 = nt->NtTauFakeable->at(1).againstElectronTightMVA6();
     }   
}

void Sync::fill(Ntuple *nt,int sync)
{
   if( sync == 1 ) m_tree->Fill();
   else
     {
	int nElectronFakeable = nt->NtElectronFakeable->size();
	int nMuonFakeable = nt->NtMuonFakeable->size();
	int nTauFakeable = nt->NtTauFakeable->size();
	int nLepFakeable = nElectronFakeable+nMuonFakeable;
	int nLepTight = nt->NtElectronTight->size()+nt->NtMuonTight->size();
	bool pass_ele_pt = (nElectronFakeable > 0) ? (nt->NtElectronFakeable->at(0).pt() > 30.) : 0;
	bool pass_muo_pt = (nMuonFakeable > 0) ? (nt->NtMuonFakeable->at(0).pt() > 25.) : 0;
	bool pass_lep_pt = (pass_ele_pt || pass_muo_pt);
	bool pass_ele_eta = (nElectronFakeable > 0) ? (fabs(nt->NtElectronFakeable->at(0).eta()) < 2.1) : 0;
	bool pass_muo_eta = (nMuonFakeable > 0) ? (fabs(nt->NtMuonFakeable->at(0).eta()) < 2.1) : 0;
	bool pass_lep_eta = (pass_ele_eta || pass_muo_eta);
	bool pass_tau_pt = (nTauFakeable >= 2) ? (nt->NtTauFakeable->at(0).pt() > 30. && nt->NtTauFakeable->at(1).pt() > 20.) : 0;
	int nJetLoose = nt->NtJetLoose->size();
	int nJetLooseBL = 0;
	int nJetLooseBM = 0;
	for(int ij=0;ij<nJetLoose;ij++)
	  {
	     if( nt->NtJetLoose->at(ij).isMediumBTag() ) nJetLooseBM++;
	     if( nt->NtJetLoose->at(ij).isLooseBTag() ) nJetLooseBL++;
	  }	
	
	bool pass_mee = 1;
	for(int ie=0;ie<nt->NtElectronLoose->size();ie++ )
	  {
	     TLorentzVector *l1 = new TLorentzVector();
	     l1->SetPtEtaPhiE(nt->NtElectronLoose->at(ie).pt(),
			      nt->NtElectronLoose->at(ie).eta(),
			      nt->NtElectronLoose->at(ie).phi(),
			      nt->NtElectronLoose->at(ie).E());
	     
	     for(int iee=ie+1;iee<nt->NtElectronLoose->size();iee++ )
	       {
		  TLorentzVector *l2 = new TLorentzVector();
		  l2->SetPtEtaPhiE(nt->NtElectronLoose->at(iee).pt(),
				   nt->NtElectronLoose->at(iee).eta(),
				   nt->NtElectronLoose->at(iee).phi(),
				   nt->NtElectronLoose->at(iee).E());
		  
		  float mll = (*l1+*l2).M();

		  if( mll < 12 ) pass_mee = 0;
		  
		  delete l2;
	       }	    
	     
	     delete l1;
	  }		

	bool pass_mmm = 1;
	for(int im=0;im<nt->NtMuonLoose->size();im++ )
	  {
	     TLorentzVector *l1 = new TLorentzVector();
	     l1->SetPtEtaPhiE(nt->NtMuonLoose->at(im).pt(),
			      nt->NtMuonLoose->at(im).eta(),
			      nt->NtMuonLoose->at(im).phi(),
			      nt->NtMuonLoose->at(im).E());
	     
	     for(int imm=im+1;imm<nt->NtMuonLoose->size();imm++ )
	       {
		  TLorentzVector *l2 = new TLorentzVector();
		  l2->SetPtEtaPhiE(nt->NtMuonLoose->at(imm).pt(),
				   nt->NtMuonLoose->at(imm).eta(),
				   nt->NtMuonLoose->at(imm).phi(),
				   nt->NtMuonLoose->at(imm).E());
		  
		  float mll = (*l1+*l2).M();
		  if( mll < 12 ) pass_mmm = 0;
		  
		  delete l2;
	       }	    
	     
	     delete l1;
	  }		

	bool pass_mem = 1;
	for(int ie=0;ie<nt->NtElectronLoose->size();ie++ )
	  {
	     TLorentzVector *l1 = new TLorentzVector();
	     l1->SetPtEtaPhiE(nt->NtElectronLoose->at(ie).pt(),
			      nt->NtElectronLoose->at(ie).eta(),
			      nt->NtElectronLoose->at(ie).phi(),
			      nt->NtElectronLoose->at(ie).E());
	     
	     for(int im=0;im<nt->NtMuonLoose->size();im++ )
	       {
		  TLorentzVector *l2 = new TLorentzVector();
		  l2->SetPtEtaPhiE(nt->NtMuonLoose->at(im).pt(),
				   nt->NtMuonLoose->at(im).eta(),
				   nt->NtMuonLoose->at(im).phi(),
				   nt->NtMuonLoose->at(im).E());
		  
		  float mll = (*l1+*l2).M();
		  if( mll < 12 ) pass_mem = 0;
		  
		  delete l2;
	       }	    
	     
	     delete l1;
	  }		
	
	bool pass_mll = (pass_mee && pass_mmm && pass_mem);
	
	if( nLepFakeable > 0 && pass_lep_pt && pass_lep_eta && nLepTight <= 1 && nTauFakeable >= 2 && pass_tau_pt &&
	    nJetLoose >= 3 && (nJetLooseBL >= 2 || nJetLooseBM >= 1) && pass_mll )
	  {
	     if( ((nElectronFakeable > 0 && nt->NtElectronLoose->at(0).isTightTTH()) ||
		  (nMuonFakeable > 0 && nt->NtMuonLoose->at(0).isTightTTH())) &&
		 (nt->NtTauFakeable->at(0).byMediumIsolationMVArun2v1DBdR03oldDMwLT() && 
		     nt->NtTauFakeable->at(1).byMediumIsolationMVArun2v1DBdR03oldDMwLT()) &&
		 (nt->NtTauFakeable->at(0).charge()*nt->NtTauFakeable->at(1).charge() < 0) )
	       {
		  m_tree_1l2tau_SR->Fill();
	       }
	     if( ((nElectronFakeable > 0 && !nt->NtElectronLoose->at(0).isTightTTH()) &&
		  (nMuonFakeable > 0 && !nt->NtMuonLoose->at(0).isTightTTH())) ||
		 (nt->NtTauFakeable->at(0).byMediumIsolationMVArun2v1DBdR03oldDMwLT() ||
		     nt->NtTauFakeable->at(1).byMediumIsolationMVArun2v1DBdR03oldDMwLT()) &&
		 (nt->NtTauFakeable->at(0).charge()*nt->NtTauFakeable->at(1).charge() < 0) )
	       {
		  m_tree_1l2tau_Fake->Fill();
	       }	     
	  }	
     }   
}
