#include "include/Sync.h"

Sync::Sync(std::string fname_out)
{
   _fname_out = fname_out;
}

Sync::~Sync()
{   
   m_file->Write();
   m_file->Close();
}

void Sync::Init()
{
   m_file = new TFile(_fname_out.c_str(),"RECREATE");
   m_tree = new TTree("Nt","Sync Ntuple");
}

void Sync::setBranchAddress()
{
   m_tree->Branch("nEvent",&nEvent,"nEvent/I");
   m_tree->Branch("ls",&ls,"ls/I");
   m_tree->Branch("run",&run,"run/I");
   m_tree->Branch("n_presel_mu",&n_presel_mu,"n_presel_mu/I");
   m_tree->Branch("n_fakeablesel_mu",&n_fakeablesel_mu,"n_fakeablesel_mu/I");
   m_tree->Branch("n_mvasel_mu",&n_mvasel_mu,"n_mvasel_mu/I");
   m_tree->Branch("n_presel_ele",&n_presel_ele,"n_presel_ele/I");
   m_tree->Branch("n_fakeablesel_ele",&n_fakeablesel_ele,"n_fakeablesel_ele/I");
   m_tree->Branch("n_mvasel_ele",&n_mvasel_ele,"n_mvasel_ele/I");
   m_tree->Branch("n_presel_tau",&n_presel_tau,"n_presel_tau/I");
   m_tree->Branch("n_presel_jet",&n_presel_jet,"n_presel_jet/I");
   
   m_tree->Branch("mu1_pt",&mu1_pt,"mu1_pt/F");
   m_tree->Branch("mu1_conept",&mu1_conept,"mu1_conept/F");
   m_tree->Branch("mu1_eta",&mu1_eta,"mu1_eta/F");
   m_tree->Branch("mu1_phi",&mu1_phi,"mu1_phi/F");
   m_tree->Branch("mu1_E",&mu1_E,"mu1_E/F");
   m_tree->Branch("mu1_charge",&mu1_charge,"mu1_charge/I");
   m_tree->Branch("mu1_miniRelIso",&mu1_miniRelIso,"mu1_miniRelIso/F");
   m_tree->Branch("mu1_miniIsoCharged",&mu1_miniIsoCharged,"mu1_miniIsoCharged/F");
   m_tree->Branch("mu1_miniIsoNeutral",&mu1_miniIsoNeutral,"mu1_miniIsoNeutral/F");
   m_tree->Branch("mu1_PFRelIso04",&mu1_PFRelIso04,"mu1_PFRelIso04/F");
   m_tree->Branch("mu1_jetNDauChargedMVASel",&mu1_jetNDauChargedMVASel,"mu1_jetNDauChargedMVASel/F");
   m_tree->Branch("mu1_jetPtRel",&mu1_jetPtRel,"mu1_jetPtRel/F");
   m_tree->Branch("mu1_jetPtRatio",&mu1_jetPtRatio,"mu1_jetPtRatio/F");
   m_tree->Branch("mu1_jetCSV",&mu1_jetCSV,"mu1_jetCSV/F");
   m_tree->Branch("mu1_sip3D",&mu1_sip3D,"mu1_sip3D/F");
   m_tree->Branch("mu1_dxy",&mu1_dxy,"mu1_dxy/F");
   m_tree->Branch("mu1_dxyAbs",&mu1_dxyAbs,"mu1_dxyAbs/F");
   m_tree->Branch("mu1_dz",&mu1_dz,"mu1_dz/F");
   m_tree->Branch("mu1_segmentCompatibility",&mu1_segmentCompatibility,"mu1_segmentCompatibility/F");
   m_tree->Branch("mu1_leptonMVA",&mu1_leptonMVA,"mu1_leptonMVA/F");
   m_tree->Branch("mu1_mediumID",&mu1_mediumID,"mu1_mediumID/F");
   m_tree->Branch("mu1_dpt_div_pt",&mu1_dpt_div_pt,"mu1_dpt_div_pt/F");
   m_tree->Branch("mu1_isfakeablesel",&mu1_isfakeablesel,"mu1_isfakeablesel/I");
   m_tree->Branch("mu1_ismvasel",&mu1_ismvasel,"mu1_ismvasel/I");

   m_tree->Branch("mu2_pt",&mu2_pt,"mu2_pt/F");
   m_tree->Branch("mu2_conept",&mu2_conept,"mu2_conept/F");
   m_tree->Branch("mu2_eta",&mu2_eta,"mu2_eta/F");
   m_tree->Branch("mu2_phi",&mu2_phi,"mu2_phi/F");
   m_tree->Branch("mu2_E",&mu2_E,"mu2_E/F");
   m_tree->Branch("mu2_charge",&mu2_charge,"mu2_charge/I");
   m_tree->Branch("mu2_miniRelIso",&mu2_miniRelIso,"mu2_miniRelIso/F");
   m_tree->Branch("mu2_miniIsoCharged",&mu2_miniIsoCharged,"mu2_miniIsoCharged/F");
   m_tree->Branch("mu2_miniIsoNeutral",&mu2_miniIsoNeutral,"mu2_miniIsoNeutral/F");
   m_tree->Branch("mu2_PFRelIso04",&mu2_PFRelIso04,"mu2_PFRelIso04/F");
   m_tree->Branch("mu2_jetNDauChargedMVASel",&mu2_jetNDauChargedMVASel,"mu2_jetNDauChargedMVASel/F");
   m_tree->Branch("mu2_jetPtRel",&mu2_jetPtRel,"mu2_jetPtRel/F");
   m_tree->Branch("mu2_jetPtRatio",&mu2_jetPtRatio,"mu2_jetPtRatio/F");
   m_tree->Branch("mu2_jetCSV",&mu2_jetCSV,"mu2_jetCSV/F");
   m_tree->Branch("mu2_sip3D",&mu2_sip3D,"mu2_sip3D/F");
   m_tree->Branch("mu2_dxy",&mu2_dxy,"mu2_dxy/F");
   m_tree->Branch("mu2_dxyAbs",&mu2_dxyAbs,"mu2_dxyAbs/F");
   m_tree->Branch("mu2_dz",&mu2_dz,"mu2_dz/F");
   m_tree->Branch("mu2_segmentCompatibility",&mu2_segmentCompatibility,"mu2_segmentCompatibility/F");
   m_tree->Branch("mu2_leptonMVA",&mu2_leptonMVA,"mu2_leptonMVA/F");
   m_tree->Branch("mu2_mediumID",&mu2_mediumID,"mu2_mediumID/F");
   m_tree->Branch("mu2_dpt_div_pt",&mu2_dpt_div_pt,"mu2_dpt_div_pt/F");
   m_tree->Branch("mu2_isfakeablesel",&mu2_isfakeablesel,"mu2_isfakeablesel/I");
   m_tree->Branch("mu2_ismvasel",&mu2_ismvasel,"mu2_ismvasel/I");
   
   m_tree->Branch("ele1_pt",&ele1_pt,"ele1_pt/F");
   m_tree->Branch("ele1_conept",&ele1_conept,"ele1_conept/F");
   m_tree->Branch("ele1_eta",&ele1_eta,"ele1_eta/F");
   m_tree->Branch("ele1_phi",&ele1_phi,"ele1_phi/F");
   m_tree->Branch("ele1_E",&ele1_E,"ele1_E/F");
   m_tree->Branch("ele1_charge",&ele1_charge,"ele1_charge/I");
   m_tree->Branch("ele1_miniRelIso",&ele1_miniRelIso,"ele1_miniRelIso/F");
   m_tree->Branch("ele1_miniIsoCharged",&ele1_miniIsoCharged,"ele1_miniIsoCharged/F");
   m_tree->Branch("ele1_miniIsoNeutral",&ele1_miniIsoNeutral,"ele1_miniIsoNeutral/F");
   m_tree->Branch("ele1_PFRelIso04",&ele1_PFRelIso04,"ele1_PFRelIso04/F");
   m_tree->Branch("ele1_jetNDauChargedMVASel",&ele1_jetNDauChargedMVASel,"ele1_jetNDauChargedMVASel/F");
   m_tree->Branch("ele1_jetPtRel",&ele1_jetPtRel,"ele1_jetPtRel/F");
   m_tree->Branch("ele1_jetPtRatio",&ele1_jetPtRatio,"ele1_jetPtRatio/F");
   m_tree->Branch("ele1_jetCSV",&ele1_jetCSV,"ele1_jetCSV/F");
   m_tree->Branch("ele1_sip3D",&ele1_sip3D,"ele1_sip3D/F");
   m_tree->Branch("ele1_dxy",&ele1_dxy,"ele1_dxy/F");
   m_tree->Branch("ele1_dxyAbs",&ele1_dxyAbs,"ele1_dxyAbs/F");
   m_tree->Branch("ele1_dz",&ele1_dz,"ele1_dz/F");
   m_tree->Branch("ele1_ntMVAeleID",&ele1_ntMVAeleID,"ele1_ntMVAeleID/F");
   m_tree->Branch("ele1_leptonMVA",&ele1_leptonMVA,"ele1_leptonMVA/F");
   m_tree->Branch("ele1_isChargeConsistent",&ele1_isChargeConsistent,"ele1_isChargeConsistent/I");
   m_tree->Branch("ele1_passesConversionVeto",&ele1_passesConversionVeto,"ele1_passesConversionVeto/F");
   m_tree->Branch("ele1_nMissingHits",&ele1_nMissingHits,"ele1_nMissingHits/F");
   m_tree->Branch("ele1_sigmaEtaEta",&ele1_sigmaEtaEta,"ele1_sigmaEtaEta/F");
   m_tree->Branch("ele1_HoE",&ele1_HoE,"ele1_HoE/F");
   m_tree->Branch("ele1_deltaEta",&ele1_deltaEta,"ele1_deltaEta/F");
   m_tree->Branch("ele1_deltaPhi",&ele1_deltaPhi,"ele1_deltaPhi/F");
   m_tree->Branch("ele1_OoEminusOoP",&ele1_OoEminusOoP,"ele1_OoEminusOoP/F");
   m_tree->Branch("ele1_isfakeablesel",&ele1_isfakeablesel,"ele1_isfakeablesel/I");
   m_tree->Branch("ele1_ismvasel",&ele1_ismvasel,"ele1_ismvasel/I");

   m_tree->Branch("ele2_pt",&ele2_pt,"ele2_pt/F");
   m_tree->Branch("ele2_conept",&ele2_conept,"ele2_conept/F");
   m_tree->Branch("ele2_eta",&ele2_eta,"ele2_eta/F");
   m_tree->Branch("ele2_phi",&ele2_phi,"ele2_phi/F");
   m_tree->Branch("ele2_E",&ele2_E,"ele2_E/F");
   m_tree->Branch("ele2_charge",&ele2_charge,"ele2_charge/I");
   m_tree->Branch("ele2_miniRelIso",&ele2_miniRelIso,"ele2_miniRelIso/F");
   m_tree->Branch("ele2_miniIsoCharged",&ele2_miniIsoCharged,"ele2_miniIsoCharged/F");
   m_tree->Branch("ele2_miniIsoNeutral",&ele2_miniIsoNeutral,"ele2_miniIsoNeutral/F");
   m_tree->Branch("ele2_PFRelIso04",&ele2_PFRelIso04,"ele2_PFRelIso04/F");
   m_tree->Branch("ele2_jetNDauChargedMVASel",&ele2_jetNDauChargedMVASel,"ele2_jetNDauChargedMVASel/F");
   m_tree->Branch("ele2_jetPtRel",&ele2_jetPtRel,"ele2_jetPtRel/F");
   m_tree->Branch("ele2_jetPtRatio",&ele2_jetPtRatio,"ele2_jetPtRatio/F");
   m_tree->Branch("ele2_jetCSV",&ele2_jetCSV,"ele2_jetCSV/F");
   m_tree->Branch("ele2_sip3D",&ele2_sip3D,"ele2_sip3D/F");
   m_tree->Branch("ele2_dxy",&ele2_dxy,"ele2_dxy/F");
   m_tree->Branch("ele2_dxyAbs",&ele2_dxyAbs,"ele2_dxyAbs/F");
   m_tree->Branch("ele2_dz",&ele2_dz,"ele2_dz/F");
   m_tree->Branch("ele2_ntMVAeleID",&ele2_ntMVAeleID,"ele2_ntMVAeleID/F");
   m_tree->Branch("ele2_leptonMVA",&ele2_leptonMVA,"ele2_leptonMVA/F");
   m_tree->Branch("ele2_isChargeConsistent",&ele2_isChargeConsistent,"ele2_isChargeConsistent/I");
   m_tree->Branch("ele2_passesConversionVeto",&ele2_passesConversionVeto,"ele2_passesConversionVeto/F");
   m_tree->Branch("ele2_nMissingHits",&ele2_nMissingHits,"ele2_nMissingHits/F");
   m_tree->Branch("ele2_sigmaEtaEta",&ele2_sigmaEtaEta,"ele2_sigmaEtaEta/F");
   m_tree->Branch("ele2_HoE",&ele2_HoE,"ele2_HoE/F");
   m_tree->Branch("ele2_deltaEta",&ele2_deltaEta,"ele2_deltaEta/F");
   m_tree->Branch("ele2_deltaPhi",&ele2_deltaPhi,"ele2_deltaPhi/F");
   m_tree->Branch("ele2_OoEminusOoP",&ele2_OoEminusOoP,"ele2_OoEminusOoP/F");
   m_tree->Branch("ele2_isfakeablesel",&ele2_isfakeablesel,"ele2_isfakeablesel/I");
   m_tree->Branch("ele2_ismvasel",&ele2_ismvasel,"ele2_ismvasel/I");

   m_tree->Branch("tau1_pt",&tau1_pt,"tau1_pt/F");
   m_tree->Branch("tau1_eta",&tau1_eta,"tau1_eta/F");
   m_tree->Branch("tau1_phi",&tau1_phi,"tau1_phi/F");
   m_tree->Branch("tau1_E",&tau1_E,"tau1_E/F");
   m_tree->Branch("tau1_charge",&tau1_charge,"tau1_charge/I");
   m_tree->Branch("tau1_dxy",&tau1_dxy,"tau1_dxy/F");
   m_tree->Branch("tau1_dz",&tau1_dz,"tau1_dz/F");
   m_tree->Branch("tau1_decayModeFindingOldDMs",&tau1_decayModeFindingOldDMs,"tau1_decayModeFindingOldDMs/F");
   m_tree->Branch("tau1_decayModeFindingNewDMs",&tau1_decayModeFindingNewDMs,"tau1_decayModeFindingNewDMs/F");
   m_tree->Branch("tau1_byCombinedIsolationDeltaBetaCorr3Hits",&tau1_byCombinedIsolationDeltaBetaCorr3Hits,"tau1_byCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits",&tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits,"tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits",&tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits,"tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3Hits",&tau1_byTightCombinedIsolationDeltaBetaCorr3Hits,"tau1_byTightCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",&tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03,"tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT",&tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT,"tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_rawMVArun2v1DBdR03oldDMwLT",&tau1_rawMVArun2v1DBdR03oldDMwLT,"tau1_rawMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau1_againstMuonLoose3",&tau1_againstMuonLoose3,"tau1_againstMuonLoose3/F");
   m_tree->Branch("tau1_againstMuonTight3",&tau1_againstMuonTight3,"tau1_againstMuonTight3/F");
   m_tree->Branch("tau1_againstElectronVLooseMVA6",&tau1_againstElectronVLooseMVA6,"tau1_againstElectronVLooseMVA6/F");
   m_tree->Branch("tau1_againstElectronLooseMVA6",&tau1_againstElectronLooseMVA6,"tau1_againstElectronLooseMVA6/F");
   m_tree->Branch("tau1_againstElectronMediumMVA6",&tau1_againstElectronMediumMVA6,"tau1_againstElectronMediumMVA6/F");
   m_tree->Branch("tau1_againstElectronTightMVA6",&tau1_againstElectronTightMVA6,"tau1_againstElectronTightMVA6/F");

   m_tree->Branch("tau2_pt",&tau2_pt,"tau2_pt/F");
   m_tree->Branch("tau2_eta",&tau2_eta,"tau2_eta/F");
   m_tree->Branch("tau2_phi",&tau2_phi,"tau2_phi/F");
   m_tree->Branch("tau2_E",&tau2_E,"tau2_E/F");
   m_tree->Branch("tau2_charge",&tau2_charge,"tau2_charge/I");
   m_tree->Branch("tau2_dxy",&tau2_dxy,"tau2_dxy/F");
   m_tree->Branch("tau2_dz",&tau2_dz,"tau2_dz/F");
   m_tree->Branch("tau2_decayModeFindingOldDMs",&tau2_decayModeFindingOldDMs,"tau2_decayModeFindingOldDMs/F");
   m_tree->Branch("tau2_decayModeFindingNewDMs",&tau2_decayModeFindingNewDMs,"tau2_decayModeFindingNewDMs/F");
   m_tree->Branch("tau2_byCombinedIsolationDeltaBetaCorr3Hits",&tau2_byCombinedIsolationDeltaBetaCorr3Hits,"tau2_byCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits",&tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits,"tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits",&tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits,"tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau2_byTightCombinedIsolationDeltaBetaCorr3Hits",&tau2_byTightCombinedIsolationDeltaBetaCorr3Hits,"tau2_byTightCombinedIsolationDeltaBetaCorr3Hits/F");
   m_tree->Branch("tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",&tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03,"tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03/F");
   m_tree->Branch("tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT",&tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT,"tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_rawMVArun2v1DBdR03oldDMwLT",&tau2_rawMVArun2v1DBdR03oldDMwLT,"tau2_rawMVArun2v1DBdR03oldDMwLT/F");
   m_tree->Branch("tau2_againstMuonLoose3",&tau2_againstMuonLoose3,"tau2_againstMuonLoose3/F");
   m_tree->Branch("tau2_againstMuonTight3",&tau2_againstMuonTight3,"tau2_againstMuonTight3/F");
   m_tree->Branch("tau2_againstElectronVLooseMVA6",&tau2_againstElectronVLooseMVA6,"tau2_againstElectronVLooseMVA6/F");
   m_tree->Branch("tau2_againstElectronLooseMVA6",&tau2_againstElectronLooseMVA6,"tau2_againstElectronLooseMVA6/F");
   m_tree->Branch("tau2_againstElectronMediumMVA6",&tau2_againstElectronMediumMVA6,"tau2_againstElectronMediumMVA6/F");
   m_tree->Branch("tau2_againstElectronTightMVA6",&tau2_againstElectronTightMVA6,"tau2_againstElectronTightMVA6/F");
   
   m_tree->Branch("jet1_pt",&jet1_pt,"jet1_pt/F");
   m_tree->Branch("jet1_eta",&jet1_eta,"jet1_eta/F");
   m_tree->Branch("jet1_phi",&jet1_phi,"jet1_phi/F");
   m_tree->Branch("jet1_E",&jet1_E,"jet1_E/F");
   m_tree->Branch("jet1_CSV",&jet1_CSV,"jet1_CSV/F");

   m_tree->Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
   m_tree->Branch("jet2_eta",&jet2_eta,"jet2_eta/F");
   m_tree->Branch("jet2_phi",&jet2_phi,"jet2_phi/F");
   m_tree->Branch("jet2_E",&jet2_E,"jet2_E/F");
   m_tree->Branch("jet2_CSV",&jet2_CSV,"jet2_CSV/F");

   m_tree->Branch("jet3_pt",&jet3_pt,"jet3_pt/F");
   m_tree->Branch("jet3_eta",&jet3_eta,"jet3_eta/F");
   m_tree->Branch("jet3_phi",&jet3_phi,"jet3_phi/F");
   m_tree->Branch("jet3_E",&jet3_E,"jet3_E/F");
   m_tree->Branch("jet3_CSV",&jet3_CSV,"jet3_CSV/F");

   m_tree->Branch("jet4_pt",&jet4_pt,"jet4_pt/F");
   m_tree->Branch("jet4_eta",&jet4_eta,"jet4_eta/F");
   m_tree->Branch("jet4_phi",&jet4_phi,"jet4_phi/F");
   m_tree->Branch("jet4_E",&jet4_E,"jet4_E/F");
   m_tree->Branch("jet4_CSV",&jet4_CSV,"jet4_CSV/F");
   
   m_tree->Branch("PFMET",&PFMET,"PFMET/F");
   m_tree->Branch("PFMETphi",&PFMETphi,"PFMETphi/F");
   m_tree->Branch("MHT",&MHT,"MHT/F");
   m_tree->Branch("metLD",&metLD,"metLD/F");
   m_tree->Branch("isGenMatched",&isGenMatched,"isGenMatched/I");
   
   m_tree->Branch("lep1_conePt",&lep1_conePt,"lep1_conePt/F");
   m_tree->Branch("lep2_conePt",&lep2_conePt,"lep2_conePt/F");
   m_tree->Branch("lep3_conePt",&lep3_conePt,"lep3_conePt/F");
   
   m_tree->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
   m_tree->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
   m_tree->Branch("mindr_lep3_jet",&mindr_lep3_jet,"mindr_lep3_jet/F");

   m_tree->Branch("mindr_tau1_jet",&mindr_tau1_jet,"mindr_tau1_jet/F");
   m_tree->Branch("mindr_tau2_jet",&mindr_tau2_jet,"mindr_tau2_jet/F");
   m_tree->Branch("mindr_tau3_jet",&mindr_tau3_jet,"mindr_tau3_jet/F");
   
   m_tree->Branch("avg_dr_jet",&avg_dr_jet,"avg_dr_jet/F");
   m_tree->Branch("avr_dr_lep_tau",&avr_dr_lep_tau,"avr_dr_lep_tau/F");
   m_tree->Branch("max_dr_lep_tau",&max_dr_lep_tau,"max_dr_lep_tau/F");
   m_tree->Branch("mindr_tau_jet",&mindr_tau_jet,"mindr_tau_jet/F");
   m_tree->Branch("min_dr_lep_tau",&min_dr_lep_tau,"min_dr_lep_tau/F");
   m_tree->Branch("min_dr_lep_jet",&min_dr_lep_jet,"min_dr_lep_jet/F");
   m_tree->Branch("dr_leps",&dr_leps,"dr_leps/F");
   m_tree->Branch("dr_taus",&dr_taus,"dr_taus/F");
   m_tree->Branch("dR_lep_tau_ss",&dR_lep_tau_ss,"dR_lep_tau_ss/F");
   m_tree->Branch("dr_lep1_tau",&dr_lep1_tau,"dr_lep1_tau/F");
   m_tree->Branch("dr_lep2_tau",&dr_lep2_tau,"dr_lep2_tau/F");
   m_tree->Branch("max_lep_eta",&max_lep_eta,"max_lep_eta/F");
   m_tree->Branch("mT_lep1",&mT_lep1,"mT_lep1/F");
   m_tree->Branch("mT_lep2",&mT_lep2,"mT_lep2/F");
   m_tree->Branch("mTauTauVis",&mTauTauVis,"mTauTauVis/F");
   m_tree->Branch("mTauTauVis1",&mTauTauVis1,"mTauTauVis1/F");
   m_tree->Branch("mTauTauVis2",&mTauTauVis2,"mTauTauVis2/F");
   m_tree->Branch("mbb",&mbb,"mbb/F");
   m_tree->Branch("mbb_loose",&mbb_loose,"mbb_loose/F");
   m_tree->Branch("cosThetaS_hadTau",&cosThetaS_hadTau,"cosThetaS_hadTau/F");
   m_tree->Branch("HTT",&HTT,"HTT/F");
   m_tree->Branch("HadTop_pt",&HadTop_pt,"HadTop_pt/F");
   m_tree->Branch("Hj_tagger",&Hj_tagger,"Hj_tagger/F");
   m_tree->Branch("nBJetLoose",&nBJetLoose,"nBJetLoose/I");
   
   m_tree->Branch("mvaOutput_plainKin_ttV",&mvaOutput_plainKin_ttV,"mvaOutput_plainKin_ttV/F");
   m_tree->Branch("mvaOutput_plainKin_ttbar",&mvaOutput_plainKin_ttbar,"mvaOutput_plainKin_ttbar/F");
   m_tree->Branch("mvaOutput_1l_2tau_HTT_SUM_VT",&mvaOutput_1l_2tau_HTT_SUM_VT,"mvaOutput_1l_2tau_HTT_SUM_VT/F");
   m_tree->Branch("mvaOutput_2l_2tau_plainKin_1B_VT",&mvaOutput_2l_2tau_plainKin_1B_VT,"mvaOutput_2l_2tau_plainKin_1B_VT/F");
   m_tree->Branch("mvaOutput_2l_2tau_plainKin_SUM_VT",&mvaOutput_2l_2tau_plainKin_SUM_VT,"mvaOutput_2l_2tau_plainKin_SUM_VT/F");
   m_tree->Branch("mvaOutput_2lss_ttV",&mvaOutput_2lss_ttV,"mvaOutput_2lss_ttV/F");
   m_tree->Branch("mvaOutput_2lss_ttbar",&mvaOutput_2lss_ttbar,"mvaOutput_2lss_ttbar/F");
   m_tree->Branch("mvaOutput_2lss_1tau_plainKin_ttbar",&mvaOutput_2lss_1tau_plainKin_ttbar,"mvaOutput_2lss_1tau_plainKin_ttbar/F");
   m_tree->Branch("mvaOutput_2lss_1tau_plainKin_ttV",&mvaOutput_2lss_1tau_plainKin_ttV,"mvaOutput_2lss_1tau_plainKin_ttV/F");
   m_tree->Branch("mvaOutput_2lss_1tau_plainKin_1B_M",&mvaOutput_2lss_1tau_plainKin_1B_M,"mvaOutput_2lss_1tau_plainKin_1B_M/F");
   m_tree->Branch("mvaOutput_2lss_1tau_plainKin_SUM_M",&mvaOutput_2lss_1tau_plainKin_SUM_M,"mvaOutput_2lss_1tau_plainKin_SUM_M/F");
   m_tree->Branch("mvaOutput_2lss_1tau_HTT_SUM_M",&mvaOutput_2lss_1tau_HTT_SUM_M,"mvaOutput_2lss_1tau_HTT_SUM_M/F");
   m_tree->Branch("mvaOutput_2lss_1tau_HTTMEM_SUM_M",&mvaOutput_2lss_1tau_HTTMEM_SUM_M,"mvaOutput_2lss_1tau_HTTMEM_SUM_M/F");
   m_tree->Branch("mvaOutput_3l_ttV",&mvaOutput_3l_ttV,"mvaOutput_3l_ttV/F");
   m_tree->Branch("mvaOutput_3l_ttbar",&mvaOutput_3l_ttbar,"mvaOutput_3l_ttbar/F");
   m_tree->Branch("mvaOutput_3l_1tau_plainKin_SUM_M",&mvaOutput_3l_1tau_plainKin_SUM_M,"mvaOutput_3l_1tau_plainKin_SUM_M/F");
   m_tree->Branch("mvaOutput_3l_1tau_plainKin_1B_M",&mvaOutput_3l_1tau_plainKin_1B_M,"mvaOutput_3l_1tau_plainKin_1B_M/F");
   
   m_tree->Branch("FR_weight",&FR_weight,"FR_weight/F");
   m_tree->Branch("triggerSF_weight",&triggerSF_weight,"triggerSF_weight/F");
   m_tree->Branch("leptonSF_weight",&leptonSF_weight,"leptonSF_weight/F");
   m_tree->Branch("tauSF_weight",&tauSF_weight,"tauSF_weight/F");
   m_tree->Branch("bTagSF_weight",&bTagSF_weight,"bTagSF_weight/F");
   m_tree->Branch("PU_weight",&PU_weight,"PU_weight/F");
   m_tree->Branch("MC_weight",&MC_weight,"MC_weight/F");
   
   m_tree->Branch("Integral_ttH",&Integral_ttH,"Integral_ttH/F");
   m_tree->Branch("Integral_ttZ",&Integral_ttZ,"Integral_ttZ/F");
   m_tree->Branch("Integral_ttZ_Zll",&Integral_ttZ_Zll,"Integral_ttZ_Zll/F");
   m_tree->Branch("Integral_ttbar",&Integral_ttbar,"Integral_ttbar/F");
   m_tree->Branch("integration_type",&integration_type,"integration_type/F");
   m_tree->Branch("memOutput_LR",&memOutput_LR,"memOutput_LR/F");
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

void Sync::fill()
{
   m_tree->Fill();
}
