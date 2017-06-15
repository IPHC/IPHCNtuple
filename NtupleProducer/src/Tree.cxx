#include "include/Tree.h"
#include "include/NtupleProducer.h"

Tree::Tree(TChain *ch,std::string fname,std::string treename)
{
    ch = new TChain(treename.c_str());

    std::ifstream infile;
    infile.open(fname.c_str());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        //	std::string fnameStr = "rfio://"+std::string(ifile);

        ch->Add(fnameStr.c_str());

        std::cout << "file: " << fnameStr << std::endl;
    }   
    infile.close();
    Init(ch);
}

Tree::~Tree()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t Tree::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

/*Long64_t Tree::LoadTree(Long64_t entry)
  {
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
  fCurrent = chain->GetTreeNumber();
  }
  return centry;
  }*/

void Tree::registerInputBranches(TChain *ch)
{
    ch->SetBranchStatus("*",1); 
    std::cout << "Successfully initialized input branches" << std::endl;
}

void Tree::Init(TChain *ch)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    trigger      = 0;
    trigger_pass = 0;
    trigger_name = 0;

    // Set object pointer
    nvertex = 0;
    pv_x = 0;
    pv_y = 0;
    pv_z = 0;
    pv_zError = 0;
    
    mc_pdfweights = 0;
    mc_pdfweightIds = 0;

    // ####################################
    // #   ____  _ _                      #
    // #  |  _ \(_) | ___   _   _ _ __    #
    // #  | |_) | | |/ _ \ | | | | '_ \   #
    // #  |  __/| | |  __/ | |_| | |_) |  #
    // #  |_|   |_|_|\___|  \__,_| .__/   #
    // #                         |_|      #
    // #                                  #                      
    // ####################################

    mc_pu_Nzpositions = 0;
    mc_pu_BunchCrossing = 0;
    mc_pu_zpositions = 0;
    mc_pu_sumpT_lowpT = 0;
    mc_pu_sumpT_highpT = 0;
    mc_pu_ntrks_lowpT = 0;
    mc_pu_ntrks_highpT = 0;

    // #################################################
    // #   _____ _           _                         #
    // #  | ____| | ___  ___| |_ _ __ ___  _ __  ___   #
    // #  |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|  #
    // #  | |___| |  __/ (__| |_| | | (_) | | | \__ \  #
    // #  |_____|_|\___|\___|\__|_|  \___/|_| |_|___/  #
    // #                                               #
    // #################################################

    el_pt = 0;
    el_eta = 0;
    el_phi = 0;
    el_m = 0;
    el_E = 0;
    el_looseCBId = 0;
    el_mediumCBId = 0;
    el_numberOfLostHits = 0;
    el_gsfTrack_PV_dxy = 0;
    el_gsfTrack_PV_dz = 0;
    el_ip3d = 0;
    el_ip3dErr = 0;
    el_miniIso = 0;
    el_miniIsoTTH = 0;

    el_id = 0;
    el_charge = 0;
    el_neutralHadronIso = 0;
    el_chargedHadronIso = 0;
    el_puChargedHadronIso = 0;
    el_ecalIso = 0;
    el_hcalIso = 0;
    el_particleIso = 0;
    el_photonIso = 0;
    el_trackIso = 0;
    el_isLoose = 0;
    el_isTight = 0;
    el_isRobustLoose = 0;
    el_isRobustTight = 0;
    el_isRobustHighEnergy = 0;
    el_vx = 0;
    el_vy = 0;
    el_vz = 0;
    el_isGsf = 0;
    el_dxy = 0;
    el_dz = 0;
    el_dxyError = 0;
    el_dzError = 0;
    el_mvaNonTrigV0 = 0;
    el_mvaNonTrigCat = 0;
    el_mvaPassMedium = 0;
    el_mvaPassTight = 0;
    el_numberOfHits = 0;
    el_pfIso_sumChargedHadronPt = 0;
    el_pfIso_sumNeutralHadronEt = 0;
    el_pfIso_sumPhotonEt = 0;
    el_pfIso_sumPUPt = 0;
    el_lepMVA = 0;
    el_lepMVA_miniRelIsoCharged = 0;
    el_lepMVA_miniRelIsoNeutral = 0;
    el_lepMVA_jetPtRelv2 = 0;
    el_lepMVA_neuRelIso = 0;
    el_lepMVA_chRelIso = 0;
    el_lepMVA_jetDR = 0;
    el_lepMVA_jetPtRatio = 0;
    el_lepMVA_jetBTagCSV = 0;
    el_lepMVA_sip3d = 0;
    el_lepMVA_dxy = 0;
    el_lepMVA_dz = 0;
    el_lepMVA_mvaId = 0;
    el_lepMVA_Moriond16 = 0; 
    el_lepMVA_eta = 0;
    el_lepMVA_jetNDauChargedMVASel = 0;         
    el_isGsfCtfScPixChargeConsistent = 0;
    el_isGsfScPixChargeConsistent = 0;
    el_passConversionVeto = 0;
    el_deltaEtaSuperClusterTrackAtVtx = 0;
    el_deltaPhiSuperClusterTrackAtVtx = 0;
    el_see = 0;
    el_hadronicOverEm = 0;
    el_scleta = 0;
    el_dB3D = 0;
    el_edB3D = 0;
    el_hasMatchedConversion = 0;
    el_sigmaIetaIeta = 0;
    el_superCluster_eta = 0;
    el_correctedEcalEnergy = 0;
    el_ecalEnergy = 0;
    el_eSuperClusterOverP = 0;
    el_trackMomentumError = 0;


    // ####################################
    // #   __  __                         #
    // #  |  \/  |_   _  ___  _ __  ___   #
    // #  | |\/| | | | |/ _ \| '_ \/ __|  #
    // #  | |  | | |_| | (_) | | | \__ \  #
    // #  |_|  |_|\__,_|\___/|_| |_|___/  #
    // #                                  #                                     
    // ####################################

    mu_pt = 0;
    mu_eta = 0;
    mu_phi = 0;
    mu_m = 0;
    mu_E = 0;
    mu_id = 0;
    mu_charge = 0;
    mu_ip3d = 0;
    mu_ip3dErr = 0;
    mu_miniIso = 0;
    mu_miniIsoTTH = 0;
    mu_isLooseMuon = 0;

    mu_neutralHadronIso = 0;
    mu_chargedHadronIso = 0;
    mu_ecalIso = 0;
    mu_hcalIso = 0;
    mu_photonIso = 0;
    mu_trackIso = 0;
    mu_isGlobalMuon     = 0;
    mu_isTrackerMuon    = 0;
    mu_isStandAloneMuon = 0;
    mu_isCaloMuon       = 0;
    mu_isPFMuon         = 0;
    mu_isMediumMuon     = 0;
    mu_isTightMuon      = 0;
    mu_vx = 0;
    mu_vy = 0;
    mu_vz = 0;
    mu_segmentCompatibility = 0;
    mu_hasGlobalTrack = 0;
    mu_globalTrack_dxy = 0;
    mu_globalTrack_dz = 0;
    mu_globalTrack_dxyError = 0;
    mu_globalTrack_dzError = 0;
    mu_globalTrack_normalizedChi2 = 0;
    mu_combinedQuality_chi2LocalPosition = 0;
    mu_combinedQuality_trkKink = 0;
    mu_hasInnerTrack = 0;
    mu_innerTrack_dxy = 0;
    mu_innerTrack_dz = 0;
    mu_innerTrack_PV_dxy = 0;
    mu_innerTrack_PV_dz = 0;
    mu_innerTrack_dxyError = 0;
    mu_innerTrack_dzError = 0;
    mu_innerTrack_validFraction = 0;
    mu_bestTrack_dxy = 0;
    mu_bestTrack_dz = 0;
    mu_bestTrack_dxyError = 0;
    mu_bestTrack_dzError  = 0;
    mu_bestTrack_pt       = 0;
    mu_bestTrack_ptError  = 0;
    mu_numberOfMatches    = 0;
    mu_numberOfValidMuonHits = 0;
    mu_pfIso03_sumChargedHadronPt = 0;
    mu_pfIso03_sumNeutralHadronEt = 0;
    mu_pfIso03_sumPhotonEt = 0;
    mu_pfIso03_sumPUPt = 0;
    mu_lepMVA = 0;
    mu_lepMVA_miniRelIsoCharged = 0;
    mu_lepMVA_miniRelIsoNeutral = 0;
    mu_lepMVA_jetPtRelv2 = 0;
    mu_lepMVA_neuRelIso = 0;
    mu_lepMVA_chRelIso = 0;
    mu_lepMVA_jetDR = 0;
    mu_lepMVA_jetPtRatio = 0;
    mu_lepMVA_jetBTagCSV = 0;
    mu_lepMVA_sip3d = 0;
    mu_lepMVA_dxy = 0;
    mu_lepMVA_dz = 0;
    mu_lepMVA_mvaId = 0;
    mu_lepMVA_Moriond16 = 0; 
    mu_lepMVA_eta = 0;
    mu_lepMVA_jetNDauChargedMVASel = 0;        
    mu_innerTrack_pt = 0;
    mu_innerTrack_ptError = 0;
    mu_dB3D = 0;
    mu_edB3D = 0;

    // #########################
    // #  _____                #
    // # |_   _|_ _ _   _ ___  #
    // #   | |/ _` | | | / __| #
    // #   | | (_| | |_| \__ \ #
    // #   |_|\__,_|\__,_|___/ #
    // #                       #
    // #########################

    tau_n = 0;
    tau_E = 0;
    tau_pt = 0;
    tau_eta = 0;
    tau_phi = 0;
    tau_m = 0;
    tau_dxy = 0;
    tau_dz = 0;
    tau_leadingTrackDxy = 0;
    tau_leadingTrackDz = 0;
    tau_charge = 0;
    tau_id = 0;
    tau_decayMode = 0;
    tau_hasLeadChargedHadrCand = 0;
    tau_leadingTrackPt = 0;
    tau_leadingTrackCharge = 0;
    tau_decayModeFindingOldDMs = 0;
    tau_decayModeFindingNewDMs = 0;
    tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
    tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
    tau_byLooseIsolationMVA3newDMwLT = 0;
    tau_byMediumIsolationMVA3newDMwLT = 0;
    tau_byTightIsolationMVA3newDMwLT = 0;
    
    tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03  = 0;
    tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = 0;
    tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03  = 0;
    tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT         = 0;
    tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT        = 0;
    tau_byTightIsolationMVArun2v1DBdR03oldDMwLT         = 0;
    tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT        = 0;
    
    tau_byCombinedIsolationDeltaBetaCorrRaw3Hits        = 0;
    tau_chargedIsoPtSum = 0;
    tau_neutralIsoPtSum = 0;
    tau_puCorrPtSum = 0;
    tau_againstMuonLoose3 = 0;
    tau_againstMuonTight3 = 0;
    //AC8X
    tau_againstElectronVLooseMVA6 = 0;
    tau_againstElectronLooseMVA6 = 0;
    tau_againstElectronMediumMVA6 = 0;
    tau_againstElectronTightMVA6 = 0;
    tau_pfEssential_jet_pt = 0;
    tau_pfEssential_jet_eta = 0;
    tau_pfEssential_jet_phi = 0;
    tau_pfEssential_jet_m = 0;
    tau_pfEssential_jetCorr_pt = 0;
    tau_pfEssential_jetCorr_eta = 0;
    tau_pfEssential_jetCorr_phi = 0;
    tau_pfEssential_jetCorr_m = 0;
    tau_pfEssential_hasSV = 0;
    tau_pfEssential_sv_x = 0;
    tau_pfEssential_sv_y = 0;
    tau_pfEssential_sv_z = 0;
    tau_pfEssential_flightLengthSig = 0;
    tau_pfEssential_dxy = 0;
    tau_pfEssential_dxy_error = 0;
    tau_pfEssential_dxy_Sig = 0;

    // ##########################
    // #       _      _         #
    // #      | | ___| |_ ___   #
    // #   _  | |/ _ \ __/ __|  #
    // #  | |_| |  __/ |_\__ \  #
    // #   \___/ \___|\__|___/  #
    // #                        #                        
    // ##########################

    jet_pt = 0;
    jet_eta = 0;
    jet_phi = 0;
    jet_m = 0;
    jet_E = 0;
    jet_qgtag = 0;
    jet_ntrk = 0;
    jet_CSVv2 = 0;
    jet_cMVAv2 = 0;
    jet_DeepCSVProbudsg = 0;
    jet_DeepCSVProbb = 0;
    jet_DeepCSVProbbb = 0;
    jet_DeepCSVProbc = 0;
    jet_DeepCSVProbcc = 0;
    jet_looseJetID = 0;
    jet_partonFlavour = 0;
    jet_hadronFlavour = 0;
    jet_neutralHadronEnergy = 0;
    jet_neutralEmEnergy = 0;
    jet_chargedHadronEnergy = 0;
    jet_chargedEmEnergy = 0;
    jet_electronEnergy = 0;
    jet_muonEnergy = 0;
    jet_photonEnergy = 0;
    jet_genJet_pt = 0;
    //jet_genJet_eta = 0;
    //jet_genJet_phi = 0;
    //jet_genJet_m = 0;
    jet_genJet_E = 0;
    //jet_genJet_status = 0;
    //jet_genJet_id = 0;
    jet_genParton_pt = 0;
    //jet_genParton_eta = 0;
    //jet_genParton_phi = 0;
    //jet_genParton_m = 0;
    jet_genParton_E = 0;
    //jet_genParton_status = 0;
    jet_genParton_id = 0;
    jet_pileupJetId = 0;

    // #########################################
    // #                                       #
    // #   __ _  ___ _ __     (_) ___| |_ ___  #
    // #  / _` |/ _ \ '_ \    | |/ _ \ __/ __| #
    // # | (_| |  __/ | | |   | |  __/ |_\__ \ #
    // #  \__, |\___|_| |_|  _/ |\___|\__|___/ #
    // #  |___/             |__/               #
    // #                                       #            
    // #########################################

    genJet_n = 0;
    genJet_pt = 0;
    genJet_eta = 0;
    genJet_phi = 0;
    genJet_m = 0;
    genJet_E = 0;

    metGen_px = 0;
    metGen_py = 0;
    metGen_pt = 0;
    metGen_phi = 0;
    metGen_sumet = 0;
    metGen_MuonEt = 0;

    // #####################
    // #   __  __  ____    #
    // #  |  \/  |/ ___|   #
    // #  | |\/| | |       #
    // #  | |  | | |___    #
    // #  |_|  |_|\____|   #
    // #                   #
    // #####################
    gen_n = 0;
    gen_pt = 0;
    gen_eta = 0;
    gen_phi = 0;
    gen_m = 0;
    gen_id = 0;
    gen_status = 0;
    gen_index = 0;
    gen_mother_index = 0;
    gen_PVz = 0;

    // ##################################################
    // #   __  __  ____     _____           _   _       #
    // #  |  \/  |/ ___|   |_   _| __ _   _| |_| |__    #
    // #  | |\/| | |         | || '__| | | | __| '_ \   #
    // #  | |  | | |___      | || |  | |_| | |_| | | |  #
    // #  |_|  |_|\____|     |_||_|   \__,_|\__|_| |_|  #
    // #                                                #
    // ##################################################

    mc_truth_h0_id = 0;
    mc_truth_h0W1_id = 0;
    mc_truth_h0Wl1_id = 0;
    mc_truth_h0Wtau1_id = 0;
    mc_truth_h0Wtaul1_id = 0;
    mc_truth_h0Wq11_id = 0;
    mc_truth_h0Wq21_id = 0;
    mc_truth_h0W2_id = 0;
    mc_truth_h0Wl2_id = 0;
    mc_truth_h0Wtau2_id = 0;
    mc_truth_h0Wtaul2_id = 0;
    mc_truth_h0Wq12_id = 0;
    mc_truth_h0Wq22_id = 0;
    mc_truth_h0Z1_id = 0;
    mc_truth_h0Zl11_id = 0;
    mc_truth_h0Zl21_id = 0;
    mc_truth_h0Zq11_id = 0;
    mc_truth_h0Zq21_id = 0;
    mc_truth_h0Z2_id = 0;
    mc_truth_h0Zl12_id = 0;
    mc_truth_h0Zl22_id = 0;
    mc_truth_h0Zq12_id = 0;
    mc_truth_h0Zq22_id = 0;
    mc_truth_h0tau1_id = 0;
    mc_truth_h0tau2_id = 0;
    mc_truth_h0taul1_id = 0;
    mc_truth_h0taul2_id = 0;
    mc_truth_h0b1_id = 0;
    mc_truth_h0b2_id = 0;

    mc_truth_t1_id = 0;
    mc_truth_t2_id = 0;
    mc_truth_tb1_id = 0;
    mc_truth_tb2_id = 0;
    mc_truth_tW1_id = 0;
    mc_truth_tW2_id = 0;
    mc_truth_tWl1_id = 0;
    mc_truth_tWl2_id = 0;
    mc_truth_tWtau1_id = 0;
    mc_truth_tWtau2_id = 0;  
    mc_truth_tWtaul1_id = 0;
    mc_truth_tWtaul2_id = 0;
    mc_truth_tWq11_id = 0;
    mc_truth_tWq21_id = 0;
    mc_truth_tWq12_id = 0;
    mc_truth_tWq22_id = 0;

    mc_truth_t_id = 0;
    mc_truth_tb_id = 0;
    mc_truth_tW_id = 0;
    mc_truth_tWl_id = 0;
    mc_truth_tWq1_id = 0;
    mc_truth_tWq2_id = 0;

    mc_truth_W_id = 0;
    mc_truth_Wl_id = 0;  
    mc_truth_Wtau_id = 0; 
    mc_truth_Wtaul_id = 0; 
    mc_truth_Wq1_id = 0;
    mc_truth_Wq2_id = 0;
    mc_truth_Z_id = 0;
    mc_truth_Zl1_id = 0;
    mc_truth_Zl2_id = 0;
    mc_truth_Ztau1_id = 0;
    mc_truth_Ztau2_id = 0;
    mc_truth_Ztaul1_id = 0;
    mc_truth_Ztaul2_id = 0;   
    mc_truth_Zq1_id = 0;
    mc_truth_Zq2_id = 0;
    mc_truth_gammal1_id = 0;
    mc_truth_gammal2_id = 0;
    mc_truth_gammatau1_id = 0;
    mc_truth_gammatau2_id = 0;
    mc_truth_gammataul1_id = 0;
    mc_truth_gammataul2_id = 0;
    mc_truth_gamma_id = 0;

    mc_truth_t1_pt = 0;
    mc_truth_t2_pt = 0;
    mc_truth_tb1_pt = 0;
    mc_truth_tb2_pt = 0;
    mc_truth_tW1_pt = 0;
    mc_truth_tW2_pt = 0;
    mc_truth_tWl1_pt = 0;
    mc_truth_tWl2_pt = 0;
    mc_truth_tWtau1_pt = 0;
    mc_truth_tWtau2_pt = 0;  
    mc_truth_tWtaul1_pt = 0;
    mc_truth_tWtaul2_pt = 0;
    mc_truth_tWq11_pt = 0;
    mc_truth_tWq21_pt = 0;
    mc_truth_tWq12_pt = 0;
    mc_truth_tWq22_pt = 0;

    mc_truth_t_pt = 0;
    mc_truth_tb_pt = 0;
    mc_truth_tW_pt = 0;
    mc_truth_tWl_pt = 0;
    mc_truth_tWq1_pt = 0;
    mc_truth_tWq2_pt = 0;

    mc_truth_h0_pt = 0;
    mc_truth_h0W1_pt = 0;
    mc_truth_h0Wl1_pt = 0;
    mc_truth_h0Wtau1_pt = 0;
    mc_truth_h0Wq11_pt = 0;
    mc_truth_h0Wq21_pt = 0;
    mc_truth_h0W2_pt = 0;
    mc_truth_h0Wl2_pt = 0;
    mc_truth_h0Wq12_pt = 0;
    mc_truth_h0Wq22_pt = 0;
    mc_truth_h0Z1_pt = 0;
    mc_truth_h0Zl11_pt = 0;
    mc_truth_h0Zl21_pt = 0;
    mc_truth_h0Zq11_pt = 0;
    mc_truth_h0Zq21_pt = 0;
    mc_truth_h0Z2_pt = 0;
    mc_truth_h0Zl12_pt = 0;
    mc_truth_h0Zl22_pt = 0;
    mc_truth_h0Zq12_pt = 0;
    mc_truth_h0Zq22_pt = 0;
    mc_truth_h0tau1_pt = 0;
    mc_truth_h0tau2_pt = 0;
    mc_truth_h0b1_pt = 0;
    mc_truth_h0b2_pt = 0;

    mc_truth_W_pt = 0;  
    mc_truth_Wl_pt = 0;  
    mc_truth_Wtau_pt = 0; 
    mc_truth_Wtaul_pt = 0; 
    mc_truth_Wq1_pt = 0; 
    mc_truth_Wq2_pt = 0;
    mc_truth_Z_pt = 0;
    mc_truth_Zl1_pt = 0;
    mc_truth_Zl2_pt = 0; 
    mc_truth_Ztau1_pt = 0;
    mc_truth_Ztau2_pt = 0;
    mc_truth_Ztaul1_pt = 0;
    mc_truth_Ztaul2_pt = 0;  
    mc_truth_Zq1_pt = 0;
    mc_truth_Zq2_pt = 0;
    mc_truth_gammal1_pt = 0;
    mc_truth_gammal2_pt = 0;
    mc_truth_gammatau1_pt = 0;
    mc_truth_gammatau2_pt= 0;
    mc_truth_gammataul1_pt = 0;
    mc_truth_gammataul2_pt= 0;
    mc_truth_gamma_pt= 0;

    mc_truth_t1_eta = 0;
    mc_truth_t2_eta = 0;
    mc_truth_tb1_eta = 0;
    mc_truth_tb2_eta = 0;
    mc_truth_tW1_eta = 0;
    mc_truth_tW2_eta = 0;
    mc_truth_tWl1_eta = 0;
    mc_truth_tWl2_eta = 0;
    mc_truth_tWtau1_eta = 0;
    mc_truth_tWtau2_eta = 0;  
    mc_truth_tWtaul1_eta = 0;
    mc_truth_tWtaul2_eta = 0;
    mc_truth_tWq11_eta = 0;
    mc_truth_tWq21_eta = 0;
    mc_truth_tWq12_eta = 0;
    mc_truth_tWq22_eta = 0;

    mc_truth_t_eta = 0;
    mc_truth_tb_eta = 0;
    mc_truth_tW_eta = 0;
    mc_truth_tWl_eta = 0;
    mc_truth_tWq1_eta = 0;
    mc_truth_tWq2_eta = 0;

    mc_truth_h0_eta = 0;
    mc_truth_h0W1_eta = 0;
    mc_truth_h0Wl1_eta = 0; 
    mc_truth_h0Wtau1_eta = 0;
    mc_truth_h0Wtaul1_eta = 0;
    mc_truth_h0Wq11_eta = 0;
    mc_truth_h0Wq21_eta = 0;
    mc_truth_h0W2_eta = 0;
    mc_truth_h0Wl2_eta = 0; 
    mc_truth_h0Wtau2_eta = 0;
    mc_truth_h0Wtaul2_eta = 0;
    mc_truth_h0Wq12_eta = 0;
    mc_truth_h0Wq22_eta = 0;
    mc_truth_h0Z1_eta = 0;
    mc_truth_h0Zl11_eta = 0;
    mc_truth_h0Zl21_eta = 0;
    mc_truth_h0Zq11_eta = 0;
    mc_truth_h0Zq21_eta = 0;
    mc_truth_h0Z2_eta = 0;
    mc_truth_h0Zl12_eta = 0;
    mc_truth_h0Zl22_eta = 0;
    mc_truth_h0Zq12_eta = 0;
    mc_truth_h0Zq22_eta = 0;
    mc_truth_h0tau1_eta = 0;
    mc_truth_h0tau2_eta = 0;
    mc_truth_h0b1_eta = 0;
    mc_truth_h0b2_eta = 0;

    mc_truth_W_eta = 0;
    mc_truth_Wl_eta = 0;   
    mc_truth_Wtau_eta = 0; 
    mc_truth_Wtaul_eta = 0; 
    mc_truth_Wq1_eta = 0;  
    mc_truth_Wq2_eta = 0;
    mc_truth_Z_eta = 0;
    mc_truth_Zl1_eta = 0;
    mc_truth_Zl2_eta = 0;
    mc_truth_Ztau1_eta = 0;
    mc_truth_Ztau2_eta = 0;
    mc_truth_Ztaul1_eta = 0;
    mc_truth_Ztaul2_eta = 0;  
    mc_truth_Zq1_eta = 0;
    mc_truth_Zq2_eta = 0;
    mc_truth_gammal1_eta = 0;
    mc_truth_gammal2_eta = 0;
    mc_truth_gammatau1_eta = 0;
    mc_truth_gammatau2_eta= 0;
    mc_truth_gammataul1_eta = 0;
    mc_truth_gammataul2_eta= 0;
    mc_truth_gamma_eta= 0;

    mc_truth_t1_phi = 0;
    mc_truth_t2_phi = 0;
    mc_truth_tb1_phi = 0;
    mc_truth_tb2_phi = 0;
    mc_truth_tW1_phi = 0;
    mc_truth_tW2_phi = 0;
    mc_truth_tWtau1_phi = 0;
    mc_truth_tWtau2_phi = 0;  
    mc_truth_tWtaul1_phi = 0;
    mc_truth_tWtaul2_phi = 0;
    mc_truth_tWl1_phi = 0;
    mc_truth_tWl2_phi = 0;
    mc_truth_tWq11_phi = 0;
    mc_truth_tWq21_phi = 0;
    mc_truth_tWq12_phi = 0;
    mc_truth_tWq22_phi = 0;

    mc_truth_t_phi = 0;
    mc_truth_tb_phi = 0;
    mc_truth_tW_phi = 0;
    mc_truth_tWl_phi = 0;
    mc_truth_tWq1_phi = 0;
    mc_truth_tWq2_phi = 0;

    mc_truth_h0_phi = 0;
    mc_truth_h0W1_phi = 0;
    mc_truth_h0Wl1_phi = 0;
    mc_truth_h0Wtau1_phi = 0;
    mc_truth_h0Wtaul1_phi = 0;
    mc_truth_h0Wq11_phi = 0;
    mc_truth_h0Wq21_phi = 0;
    mc_truth_h0W2_phi = 0;
    mc_truth_h0Wl2_phi = 0; 
    mc_truth_h0Wtau2_phi = 0;
    mc_truth_h0Wtaul2_phi = 0;
    mc_truth_h0Wq12_phi = 0;
    mc_truth_h0Wq22_phi = 0;
    mc_truth_h0Z1_phi = 0;
    mc_truth_h0Zl11_phi = 0;
    mc_truth_h0Zl21_phi = 0;
    mc_truth_h0Zq11_phi = 0;
    mc_truth_h0Zq21_phi = 0;
    mc_truth_h0Z2_phi = 0;
    mc_truth_h0Zl12_phi = 0;
    mc_truth_h0Zl22_phi = 0;
    mc_truth_h0Zq12_phi = 0;
    mc_truth_h0Zq22_phi = 0;
    mc_truth_h0tau1_phi = 0;
    mc_truth_h0tau2_phi = 0;
    mc_truth_h0b1_phi = 0;
    mc_truth_h0b2_phi = 0;

    mc_truth_W_phi = 0;
    mc_truth_Wl_phi = 0; 
    mc_truth_Wtau_phi = 0;
    mc_truth_Wtaul_phi = 0;
    mc_truth_Wq1_phi = 0;
    mc_truth_Wq2_phi = 0;
    mc_truth_Z_phi = 0;
    mc_truth_Zl1_phi = 0;
    mc_truth_Zl2_phi = 0;
    mc_truth_Ztau1_phi = 0;
    mc_truth_Ztau2_phi = 0;
    mc_truth_Ztaul1_phi = 0;
    mc_truth_Ztaul2_phi = 0;  
    mc_truth_Zq1_phi = 0;
    mc_truth_Zq2_phi = 0;
    mc_truth_gammal1_phi = 0;
    mc_truth_gammal2_phi = 0;
    mc_truth_gammatau1_phi = 0;
    mc_truth_gammatau2_phi= 0;
    mc_truth_gammataul1_phi = 0;
    mc_truth_gammataul2_phi= 0;
    mc_truth_gamma_phi= 0;

    mc_truth_t1_E = 0;
    mc_truth_t2_E = 0;
    mc_truth_tb1_E = 0;
    mc_truth_tb2_E = 0;
    mc_truth_tW1_E = 0;
    mc_truth_tW2_E = 0;
    mc_truth_tWl1_E = 0;
    mc_truth_tWl2_E = 0; 
    mc_truth_tWtau1_E = 0;
    mc_truth_tWtau2_E = 0;  
    mc_truth_tWtaul1_E = 0;
    mc_truth_tWtaul2_E = 0;
    mc_truth_tWq11_E = 0;
    mc_truth_tWq21_E = 0;
    mc_truth_tWq12_E = 0;
    mc_truth_tWq22_E = 0;

    mc_truth_t_E = 0;
    mc_truth_tb_E = 0;
    mc_truth_tW_E = 0;
    mc_truth_tWl_E = 0;
    mc_truth_tWq1_E = 0;
    mc_truth_tWq2_E = 0;

    mc_truth_h0_E = 0;
    mc_truth_h0W1_E = 0;
    mc_truth_h0Wl1_E = 0;
    mc_truth_h0Wtau1_E = 0;
    mc_truth_h0Wtaul1_E = 0;
    mc_truth_h0Wq11_E = 0;
    mc_truth_h0Wq21_E = 0;
    mc_truth_h0W2_E = 0;
    mc_truth_h0Wl2_E = 0; 
    mc_truth_h0Wtau2_E = 0;
    mc_truth_h0Wtaul2_E = 0;
    mc_truth_h0Wq12_E = 0;
    mc_truth_h0Wq22_E = 0;
    mc_truth_h0Z1_E = 0;
    mc_truth_h0Zl11_E = 0;
    mc_truth_h0Zl21_E = 0;
    mc_truth_h0Zq11_E = 0;
    mc_truth_h0Zq21_E = 0;
    mc_truth_h0Z2_E = 0;
    mc_truth_h0Zl12_E = 0;
    mc_truth_h0Zl22_E = 0;
    mc_truth_h0Zq12_E = 0;
    mc_truth_h0Zq22_E = 0;
    mc_truth_h0tau1_E = 0;
    mc_truth_h0tau2_E = 0;
    mc_truth_h0b1_E = 0;
    mc_truth_h0b2_E = 0;

    mc_truth_W_E = 0;
    mc_truth_Wl_E = 0;  
    mc_truth_Wtau_E = 0;
    mc_truth_Wtaul_E = 0; 
    mc_truth_Wq1_E = 0;
    mc_truth_Wq2_E = 0;
    mc_truth_Z_E = 0;
    mc_truth_Zl1_E = 0;
    mc_truth_Zl2_E = 0; 
    mc_truth_Ztau1_E = 0;
    mc_truth_Ztau2_E = 0;
    mc_truth_Ztaul1_E = 0;
    mc_truth_Ztaul2_E = 0;  
    mc_truth_Zq1_E = 0;
    mc_truth_Zq2_E = 0;
    mc_truth_gammal1_E = 0;
    mc_truth_gammal2_E = 0;
    mc_truth_gammatau1_E = 0;
    mc_truth_gammatau2_E = 0;
    mc_truth_gammataul1_E = 0;
    mc_truth_gammataul2_E = 0;
    mc_truth_gamma_E = 0;

    triggerobject_n = 0;
    triggerobject_pt = 0;
    triggerobject_eta = 0;
    triggerobject_phi = 0;  
    triggerobject_collection = 0;
    triggerobject_filterIds_n = 0;
    triggerobject_filterLabels_n = 0;
    triggerobject_pathNamesAll_n = 0;

    triggerobject_pathNamesAll = 0;
    triggerobject_pathNamesAll_isL3 = 0;
    triggerobject_pathNamesAll_isLF = 0;
    triggerobject_pathNamesAll_isBoth = 0;
    triggerobject_pathNamesAll_isNone = 0;

    triggerobject_filterIds = 0;
    triggerobject_filterLabels = 0;


    // Set branch addresses and branch pointers
    if (!ch) return;
    fChain = ch;
    fCurrent = -1;
    //   fChain->SetMakeClass(1);

    fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
    fChain->SetBranchAddress("ev_id", &ev_id, &b_ev_id);
    fChain->SetBranchAddress("ev_lumi", &ev_lumi, &b_ev_lumi);
    fChain->SetBranchAddress("ev_rho", &ev_rho, &b_ev_rho);

    fChain->SetBranchAddress("trigger",      &trigger,      &b_trigger     );
    fChain->SetBranchAddress("trigger_pass", &trigger_pass, &b_trigger_pass);
    fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);

    fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
    fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
    fChain->SetBranchAddress("met_sumet", &met_sumet, &b_met_sumet);
    fChain->SetBranchAddress("met_cov00", &met_cov00, &b_met_cov00);
    fChain->SetBranchAddress("met_cov01", &met_cov01, &b_met_cov01);
    fChain->SetBranchAddress("met_cov10", &met_cov10, &b_met_cov10);
    fChain->SetBranchAddress("met_cov11", &met_cov11, &b_met_cov11);

    fChain->SetBranchAddress("metNoHF_pt",    &metNoHF_pt,    &b_metNoHF_pt);
    fChain->SetBranchAddress("metNoHF_phi",   &metNoHF_phi,   &b_metNoHF_phi);
    fChain->SetBranchAddress("metNoHF_sumet", &metNoHF_sumet, &b_metNoHF_sumet);

    fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
    fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
    fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
    fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
    fChain->SetBranchAddress("pv_zError", &pv_zError, &b_pv_zError);

    fChain->SetBranchAddress("weight_scale_muF0p5", &weight_scale_muF0p5, &b_weight_scale_muF0p5);
    fChain->SetBranchAddress("weight_scale_muF2",   &weight_scale_muF2,   &b_weight_scale_muF2  );
    fChain->SetBranchAddress("weight_scale_muR0p5", &weight_scale_muR0p5, &b_weight_scale_muR0p5);
    fChain->SetBranchAddress("weight_scale_muR2",   &weight_scale_muR2,   &b_weight_scale_muR2  );
    fChain->SetBranchAddress("mc_id", &mc_id, &b_mc_id);
    fChain->SetBranchAddress("mc_f1", &mc_f1, &b_mc_f1);
    fChain->SetBranchAddress("mc_f2", &mc_f2, &b_mc_f2);
    fChain->SetBranchAddress("mc_x1", &mc_x1, &b_mc_x1);
    fChain->SetBranchAddress("mc_x2", &mc_x2, &b_mc_x2);
    fChain->SetBranchAddress("mc_scale", &mc_scale, &b_mc_scale);
    fChain->SetBranchAddress("mc_ptHat", &mc_ptHat, &b_mc_ptHat);
    fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight); 
    fChain->SetBranchAddress("mc_pdfweights", &mc_pdfweights, &b_mc_pdfweights);
    fChain->SetBranchAddress("mc_pdfweightIds", &mc_pdfweightIds, &b_mc_pdfweightIds);
   
    fChain->SetBranchAddress("mc_pu_intime_NumInt", &mc_pu_intime_NumInt, &b_mc_pu_intime_NumInt);
    fChain->SetBranchAddress("mc_pu_trueNumInt", &mc_pu_trueNumInt, &b_mc_pu_trueNumInt);
    fChain->SetBranchAddress("mc_pu_before_npu", &mc_pu_before_npu, &b_mc_pu_before_npu);
    fChain->SetBranchAddress("mc_pu_after_npu", &mc_pu_after_npu, &b_mc_pu_after_npu);
    fChain->SetBranchAddress("mc_pu_Npvi", &mc_pu_Npvi, &b_mc_pu_Npvi);
    fChain->SetBranchAddress("mc_pu_Nzpositions", &mc_pu_Nzpositions, &b_mc_pu_Nzpositions);
    fChain->SetBranchAddress("mc_pu_BunchCrossing", &mc_pu_BunchCrossing, &b_mc_pu_BunchCrossing);
    fChain->SetBranchAddress("mc_pu_zpositions", &mc_pu_zpositions, &b_mc_pu_zpositions);
    fChain->SetBranchAddress("mc_pu_sumpT_lowpT", &mc_pu_sumpT_lowpT, &b_mc_pu_sumpT_lowpT);
    fChain->SetBranchAddress("mc_pu_sumpT_highpT", &mc_pu_sumpT_highpT, &b_mc_pu_sumpT_highpT);
    fChain->SetBranchAddress("mc_pu_ntrks_lowpT", &mc_pu_ntrks_lowpT, &b_mc_pu_ntrks_lowpT);
    fChain->SetBranchAddress("mc_pu_ntrks_highpT", &mc_pu_ntrks_highpT, &b_mc_pu_ntrks_highpT);

    // #################################################
    // #   _____ _           _                         #
    // #  | ____| | ___  ___| |_ _ __ ___  _ __  ___   #
    // #  |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|  #
    // #  | |___| |  __/ (__| |_| | | (_) | | | \__ \  #
    // #  |_____|_|\___|\___|\__|_|  \___/|_| |_|___/  #
    // #                                               #
    // #################################################

    fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
    fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
    fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
    fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
    fChain->SetBranchAddress("el_m", &el_m, &b_el_m);
    fChain->SetBranchAddress("el_E", &el_E, &b_el_E);
    fChain->SetBranchAddress("el_looseCBId", &el_looseCBId, &b_el_looseCBId);
    fChain->SetBranchAddress("el_mediumCBId", &el_mediumCBId, &b_el_mediumCBId);
    fChain->SetBranchAddress("el_numberOfLostHits", &el_numberOfLostHits, &b_el_numberOfLostHits);
    fChain->SetBranchAddress("el_gsfTrack_PV_dxy", &el_gsfTrack_PV_dxy, &b_el_gsfTrack_PV_dxy);
    fChain->SetBranchAddress("el_gsfTrack_PV_dz", &el_gsfTrack_PV_dz, &b_el_gsfTrack_PV_dz);
    fChain->SetBranchAddress("el_ip3d", &el_ip3d, &b_el_ip3d);
    fChain->SetBranchAddress("el_ip3dErr", &el_ip3dErr, &b_el_ip3dErr);
    fChain->SetBranchAddress("el_miniIso", &el_miniIso, &b_el_miniIso);
    fChain->SetBranchAddress("el_miniIsoTTH", &el_miniIsoTTH, &b_el_miniIsoTTH);

    fChain->SetBranchAddress("el_id", &el_id, &b_el_id);
    fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
    fChain->SetBranchAddress("el_neutralHadronIso", &el_neutralHadronIso, &b_el_neutralHadronIso);
    fChain->SetBranchAddress("el_chargedHadronIso", &el_chargedHadronIso, &b_el_chargedHadronIso);
    fChain->SetBranchAddress("el_puChargedHadronIso", &el_puChargedHadronIso, &b_el_puChargedHadronIso);
    fChain->SetBranchAddress("el_ecalIso", &el_ecalIso, &b_el_ecalIso);
    fChain->SetBranchAddress("el_hcalIso", &el_hcalIso, &b_el_hcalIso);
    fChain->SetBranchAddress("el_particleIso", &el_particleIso, &b_el_particleIso);
    fChain->SetBranchAddress("el_photonIso", &el_photonIso, &b_el_photonIso);
    fChain->SetBranchAddress("el_trackIso", &el_trackIso, &b_el_trackIso);
    //fChain->SetBranchAddress("el_isLoose", &el_isLoose, &b_el_isLoose);
    //fChain->SetBranchAddress("el_isTight", &el_isTight, &b_el_isTight);
    //fChain->SetBranchAddress("el_isRobustLoose", &el_isRobustLoose, &b_el_isRobustLoose);
    //fChain->SetBranchAddress("el_isRobustTight", &el_isRobustTight, &b_el_isRobustTight);
    //fChain->SetBranchAddress("el_isRobustHighEnergy", &el_isRobustHighEnergy, &b_el_isRobustHighEnergy);
    fChain->SetBranchAddress("el_vx", &el_vx, &b_el_vx);
    fChain->SetBranchAddress("el_vy", &el_vy, &b_el_vy);
    fChain->SetBranchAddress("el_vz", &el_vz, &b_el_vz);
    //fChain->SetBranchAddress("el_isGsf", &el_isGsf, &b_el_isGsf);
    //fChain->SetBranchAddress("el_dxy", &el_dxy, &b_el_dxy);
    //fChain->SetBranchAddress("el_dz", &el_dz, &b_el_dz);
    //fChain->SetBranchAddress("el_dxyError", &el_dxyError, &b_el_dxyError);
    //fChain->SetBranchAddress("el_dzError", &el_dzError, &b_el_dzError);
    fChain->SetBranchAddress("el_mvaNonTrigV0", &el_mvaNonTrigV0, &b_el_mvaNonTrigV0);
    fChain->SetBranchAddress("el_mvaNonTrigCat", &el_mvaNonTrigCat, &b_el_mvaNonTrigCat);
    //fChain->SetBranchAddress("el_mvaPassMedium", &el_mvaPassMedium, &b_el_mvaPassMedium);
    //fChain->SetBranchAddress("el_mvaPassTight", &el_mvaPassTight, &b_el_mvaPassTight);
    fChain->SetBranchAddress("el_numberOfHits", &el_numberOfHits, &b_el_numberOfHits);
    fChain->SetBranchAddress("el_pfIso_sumChargedHadronPt", &el_pfIso_sumChargedHadronPt, &b_el_pfIso_sumChargedHadronPt);
    fChain->SetBranchAddress("el_pfIso_sumNeutralHadronEt", &el_pfIso_sumNeutralHadronEt, &b_el_pfIso_sumNeutralHadronEt);
    fChain->SetBranchAddress("el_pfIso_sumPhotonEt", &el_pfIso_sumPhotonEt, &b_el_pfIso_sumPhotonEt);
    fChain->SetBranchAddress("el_pfIso_sumPUPt", &el_pfIso_sumPUPt, &b_el_pfIso_sumPUPt);
    fChain->SetBranchAddress("el_lepMVA", &el_lepMVA, &b_el_lepMVA);
    fChain->SetBranchAddress("el_lepMVA_miniRelIsoCharged", &el_lepMVA_miniRelIsoCharged, &b_el_lepMVA_miniRelIsoCharged);
    fChain->SetBranchAddress("el_lepMVA_miniRelIsoNeutral", &el_lepMVA_miniRelIsoNeutral, &b_el_lepMVA_miniRelIsoNeutral);
    fChain->SetBranchAddress("el_lepMVA_jetPtRelv2", &el_lepMVA_jetPtRelv2, &b_el_lepMVA_jetPtRelv2);
    fChain->SetBranchAddress("el_lepMVA_neuRelIso", &el_lepMVA_neuRelIso, &b_el_lepMVA_neuRelIso);
    fChain->SetBranchAddress("el_lepMVA_chRelIso", &el_lepMVA_chRelIso, &b_el_lepMVA_chRelIso);
    fChain->SetBranchAddress("el_lepMVA_jetDR", &el_lepMVA_jetDR, &b_el_lepMVA_jetDR);
    fChain->SetBranchAddress("el_lepMVA_jetPtRatio", &el_lepMVA_jetPtRatio, &b_el_lepMVA_jetPtRatio);
    fChain->SetBranchAddress("el_lepMVA_jetBTagCSV", &el_lepMVA_jetBTagCSV, &b_el_lepMVA_jetBTagCSV);
    fChain->SetBranchAddress("el_lepMVA_sip3d", &el_lepMVA_sip3d, &b_el_lepMVA_sip3d);
    fChain->SetBranchAddress("el_lepMVA_dxy", &el_lepMVA_dxy, &b_el_lepMVA_dxy);
    fChain->SetBranchAddress("el_lepMVA_dz", &el_lepMVA_dz, &b_el_lepMVA_dz);
    fChain->SetBranchAddress("el_lepMVA_mvaId", &el_lepMVA_mvaId, &b_el_lepMVA_mvaId);
    fChain->SetBranchAddress("el_lepMVA_Moriond16", &el_lepMVA_Moriond16, &b_el_lepMVA_Moriond16);
    fChain->SetBranchAddress("el_lepMVA_eta", &el_lepMVA_eta, &b_el_lepMVA_eta);
    fChain->SetBranchAddress("el_lepMVA_jetNDauChargedMVASel", &el_lepMVA_jetNDauChargedMVASel, &b_el_lepMVA_jetNDauChargedMVASel);   
    fChain->SetBranchAddress("el_isGsfCtfScPixChargeConsistent", &el_isGsfCtfScPixChargeConsistent, &b_el_isGsfCtfScPixChargeConsistent);
    fChain->SetBranchAddress("el_isGsfScPixChargeConsistent", &el_isGsfScPixChargeConsistent, &b_el_isGsfScPixChargeConsistent);
    fChain->SetBranchAddress("el_passConversionVeto", &el_passConversionVeto, &b_el_passConversionVeto);
    fChain->SetBranchAddress("el_deltaEtaSuperClusterTrackAtVtx", &el_deltaEtaSuperClusterTrackAtVtx, &b_el_deltaEtaSuperClusterTrackAtVtx);
    fChain->SetBranchAddress("el_deltaPhiSuperClusterTrackAtVtx", &el_deltaPhiSuperClusterTrackAtVtx, &b_el_deltaPhiSuperClusterTrackAtVtx);
    fChain->SetBranchAddress("el_see", &el_see, &b_el_see);
    fChain->SetBranchAddress("el_hadronicOverEm", &el_hadronicOverEm, &b_el_hadronicOverEm);
    //fChain->SetBranchAddress("el_scleta", &el_scleta, &b_el_scleta);
    fChain->SetBranchAddress("el_dB3D", &el_dB3D, &b_el_dB3D);
    fChain->SetBranchAddress("el_edB3D", &el_edB3D, &b_el_edB3D);
    fChain->SetBranchAddress("el_hasMatchedConversion", &el_hasMatchedConversion, &b_el_hasMatchedConversion);
    fChain->SetBranchAddress("el_sigmaIetaIeta", &el_sigmaIetaIeta, &b_el_sigmaIetaIeta);
    fChain->SetBranchAddress("el_superCluster_eta", &el_superCluster_eta, &b_el_superCluster_eta);
    fChain->SetBranchAddress("el_correctedEcalEnergy", &el_correctedEcalEnergy, &b_el_correctedEcalEnergy);
    fChain->SetBranchAddress("el_ecalEnergy", &el_ecalEnergy, &b_el_ecalEnergy);
    fChain->SetBranchAddress("el_eSuperClusterOverP", &el_eSuperClusterOverP, &b_el_eSuperClusterOverP);
    fChain->SetBranchAddress("el_trackMomentumError", &el_trackMomentumError, &b_el_trackMomentumError);

    // ####################################
    // #   __  __                         #
    // #  |  \/  |_   _  ___  _ __  ___   #
    // #  | |\/| | | | |/ _ \| '_ \/ __|  #
    // #  | |  | | |_| | (_) | | | \__ \  #
    // #  |_|  |_|\__,_|\___/|_| |_|___/  #
    // #                                  #                                     
    // ####################################

    fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
    fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
    fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
    fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
    fChain->SetBranchAddress("mu_m", &mu_m, &b_mu_m);
    fChain->SetBranchAddress("mu_E", &mu_E, &b_mu_E);
    fChain->SetBranchAddress("mu_id", &mu_id, &b_mu_id);
    fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
    fChain->SetBranchAddress("mu_ip3d", &mu_ip3d, &b_mu_ip3d);
    fChain->SetBranchAddress("mu_ip3dErr", &mu_ip3dErr, &b_mu_ip3dErr);
    fChain->SetBranchAddress("mu_miniIso", &mu_miniIso, &b_mu_miniIso);
    fChain->SetBranchAddress("mu_miniIsoTTH", &mu_miniIsoTTH, &b_mu_miniIsoTTH);
    fChain->SetBranchAddress("mu_isLooseMuon", &mu_isLooseMuon, &b_mu_isLooseMuon);

    fChain->SetBranchAddress("mu_neutralHadronIso", &mu_neutralHadronIso, &b_mu_neutralHadronIso);
    fChain->SetBranchAddress("mu_chargedHadronIso", &mu_chargedHadronIso, &b_mu_chargedHadronIso);
    fChain->SetBranchAddress("mu_ecalIso", &mu_ecalIso, &b_mu_ecalIso);
    fChain->SetBranchAddress("mu_hcalIso", &mu_hcalIso, &b_mu_hcalIso);
    fChain->SetBranchAddress("mu_photonIso", &mu_photonIso, &b_mu_photonIso);
    fChain->SetBranchAddress("mu_trackIso", &mu_trackIso, &b_mu_trackIso);
    fChain->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon, &b_mu_isGlobalMuon);
    fChain->SetBranchAddress("mu_isTrackerMuon", &mu_isTrackerMuon, &b_mu_isTrackerMuon);
    fChain->SetBranchAddress("mu_isStandAloneMuon", &mu_isStandAloneMuon, &b_mu_isStandAloneMuon);
    fChain->SetBranchAddress("mu_isCaloMuon", &mu_isCaloMuon, &b_mu_isCaloMuon);
    fChain->SetBranchAddress("mu_isPFMuon", &mu_isPFMuon, &b_mu_isPFMuon);
    fChain->SetBranchAddress("mu_isMediumMuon", &mu_isMediumMuon, &b_mu_isMediumMuon );
    fChain->SetBranchAddress("mu_isTightMuon",  &mu_isTightMuon,  &b_mu_isTightMuon  );
    fChain->SetBranchAddress("mu_vx", &mu_vx, &b_mu_vx);
    fChain->SetBranchAddress("mu_vy", &mu_vy, &b_mu_vy);
    fChain->SetBranchAddress("mu_vz", &mu_vz, &b_mu_vz);
    fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
    //fChain->SetBranchAddress("mu_hasGlobalTrack", &mu_hasGlobalTrack, &b_mu_hasGlobalTrack);
    //fChain->SetBranchAddress("mu_globalTrack_dxy", &mu_globalTrack_dxy, &b_mu_globalTrack_dxy);
    //fChain->SetBranchAddress("mu_globalTrack_dz", &mu_globalTrack_dz, &b_mu_globalTrack_dz);
    //fChain->SetBranchAddress("mu_globalTrack_dxyError", &mu_globalTrack_dxyError, &b_mu_globalTrack_dxyError);
    //fChain->SetBranchAddress("mu_globalTrack_dzError", &mu_globalTrack_dzError, &b_mu_globalTrack_dzError);
    //fChain->SetBranchAddress("mu_globalTrack_normalizedChi2", &mu_globalTrack_normalizedChi2, &b_mu_globalTrack_normalizedChi2);
    fChain->SetBranchAddress("mu_combinedQuality_chi2LocalPosition", &mu_combinedQuality_chi2LocalPosition, &b_mu_combinedQuality_chi2LocalPosition);
    fChain->SetBranchAddress("mu_combinedQuality_trkKink", &mu_combinedQuality_trkKink, &b_mu_combinedQuality_trkKink);
    //fChain->SetBranchAddress("mu_hasInnerTrack", &mu_hasInnerTrack, &b_mu_hasInnerTrack);
    //fChain->SetBranchAddress("mu_innerTrack_dxy", &mu_innerTrack_dxy, &b_mu_innerTrack_dxy);
    //fChain->SetBranchAddress("mu_innerTrack_dz", &mu_innerTrack_dz, &b_mu_innerTrack_dz);
    fChain->SetBranchAddress("mu_innerTrack_PV_dxy", &mu_innerTrack_PV_dxy, &b_mu_innerTrack_PV_dxy);
    fChain->SetBranchAddress("mu_innerTrack_PV_dz", &mu_innerTrack_PV_dz, &b_mu_innerTrack_PV_dz);
    //fChain->SetBranchAddress("mu_innerTrack_dxyError", &mu_innerTrack_dxyError, &b_mu_innerTrack_dxyError);
    //fChain->SetBranchAddress("mu_innerTrack_dzError", &mu_innerTrack_dzError, &b_mu_innerTrack_dzError);
    //fChain->SetBranchAddress("mu_innerTrack_validFraction", &mu_innerTrack_validFraction, &b_mu_innerTrack_validFraction);
    //fChain->SetBranchAddress("mu_bestTrack_dxy", &mu_bestTrack_dxy, &b_mu_bestTrack_dxy);
    //fChain->SetBranchAddress("mu_bestTrack_dz", &mu_bestTrack_dz, &b_mu_bestTrack_dz);
    //fChain->SetBranchAddress("mu_bestTrack_dxyError", &mu_bestTrack_dxyError, &b_mu_bestTrack_dxyError);
    //fChain->SetBranchAddress("mu_bestTrack_dzError", &mu_bestTrack_dzError, &b_mu_bestTrack_dzError);
    fChain->SetBranchAddress("mu_bestTrack_pt",      &mu_bestTrack_pt,      &b_mu_bestTrack_pt     );
    fChain->SetBranchAddress("mu_bestTrack_ptError", &mu_bestTrack_ptError, &b_mu_bestTrack_ptError);
    fChain->SetBranchAddress("mu_numberOfMatches", &mu_numberOfMatches, &b_mu_numberOfMatches);
    //fChain->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits, &b_mu_numberOfValidMuonHits);
    fChain->SetBranchAddress("mu_pfIso03_sumChargedHadronPt", &mu_pfIso03_sumChargedHadronPt, &b_mu_pfIso03_sumChargedHadronPt);
    fChain->SetBranchAddress("mu_pfIso03_sumNeutralHadronEt", &mu_pfIso03_sumNeutralHadronEt, &b_mu_pfIso03_sumNeutralHadronEt);
    fChain->SetBranchAddress("mu_pfIso03_sumPhotonEt", &mu_pfIso03_sumPhotonEt, &b_mu_pfIso03_sumPhotonEt);
    fChain->SetBranchAddress("mu_pfIso03_sumPUPt", &mu_pfIso03_sumPUPt, &b_mu_pfIso03_sumPUPt);
    fChain->SetBranchAddress("mu_lepMVA", &mu_lepMVA, &b_mu_lepMVA);
    fChain->SetBranchAddress("mu_lepMVA_miniRelIsoCharged", &mu_lepMVA_miniRelIsoCharged, &b_mu_lepMVA_miniRelIsoCharged);
    fChain->SetBranchAddress("mu_lepMVA_miniRelIsoNeutral", &mu_lepMVA_miniRelIsoNeutral, &b_mu_lepMVA_miniRelIsoNeutral);
    fChain->SetBranchAddress("mu_lepMVA_jetPtRelv2", &mu_lepMVA_jetPtRelv2, &b_mu_lepMVA_jetPtRelv2);
    //fChain->SetBranchAddress("mu_lepMVA_neuRelIso", &mu_lepMVA_neuRelIso, &b_mu_lepMVA_neuRelIso);
    //fChain->SetBranchAddress("mu_lepMVA_chRelIso", &mu_lepMVA_chRelIso, &b_mu_lepMVA_chRelIso);
    fChain->SetBranchAddress("mu_lepMVA_jetDR", &mu_lepMVA_jetDR, &b_mu_lepMVA_jetDR);
    fChain->SetBranchAddress("mu_lepMVA_jetPtRatio", &mu_lepMVA_jetPtRatio, &b_mu_lepMVA_jetPtRatio);
    fChain->SetBranchAddress("mu_lepMVA_jetBTagCSV", &mu_lepMVA_jetBTagCSV, &b_mu_lepMVA_jetBTagCSV);
    fChain->SetBranchAddress("mu_lepMVA_sip3d", &mu_lepMVA_sip3d, &b_mu_lepMVA_sip3d);
    fChain->SetBranchAddress("mu_lepMVA_dxy", &mu_lepMVA_dxy, &b_mu_lepMVA_dxy);
    fChain->SetBranchAddress("mu_lepMVA_dz", &mu_lepMVA_dz, &b_mu_lepMVA_dz);
    fChain->SetBranchAddress("mu_lepMVA_mvaId", &mu_lepMVA_mvaId, &b_mu_lepMVA_mvaId);
    fChain->SetBranchAddress("mu_innerTrack_pt", &mu_innerTrack_pt, &b_mu_innerTrack_pt);
    fChain->SetBranchAddress("mu_innerTrack_ptError", &mu_innerTrack_ptError, &b_mu_innerTrack_ptError);
    fChain->SetBranchAddress("mu_dB3D", &mu_dB3D, &b_mu_dB3D);
    fChain->SetBranchAddress("mu_edB3D", &mu_edB3D, &b_mu_edB3D);
    fChain->SetBranchAddress("mu_lepMVA_Moriond16", &mu_lepMVA_Moriond16, &b_mu_lepMVA_Moriond16);
    fChain->SetBranchAddress("mu_lepMVA_eta", &mu_lepMVA_eta, &b_mu_lepMVA_eta);
    fChain->SetBranchAddress("mu_lepMVA_jetNDauChargedMVASel", &mu_lepMVA_jetNDauChargedMVASel, &b_mu_lepMVA_jetNDauChargedMVASel);   

    // #########################
    // #  _____                #
    // # |_   _|_ _ _   _ ___  #
    // #   | |/ _` | | | / __| #
    // #   | | (_| | |_| \__ \ #
    // #   |_|\__,_|\__,_|___/ #
    // #                       #
    // #########################

    fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
    fChain->SetBranchAddress("tau_E", &tau_E,  &b_tau_E);
    fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
    fChain->SetBranchAddress("tau_eta", &tau_eta,  &b_tau_eta);
    fChain->SetBranchAddress("tau_phi", &tau_phi,  &b_tau_phi);
    fChain->SetBranchAddress("tau_m", &tau_m,  &b_tau_m);
    fChain->SetBranchAddress("tau_dxy", &tau_dxy,  &b_tau_dxy);
    fChain->SetBranchAddress("tau_dz", &tau_dz,  &b_tau_dz);
    fChain->SetBranchAddress("tau_leadingTrackDxy", &tau_leadingTrackDxy, &b_tau_leadingTrackDxy);
    fChain->SetBranchAddress("tau_leadingTrackDz", & tau_leadingTrackDz, &b_tau_leadingTrackDz);
    fChain->SetBranchAddress("tau_charge", &tau_charge,  &b_tau_charge);
    fChain->SetBranchAddress("tau_id", &tau_id,  &b_tau_id);
    fChain->SetBranchAddress("tau_decayMode", &tau_decayMode, &b_tau_decayMode);
    fChain->SetBranchAddress("tau_hasLeadChargedHadrCand", &tau_hasLeadChargedHadrCand,  &b_tau_hasLeadChargedHadrCand);
    fChain->SetBranchAddress("tau_leadingTrackPt", &tau_leadingTrackPt,  &b_tau_leadingTrackPt);
    fChain->SetBranchAddress("tau_leadingTrackCharge", &tau_leadingTrackCharge,  &b_tau_leadingTrackCharge);
    fChain->SetBranchAddress("tau_decayModeFindingOldDMs", &tau_decayModeFindingOldDMs, &b_tau_decayModeFindingOldDMs);
    fChain->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
    fChain->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,&b_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    fChain->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
    fChain->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
    fChain->SetBranchAddress("tau_byLooseIsolationMVA3newDMwLT", &tau_byLooseIsolationMVA3newDMwLT, &b_tau_byLooseIsolationMVA3newDMwLT);
    fChain->SetBranchAddress("tau_byMediumIsolationMVA3newDMwLT", &tau_byMediumIsolationMVA3newDMwLT, &b_tau_byMediumIsolationMVA3newDMwLT);
    fChain->SetBranchAddress("tau_byTightIsolationMVA3newDMwLT", &tau_byTightIsolationMVA3newDMwLT, &b_tau_byTightIsolationMVA3newDMwLT);
    fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);

    fChain->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",  &tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03,  &b_tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03  );
    fChain->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03, &b_tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 );
    fChain->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",  &tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03,  &b_tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03  );
    fChain->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT",         &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT,         &b_tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT         );
    fChain->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT",        &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT,        &b_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT        );
    fChain->SetBranchAddress("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT",         &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT,         &b_tau_byTightIsolationMVArun2v1DBdR03oldDMwLT         ); 
    fChain->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT",        &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT,        &b_tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT        );
       
    fChain->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
    fChain->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
    fChain->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
    fChain->SetBranchAddress("tau_againstMuonLoose3", &tau_againstMuonLoose3, &b_tau_againstMuonLoose3);
    fChain->SetBranchAddress("tau_againstMuonTight3", &tau_againstMuonTight3, &b_tau_againstMuonTight3);
    //AC8X
    fChain->SetBranchAddress("tau_againstElectronVLooseMVA6", &tau_againstElectronVLooseMVA6, &b_tau_againstElectronVLooseMVA6);
    fChain->SetBranchAddress("tau_againstElectronLooseMVA6", &tau_againstElectronLooseMVA6, &b_tau_againstElectronLooseMVA6);
    fChain->SetBranchAddress("tau_againstElectronMediumMVA6", &tau_againstElectronMediumMVA6, &b_tau_againstElectronMediumMVA6);
    fChain->SetBranchAddress("tau_againstElectronTightMVA6", &tau_againstElectronTightMVA6, &b_tau_againstElectronTightMVA6);
    fChain->SetBranchAddress("tau_pfEssential_jet_pt", &tau_pfEssential_jet_pt, &b_tau_pfEssential_jet_pt);
    fChain->SetBranchAddress("tau_pfEssential_jet_eta", &tau_pfEssential_jet_eta, &b_tau_pfEssential_jet_eta);
    fChain->SetBranchAddress("tau_pfEssential_jet_phi", &tau_pfEssential_jet_phi, &b_tau_pfEssential_jet_phi);
    fChain->SetBranchAddress("tau_pfEssential_jet_m", &tau_pfEssential_jet_m, &b_tau_pfEssential_jet_m);
    fChain->SetBranchAddress("tau_pfEssential_jetCorr_pt", &tau_pfEssential_jetCorr_pt, &b_tau_pfEssential_jetCorr_pt);
    fChain->SetBranchAddress("tau_pfEssential_jetCorr_eta", &tau_pfEssential_jetCorr_eta, &b_tau_pfEssential_jetCorr_eta);
    fChain->SetBranchAddress("tau_pfEssential_jetCorr_phi", &tau_pfEssential_jetCorr_phi, &b_tau_pfEssential_jetCorr_phi);
    fChain->SetBranchAddress("tau_pfEssential_jetCorr_m", &tau_pfEssential_jetCorr_m, &b_tau_pfEssential_jetCorr_m);
    fChain->SetBranchAddress("tau_pfEssential_hasSV", &tau_pfEssential_hasSV, &b_tau_pfEssential_hasSV);
    fChain->SetBranchAddress("tau_pfEssential_sv_x", &tau_pfEssential_sv_x, &b_tau_pfEssential_sv_x);
    fChain->SetBranchAddress("tau_pfEssential_sv_y", &tau_pfEssential_sv_y, &b_tau_pfEssential_sv_y);
    fChain->SetBranchAddress("tau_pfEssential_sv_z", &tau_pfEssential_sv_z, &b_tau_pfEssential_sv_z);
    fChain->SetBranchAddress("tau_pfEssential_flightLengthSig", &tau_pfEssential_flightLengthSig, &b_tau_pfEssential_flightLengthSig);
    fChain->SetBranchAddress("tau_pfEssential_dxy", &tau_pfEssential_dxy, &b_tau_pfEssential_dxy);
    fChain->SetBranchAddress("tau_pfEssential_dxy_error", &tau_pfEssential_dxy_error, &b_tau_pfEssential_dxy_error);
    fChain->SetBranchAddress("tau_pfEssential_dxy_Sig", &tau_pfEssential_dxy_Sig, &b_tau_pfEssential_dxy_Sig);

    // ##########################
    // #       _      _         #
    // #      | | ___| |_ ___   #
    // #   _  | |/ _ \ __/ __|  #
    // #  | |_| |  __/ |_\__ \  #
    // #   \___/ \___|\__|___/  #
    // #                        #                        
    // ##########################

    fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
    fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
    fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
    fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
    fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
    fChain->SetBranchAddress("jet_E", &jet_E, &b_jet_E);
    fChain->SetBranchAddress("jet_qgtag", &jet_qgtag, &b_jet_qgtag);
    fChain->SetBranchAddress("jet_ntrk", &jet_ntrk, &b_jet_ntrk);
    fChain->SetBranchAddress("jet_CSVv2", &jet_CSVv2, &b_jet_CSVv2);
    fChain->SetBranchAddress("jet_cMVAv2", &jet_cMVAv2, &b_jet_cMVAv2);
    fChain->SetBranchAddress("jet_DeepCSVProbudsg", &jet_DeepCSVProbudsg, &b_jet_DeepCSVProbudsg);
    fChain->SetBranchAddress("jet_DeepCSVProbb", &jet_DeepCSVProbb, &b_jet_DeepCSVProbb);
    fChain->SetBranchAddress("jet_DeepCSVProbbb", &jet_DeepCSVProbb, &b_jet_DeepCSVProbbb);
    fChain->SetBranchAddress("jet_DeepCSVProbc", &jet_DeepCSVProbc, &b_jet_DeepCSVProbc);
    fChain->SetBranchAddress("jet_DeepCSVProbcc", &jet_DeepCSVProbcc, &b_jet_DeepCSVProbcc);
    fChain->SetBranchAddress("jet_looseJetID", &jet_looseJetID, &b_jet_looseJetID);
    fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
    fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
    fChain->SetBranchAddress("jet_neutralHadronEnergy", &jet_neutralHadronEnergy, &b_jet_neutralHadronEnergy);
    fChain->SetBranchAddress("jet_neutralEmEnergy", &jet_neutralEmEnergy, &b_jet_neutralEmEnergy);
    fChain->SetBranchAddress("jet_chargedHadronEnergy", &jet_chargedHadronEnergy, &b_jet_chargedHadronEnergy);
    fChain->SetBranchAddress("jet_chargedEmEnergy", &jet_chargedEmEnergy, &b_jet_chargedEmEnergy);
    fChain->SetBranchAddress("jet_electronEnergy", &jet_electronEnergy, &b_jet_electronEnergy);
    fChain->SetBranchAddress("jet_muonEnergy", &jet_muonEnergy, &b_jet_muonEnergy);
    fChain->SetBranchAddress("jet_photonEnergy", &jet_photonEnergy, &b_jet_photonEnergy);
    fChain->SetBranchAddress("jet_genJet_pt", &jet_genJet_pt, &b_jet_genJet_pt);
    fChain->SetBranchAddress("jet_genParton_pt", &jet_genParton_pt, &b_jet_genParton_pt);
    fChain->SetBranchAddress("jet_genParton_id", &jet_genParton_id, &b_jet_genParton_id);
    fChain->SetBranchAddress("jet_genParton_E", &jet_genParton_E, &b_jet_genParton_E);
    //fChain->SetBranchAddress("jet_genJet_eta", &jet_genJet_eta, &b_jet_genJet_eta);
    //fChain->SetBranchAddress("jet_genJet_phi", &jet_genJet_phi, &b_jet_genJet_phi);
    //fChain->SetBranchAddress("jet_genJet_m", &jet_genJet_m, &b_jet_genJet_m);
    fChain->SetBranchAddress("jet_genJet_E", &jet_genJet_E, &b_jet_genJet_E);
    //fChain->SetBranchAddress("jet_genJet_status", &jet_genJet_status, &b_jet_genJet_status);
    //fChain->SetBranchAddress("jet_genJet_id", &jet_genJet_id, &b_jet_genJet_id);
    fChain->SetBranchAddress("jet_pileupJetId", &jet_pileupJetId, &b_jet_pileupJetId);

    if( fChain->GetBranch("gen_n") ) fChain->SetBranchAddress("gen_n", &gen_n, &b_gen_n);
    if( fChain->GetBranch("gen_pt") ) fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
    if( fChain->GetBranch("gen_eta") ) fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
    if( fChain->GetBranch("gen_phi") ) fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
    if( fChain->GetBranch("gen_m") ) fChain->SetBranchAddress("gen_m", &gen_m, &b_gen_m);
    if( fChain->GetBranch("gen_id") ) fChain->SetBranchAddress("gen_id", &gen_id, &b_gen_id);
    if( fChain->GetBranch("gen_status") ) fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);
    if( fChain->GetBranch("gen_index") ) fChain->SetBranchAddress("gen_index", &gen_index, &b_gen_index);
    if( fChain->GetBranch("gen_mother_index") ) fChain->SetBranchAddress("gen_mother_index", &gen_mother_index, &b_gen_mother_index);
    if( fChain->GetBranch("gen_PVz") ) fChain->SetBranchAddress("gen_PVz", &gen_PVz, &b_gen_PVz);

    // #########################################
    // #                                       #
    // #   __ _  ___ _ __     (_) ___| |_ ___  #
    // #  / _` |/ _ \ '_ \    | |/ _ \ __/ __| #
    // # | (_| |  __/ | | |   | |  __/ |_\__ \ #
    // #  \__, |\___|_| |_|  _/ |\___|\__|___/ #
    // #  |___/             |__/               #
    // #                                       #            
    // #########################################

    if( fChain->GetBranch("genJet_n") ) fChain->SetBranchAddress("genJet_n", &genJet_n, &b_genJet_n);
    if( fChain->GetBranch("genJet_pt") ) fChain->SetBranchAddress("genJet_pt", &genJet_pt, &b_genJet_pt);
    if( fChain->GetBranch("genJet_eta") ) fChain->SetBranchAddress("genJet_eta", &genJet_eta, &b_genJet_eta);
    if( fChain->GetBranch("genJet_phi") ) fChain->SetBranchAddress("genJet_phi", &genJet_phi, &b_genJet_phi);
    if( fChain->GetBranch("genJet_E") ) fChain->SetBranchAddress("genJet_E", &genJet_E, &b_genJet_E);
    if( fChain->GetBranch("genJet_m") ) fChain->SetBranchAddress("genJet_m", &genJet_m, &b_genJet_m);

    if( fChain->GetBranch("metGen_px") ) fChain->SetBranchAddress("metGen_px", &metGen_px, &b_metGen_px);
    if( fChain->GetBranch("metGen_py") ) fChain->SetBranchAddress("metGen_py", &metGen_py, &b_metGen_py);
    if( fChain->GetBranch("metGen_pt") ) fChain->SetBranchAddress("metGen_pt", &metGen_pt, &b_metGen_pt);
    if( fChain->GetBranch("metGen_phi") ) fChain->SetBranchAddress("metGen_phi", &metGen_phi, &b_metGen_phi);
    if( fChain->GetBranch("metGen_sumet") ) fChain->SetBranchAddress("metGen_sumet", &metGen_sumet, &b_metGen_sumet);
    if( fChain->GetBranch("metGen_MuonEt") ) fChain->SetBranchAddress("metGen_MuonEt", &metGen_MuonEt, &b_metGen_MuonEt);

    if( fChain->GetBranch("triggerobject_n") ) fChain->SetBranchAddress("triggerobject_n", &triggerobject_n, &b_triggerobject_n);
    if( fChain->GetBranch("triggerobject_pt") ) fChain->SetBranchAddress("triggerobject_pt", &triggerobject_pt, &b_triggerobject_pt);
    if( fChain->GetBranch("triggerobject_eta") ) fChain->SetBranchAddress("triggerobject_eta", &triggerobject_eta, &b_triggerobject_eta);
    if( fChain->GetBranch("triggerobject_phi") ) fChain->SetBranchAddress("triggerobject_phi", &triggerobject_phi, &b_triggerobject_phi);
    if( fChain->GetBranch("triggerobject_collection") ) fChain->SetBranchAddress("triggerobject_collection", &triggerobject_collection, &b_triggerobject_collection);

    if( fChain->GetBranch("triggerobject_filterIds_n") ) fChain->SetBranchAddress("triggerobject_filterIds_n", &triggerobject_filterIds_n, &b_triggerobject_filterIds_n);
    if( fChain->GetBranch("triggerobject_filterLabels_n") ) fChain->SetBranchAddress("triggerobject_filterLabels_n", &triggerobject_filterLabels_n, &b_triggerobject_filterLabels_n);
    if( fChain->GetBranch("triggerobject_pathNamesAll_n") ) fChain->SetBranchAddress("triggerobject_pathNamesAll_n", &triggerobject_pathNamesAll_n, &b_triggerobject_pathNamesAll_n);

    if( fChain->GetBranch("triggerobject_pathNamesAll") ) fChain->SetBranchAddress("triggerobject_pathNamesAll", &triggerobject_pathNamesAll, &b_triggerobject_pathNamesAll);
    if( fChain->GetBranch("triggerobject_pathNamesAll_isL3") ) fChain->SetBranchAddress("triggerobject_pathNamesAll_isL3", &triggerobject_pathNamesAll_isL3, &b_triggerobject_pathNamesAll_isL3);
    if( fChain->GetBranch("triggerobject_pathNamesAll_isLF") ) fChain->SetBranchAddress("triggerobject_pathNamesAll_isLF", &triggerobject_pathNamesAll_isLF, &b_triggerobject_pathNamesAll_isLF);
    if( fChain->GetBranch("triggerobject_pathNamesAll_isBoth") ) fChain->SetBranchAddress("triggerobject_pathNamesAll_isBoth", &triggerobject_pathNamesAll_isBoth, &b_triggerobject_pathNamesAll_isBoth);
    if( fChain->GetBranch("triggerobject_pathNamesAll_isNone") ) fChain->SetBranchAddress("triggerobject_pathNamesAll_isNone", &triggerobject_pathNamesAll_isNone, &b_triggerobject_pathNamesAll_isNone);

    if( fChain->GetBranch("triggerobject_filterIds") ) fChain->SetBranchAddress("triggerobject_filterIds", &triggerobject_filterIds, &b_triggerobject_filterIds);
    if( fChain->GetBranch("triggerobject_filterLabels") ) fChain->SetBranchAddress("triggerobject_filterLabels", &triggerobject_filterLabels, &b_triggerobject_filterLabels);



    // ##################################################
    // #   __  __  ____     _____           _   _       #
    // #  |  \/  |/ ___|   |_   _| __ _   _| |_| |__    #
    // #  | |\/| | |         | || '__| | | | __| '_ \   #
    // #  | |  | | |___      | || |  | |_| | |_| | | |  #
    // #  |_|  |_|\____|     |_||_|   \__,_|\__|_| |_|  #
    // #                                                #
    // ##################################################

    if( fChain->GetBranch("mc_truth_h0_id") ) fChain->SetBranchAddress("mc_truth_h0_id", &mc_truth_h0_id, &b_mc_truth_h0_id);
    if( fChain->GetBranch("mc_truth_h0W1_id") ) fChain->SetBranchAddress("mc_truth_h0W1_id", &mc_truth_h0W1_id, &b_mc_truth_h0W1_id);
    if( fChain->GetBranch("mc_truth_h0Wl1_id") ) fChain->SetBranchAddress("mc_truth_h0Wl1_id", &mc_truth_h0Wl1_id, &b_mc_truth_h0Wl1_id);
    if( fChain->GetBranch("mc_truth_h0Wtau1_id") ) fChain->SetBranchAddress("mc_truth_h0Wtau1_id", &mc_truth_h0Wtau1_id, &b_mc_truth_h0Wtau1_id);    
    if( fChain->GetBranch("mc_truth_h0Wtaul1_id") ) fChain->SetBranchAddress("mc_truth_h0Wtaul1_id", &mc_truth_h0Wtaul1_id, &b_mc_truth_h0Wtaul1_id);    
    if( fChain->GetBranch("mc_truth_h0Wq11_id") ) fChain->SetBranchAddress("mc_truth_h0Wq11_id", &mc_truth_h0Wq11_id, &b_mc_truth_h0Wq11_id);
    if( fChain->GetBranch("mc_truth_h0Wq21_id") ) fChain->SetBranchAddress("mc_truth_h0Wq21_id", &mc_truth_h0Wq21_id, &b_mc_truth_h0Wq21_id);
    if( fChain->GetBranch("mc_truth_h0W2_id") ) fChain->SetBranchAddress("mc_truth_h0W2_id", &mc_truth_h0W2_id, &b_mc_truth_h0W2_id);
    if( fChain->GetBranch("mc_truth_h0Wl2_id") ) fChain->SetBranchAddress("mc_truth_h0Wl2_id", &mc_truth_h0Wl2_id, &b_mc_truth_h0Wl2_id);
    if( fChain->GetBranch("mc_truth_h0Wtau2_id") ) fChain->SetBranchAddress("mc_truth_h0Wtau2_id", &mc_truth_h0Wtau2_id, &b_mc_truth_h0Wtau2_id);    
    if( fChain->GetBranch("mc_truth_h0Wtaul2_id") ) fChain->SetBranchAddress("mc_truth_h0Wtaul2_id", &mc_truth_h0Wtaul2_id, &b_mc_truth_h0Wtaul2_id);      
    if( fChain->GetBranch("mc_truth_h0Wq12_id") ) fChain->SetBranchAddress("mc_truth_h0Wq12_id", &mc_truth_h0Wq12_id, &b_mc_truth_h0Wq12_id);
    if( fChain->GetBranch("mc_truth_h0Wq22_id") ) fChain->SetBranchAddress("mc_truth_h0Wq22_id", &mc_truth_h0Wq22_id, &b_mc_truth_h0Wq22_id);
    if( fChain->GetBranch("mc_truth_h0Z1_id") ) fChain->SetBranchAddress("mc_truth_h0Z1_id", &mc_truth_h0Z1_id, &b_mc_truth_h0Z1_id);
    if( fChain->GetBranch("mc_truth_h0Zl11_id") ) fChain->SetBranchAddress("mc_truth_h0Zl11_id", &mc_truth_h0Zl11_id, &b_mc_truth_h0Zl11_id);
    if( fChain->GetBranch("mc_truth_h0Zl21_id") ) fChain->SetBranchAddress("mc_truth_h0Zl21_id", &mc_truth_h0Zl21_id, &b_mc_truth_h0Zl21_id);
    if( fChain->GetBranch("mc_truth_h0Zq11_id") ) fChain->SetBranchAddress("mc_truth_h0Zq11_id", &mc_truth_h0Zq11_id, &b_mc_truth_h0Zq11_id);
    if( fChain->GetBranch("mc_truth_h0Zq21_id") ) fChain->SetBranchAddress("mc_truth_h0Zq21_id", &mc_truth_h0Zq21_id, &b_mc_truth_h0Zq21_id);
    if( fChain->GetBranch("mc_truth_h0Z2_id") ) fChain->SetBranchAddress("mc_truth_h0Z2_id", &mc_truth_h0Z2_id, &b_mc_truth_h0Z2_id);
    if( fChain->GetBranch("mc_truth_h0Zl12_id") ) fChain->SetBranchAddress("mc_truth_h0Zl12_id", &mc_truth_h0Zl12_id, &b_mc_truth_h0Zl12_id);
    if( fChain->GetBranch("mc_truth_h0Zl22_id") ) fChain->SetBranchAddress("mc_truth_h0Zl22_id", &mc_truth_h0Zl22_id, &b_mc_truth_h0Zl22_id);
    if( fChain->GetBranch("mc_truth_h0Zq12_id") ) fChain->SetBranchAddress("mc_truth_h0Zq12_id", &mc_truth_h0Zq12_id, &b_mc_truth_h0Zq12_id);
    if( fChain->GetBranch("mc_truth_h0Zq22_id") ) fChain->SetBranchAddress("mc_truth_h0Zq22_id", &mc_truth_h0Zq22_id, &b_mc_truth_h0Zq22_id);
    if( fChain->GetBranch("mc_truth_h0tau1_id") ) fChain->SetBranchAddress("mc_truth_h0tau1_id", &mc_truth_h0tau1_id, &b_mc_truth_h0tau1_id);
    if( fChain->GetBranch("mc_truth_h0tau2_id") ) fChain->SetBranchAddress("mc_truth_h0tau2_id", &mc_truth_h0tau2_id, &b_mc_truth_h0tau2_id);
    if( fChain->GetBranch("mc_truth_h0taul1_id") ) fChain->SetBranchAddress("mc_truth_h0taul1_id", &mc_truth_h0taul1_id, &b_mc_truth_h0taul1_id);
    if( fChain->GetBranch("mc_truth_h0taul2_id") ) fChain->SetBranchAddress("mc_truth_h0taul2_id", &mc_truth_h0taul2_id, &b_mc_truth_h0taul2_id);
    if( fChain->GetBranch("mc_truth_h0b1_id") ) fChain->SetBranchAddress("mc_truth_h0b1_id", &mc_truth_h0b1_id, &b_mc_truth_h0b1_id);
    if( fChain->GetBranch("mc_truth_h0b2_id") ) fChain->SetBranchAddress("mc_truth_h0b2_id", &mc_truth_h0b2_id, &b_mc_truth_h0b2_id);

    if( fChain->GetBranch("mc_truth_t1_id") ) fChain->SetBranchAddress("mc_truth_t1_id", &mc_truth_t1_id, &b_mc_truth_t1_id);
    if( fChain->GetBranch("mc_truth_t2_id") ) fChain->SetBranchAddress("mc_truth_t2_id", &mc_truth_t2_id, &b_mc_truth_t2_id);
    if( fChain->GetBranch("mc_truth_tb1_id") ) fChain->SetBranchAddress("mc_truth_tb1_id", &mc_truth_tb1_id, &b_mc_truth_tb1_id);
    if( fChain->GetBranch("mc_truth_tb2_id") ) fChain->SetBranchAddress("mc_truth_tb2_id", &mc_truth_tb2_id, &b_mc_truth_tb2_id);
    if( fChain->GetBranch("mc_truth_tW1_id") ) fChain->SetBranchAddress("mc_truth_tW1_id", &mc_truth_tW1_id, &b_mc_truth_tW1_id);
    if( fChain->GetBranch("mc_truth_tW2_id") ) fChain->SetBranchAddress("mc_truth_tW2_id", &mc_truth_tW2_id, &b_mc_truth_tW2_id);
    if( fChain->GetBranch("mc_truth_tWl1_id") ) fChain->SetBranchAddress("mc_truth_tWl1_id", &mc_truth_tWl1_id, &b_mc_truth_tWl1_id);
    if( fChain->GetBranch("mc_truth_tWl2_id") ) fChain->SetBranchAddress("mc_truth_tWl2_id", &mc_truth_tWl2_id, &b_mc_truth_tWl2_id);
    if( fChain->GetBranch("mc_truth_tWtau1_id") ) fChain->SetBranchAddress("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, &b_mc_truth_tWtau1_id);
    if( fChain->GetBranch("mc_truth_tWtau2_id") ) fChain->SetBranchAddress("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, &b_mc_truth_tWtau2_id);
    if( fChain->GetBranch("mc_truth_tWtaul1_id") ) fChain->SetBranchAddress("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, &b_mc_truth_tWtaul1_id);
    if( fChain->GetBranch("mc_truth_tWtaul2_id") ) fChain->SetBranchAddress("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, &b_mc_truth_tWtaul2_id);
    if( fChain->GetBranch("mc_truth_tWq11_id") ) fChain->SetBranchAddress("mc_truth_tWq11_id", &mc_truth_tWq11_id, &b_mc_truth_tWq11_id);
    if( fChain->GetBranch("mc_truth_tWq21_id") ) fChain->SetBranchAddress("mc_truth_tWq21_id", &mc_truth_tWq21_id, &b_mc_truth_tWq21_id);
    if( fChain->GetBranch("mc_truth_tWq12_id") ) fChain->SetBranchAddress("mc_truth_tWq12_id", &mc_truth_tWq12_id, &b_mc_truth_tWq12_id);
    if( fChain->GetBranch("mc_truth_tWq22_id") ) fChain->SetBranchAddress("mc_truth_tWq22_id", &mc_truth_tWq22_id, &b_mc_truth_tWq22_id);   

    if( fChain->GetBranch("mc_truth_t_id") ) fChain->SetBranchAddress("mc_truth_t_id", &mc_truth_t_id, &b_mc_truth_t_id);
    if( fChain->GetBranch("mc_truth_tb_id") ) fChain->SetBranchAddress("mc_truth_tb_id", &mc_truth_tb_id, &b_mc_truth_tb_id);
    if( fChain->GetBranch("mc_truth_tW_id") ) fChain->SetBranchAddress("mc_truth_tW_id", &mc_truth_tW_id, &b_mc_truth_tW_id);
    if( fChain->GetBranch("mc_truth_tWl_id") ) fChain->SetBranchAddress("mc_truth_tWl_id", &mc_truth_tWl_id, &b_mc_truth_tWl_id);
    if( fChain->GetBranch("mc_truth_tWq1_id") ) fChain->SetBranchAddress("mc_truth_tWq1_id", &mc_truth_tWq1_id, &b_mc_truth_tWq1_id);
    if( fChain->GetBranch("mc_truth_tWq2_id") ) fChain->SetBranchAddress("mc_truth_tWq2_id", &mc_truth_tWq2_id, &b_mc_truth_tWq2_id);

    if( fChain->GetBranch("mc_truth_W_id") ) fChain->SetBranchAddress("mc_truth_W_id", &mc_truth_W_id, &b_mc_truth_W_id);
    if( fChain->GetBranch("mc_truth_Wl_id") ) fChain->SetBranchAddress("mc_truth_Wl_id", &mc_truth_Wl_id, &b_mc_truth_Wl_id);
    if( fChain->GetBranch("mc_truth_Wtau_id") ) fChain->SetBranchAddress("mc_truth_Wtau_id", &mc_truth_Wtau_id, &b_mc_truth_Wtau_id);
    if( fChain->GetBranch("mc_truth_Wtaul_id") ) fChain->SetBranchAddress("mc_truth_Wtaul_id", &mc_truth_Wtaul_id, &b_mc_truth_Wtaul_id);
    if( fChain->GetBranch("mc_truth_Wq1_id") ) fChain->SetBranchAddress("mc_truth_Wq1_id", &mc_truth_Wq1_id, &b_mc_truth_Wq1_id);
    if( fChain->GetBranch("mc_truth_Wq2_id") ) fChain->SetBranchAddress("mc_truth_Wq2_id", &mc_truth_Wq2_id, &b_mc_truth_Wq2_id);
    if( fChain->GetBranch("mc_truth_Z_id") ) fChain->SetBranchAddress("mc_truth_Z_id", &mc_truth_Z_id, &b_mc_truth_Z_id);
    if( fChain->GetBranch("mc_truth_Zl1_id") ) fChain->SetBranchAddress("mc_truth_Zl1_id", &mc_truth_Zl1_id, &b_mc_truth_Zl1_id);
    if( fChain->GetBranch("mc_truth_Zl2_id") ) fChain->SetBranchAddress("mc_truth_Zl2_id", &mc_truth_Zl2_id, &b_mc_truth_Zl2_id);
    if( fChain->GetBranch("mc_truth_Ztau1_id") ) fChain->SetBranchAddress("mc_truth_Ztau1_id", &mc_truth_Ztau1_id, &b_mc_truth_Ztau1_id);
    if( fChain->GetBranch("mc_truth_Ztau2_id") ) fChain->SetBranchAddress("mc_truth_Ztau2_id", &mc_truth_Ztau2_id, &b_mc_truth_Ztau2_id);
    if( fChain->GetBranch("mc_truth_Ztaul1_id") ) fChain->SetBranchAddress("mc_truth_Ztaul1_id", &mc_truth_Ztaul1_id, &b_mc_truth_Ztaul1_id);
    if( fChain->GetBranch("mc_truth_Ztaul2_id") ) fChain->SetBranchAddress("mc_truth_Ztaul2_id", &mc_truth_Ztaul2_id, &b_mc_truth_Ztaul2_id);    
    if( fChain->GetBranch("mc_truth_Zq1_id") ) fChain->SetBranchAddress("mc_truth_Zq1_id", &mc_truth_Zq1_id, &b_mc_truth_Zq1_id);
    if( fChain->GetBranch("mc_truth_Zq2_id") ) fChain->SetBranchAddress("mc_truth_Zq2_id", &mc_truth_Zq2_id, &b_mc_truth_Zq2_id);
    if( fChain->GetBranch("mc_truth_gammal1_id") ) fChain->SetBranchAddress("mc_truth_gammal1_id", &mc_truth_gammal1_id, &b_mc_truth_gammal1_id);
    if( fChain->GetBranch("mc_truth_gammal2_id") ) fChain->SetBranchAddress("mc_truth_gammal2_id", &mc_truth_gammal2_id, &b_mc_truth_gammal2_id);
    if( fChain->GetBranch("mc_truth_gammatau1_id") ) fChain->SetBranchAddress("mc_truth_gammatau1_id", &mc_truth_gammatau1_id, &b_mc_truth_gammatau1_id);
    if( fChain->GetBranch("mc_truth_gammatau2_id") ) fChain->SetBranchAddress("mc_truth_gammatau2_id", &mc_truth_gammatau2_id, &b_mc_truth_gammatau2_id);
    if( fChain->GetBranch("mc_truth_gammataul1_id") ) fChain->SetBranchAddress("mc_truth_gammataul1_id", &mc_truth_gammataul1_id, &b_mc_truth_gammataul1_id);
    if( fChain->GetBranch("mc_truth_gammataul2_id") ) fChain->SetBranchAddress("mc_truth_gammataul2_id", &mc_truth_gammataul2_id, &b_mc_truth_gammataul2_id);
    if( fChain->GetBranch("mc_truth_gamma_id") ) fChain->SetBranchAddress("mc_truth_gamma_id", &mc_truth_gamma_id, &b_mc_truth_gamma_id);

    if( fChain->GetBranch("mc_truth_h0_pt") ) fChain->SetBranchAddress("mc_truth_h0_pt", &mc_truth_h0_pt, &b_mc_truth_h0_pt);
    if( fChain->GetBranch("mc_truth_h0W1_pt") ) fChain->SetBranchAddress("mc_truth_h0W1_pt", &mc_truth_h0W1_pt, &b_mc_truth_h0W1_pt);
    if( fChain->GetBranch("mc_truth_h0Wl1_pt") ) fChain->SetBranchAddress("mc_truth_h0Wl1_pt", &mc_truth_h0Wl1_pt, &b_mc_truth_h0Wl1_pt);
    if( fChain->GetBranch("mc_truth_h0Wtau1_pt") ) fChain->SetBranchAddress("mc_truth_h0Wtau1_pt", &mc_truth_h0Wtau1_pt, &b_mc_truth_h0Wtau1_pt);
    if( fChain->GetBranch("mc_truth_h0Wtaul1_pt") ) fChain->SetBranchAddress("mc_truth_h0Wtaul1_pt", &mc_truth_h0Wtaul1_pt, &b_mc_truth_h0Wtaul1_pt);
    if( fChain->GetBranch("mc_truth_h0Wq11_pt") ) fChain->SetBranchAddress("mc_truth_h0Wq11_pt", &mc_truth_h0Wq11_pt, &b_mc_truth_h0Wq11_pt);
    if( fChain->GetBranch("mc_truth_h0Wq21_pt") ) fChain->SetBranchAddress("mc_truth_h0Wq21_pt", &mc_truth_h0Wq21_pt, &b_mc_truth_h0Wq21_pt);
    if( fChain->GetBranch("mc_truth_h0W2_pt") ) fChain->SetBranchAddress("mc_truth_h0W2_pt", &mc_truth_h0W2_pt, &b_mc_truth_h0W2_pt);
    if( fChain->GetBranch("mc_truth_h0Wl2_pt") ) fChain->SetBranchAddress("mc_truth_h0Wl2_pt", &mc_truth_h0Wl2_pt, &b_mc_truth_h0Wl2_pt);
    if( fChain->GetBranch("mc_truth_h0Wtau2_pt") ) fChain->SetBranchAddress("mc_truth_h0Wtau2_pt", &mc_truth_h0Wtau2_pt, &b_mc_truth_h0Wtau2_pt);
    if( fChain->GetBranch("mc_truth_h0Wtaul2_pt") ) fChain->SetBranchAddress("mc_truth_h0Wtaul2_pt", &mc_truth_h0Wtaul2_pt, &b_mc_truth_h0Wtaul2_pt);
    if( fChain->GetBranch("mc_truth_h0Wq12_pt") ) fChain->SetBranchAddress("mc_truth_h0Wq12_pt", &mc_truth_h0Wq12_pt, &b_mc_truth_h0Wq12_pt);
    if( fChain->GetBranch("mc_truth_h0Wq22_pt") ) fChain->SetBranchAddress("mc_truth_h0Wq22_pt", &mc_truth_h0Wq22_pt, &b_mc_truth_h0Wq22_pt);
    if( fChain->GetBranch("mc_truth_h0Z1_pt") ) fChain->SetBranchAddress("mc_truth_h0Z1_pt", &mc_truth_h0Z1_pt, &b_mc_truth_h0Z1_pt);
    if( fChain->GetBranch("mc_truth_h0Zl11_pt") ) fChain->SetBranchAddress("mc_truth_h0Zl11_pt", &mc_truth_h0Zl11_pt, &b_mc_truth_h0Zl11_pt);
    if( fChain->GetBranch("mc_truth_h0Zl21_pt") ) fChain->SetBranchAddress("mc_truth_h0Zl21_pt", &mc_truth_h0Zl21_pt, &b_mc_truth_h0Zl21_pt);
    if( fChain->GetBranch("mc_truth_h0Zq11_pt") ) fChain->SetBranchAddress("mc_truth_h0Zq11_pt", &mc_truth_h0Zq11_pt, &b_mc_truth_h0Zq11_pt);
    if( fChain->GetBranch("mc_truth_h0Zq21_pt") ) fChain->SetBranchAddress("mc_truth_h0Zq21_pt", &mc_truth_h0Zq21_pt, &b_mc_truth_h0Zq21_pt);
    if( fChain->GetBranch("mc_truth_h0Z2_pt") ) fChain->SetBranchAddress("mc_truth_h0Z2_pt", &mc_truth_h0Z2_pt, &b_mc_truth_h0Z2_pt);
    if( fChain->GetBranch("mc_truth_h0Zl12_pt") ) fChain->SetBranchAddress("mc_truth_h0Zl12_pt", &mc_truth_h0Zl12_pt, &b_mc_truth_h0Zl12_pt);
    if( fChain->GetBranch("mc_truth_h0Zl22_pt") ) fChain->SetBranchAddress("mc_truth_h0Zl22_pt", &mc_truth_h0Zl22_pt, &b_mc_truth_h0Zl22_pt);
    if( fChain->GetBranch("mc_truth_h0Zq12_pt") ) fChain->SetBranchAddress("mc_truth_h0Zq12_pt", &mc_truth_h0Zq12_pt, &b_mc_truth_h0Zq12_pt);
    if( fChain->GetBranch("mc_truth_h0Zq22_pt") ) fChain->SetBranchAddress("mc_truth_h0Zq22_pt", &mc_truth_h0Zq22_pt, &b_mc_truth_h0Zq22_pt);
    if( fChain->GetBranch("mc_truth_h0tau1_pt") ) fChain->SetBranchAddress("mc_truth_h0tau1_pt", &mc_truth_h0tau1_pt, &b_mc_truth_h0tau1_pt);
    if( fChain->GetBranch("mc_truth_h0tau2_pt") ) fChain->SetBranchAddress("mc_truth_h0tau2_pt", &mc_truth_h0tau2_pt, &b_mc_truth_h0tau2_pt);
    if( fChain->GetBranch("mc_truth_h0taul1_pt") ) fChain->SetBranchAddress("mc_truth_h0taul1_pt", &mc_truth_h0taul1_pt, &b_mc_truth_h0taul1_pt);
    if( fChain->GetBranch("mc_truth_h0taul2_pt") ) fChain->SetBranchAddress("mc_truth_h0taul2_pt", &mc_truth_h0taul2_pt, &b_mc_truth_h0taul2_pt);
    if( fChain->GetBranch("mc_truth_h0b1_pt") ) fChain->SetBranchAddress("mc_truth_h0b1_pt", &mc_truth_h0b1_pt, &b_mc_truth_h0b1_pt);
    if( fChain->GetBranch("mc_truth_h0b2_pt") ) fChain->SetBranchAddress("mc_truth_h0b2_pt", &mc_truth_h0b2_pt, &b_mc_truth_h0b2_pt);

    if( fChain->GetBranch("mc_truth_t1_pt") ) fChain->SetBranchAddress("mc_truth_t1_pt", &mc_truth_t1_pt, &b_mc_truth_t1_pt);
    if( fChain->GetBranch("mc_truth_t2_pt") ) fChain->SetBranchAddress("mc_truth_t2_pt", &mc_truth_t2_pt, &b_mc_truth_t2_pt);
    if( fChain->GetBranch("mc_truth_tb1_pt") ) fChain->SetBranchAddress("mc_truth_tb1_pt", &mc_truth_tb1_pt, &b_mc_truth_tb1_pt);
    if( fChain->GetBranch("mc_truth_tb2_pt") ) fChain->SetBranchAddress("mc_truth_tb2_pt", &mc_truth_tb2_pt, &b_mc_truth_tb2_pt);
    if( fChain->GetBranch("mc_truth_tW1_pt") ) fChain->SetBranchAddress("mc_truth_tW1_pt", &mc_truth_tW1_pt, &b_mc_truth_tW1_pt);
    if( fChain->GetBranch("mc_truth_tW2_pt") ) fChain->SetBranchAddress("mc_truth_tW2_pt", &mc_truth_tW2_pt, &b_mc_truth_tW2_pt);
    if( fChain->GetBranch("mc_truth_tWl1_pt") ) fChain->SetBranchAddress("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, &b_mc_truth_tWl1_pt);
    if( fChain->GetBranch("mc_truth_tWl2_pt") ) fChain->SetBranchAddress("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, &b_mc_truth_tWl2_pt);
    if( fChain->GetBranch("mc_truth_tWtau1_pt") ) fChain->SetBranchAddress("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, &b_mc_truth_tWtau1_pt);
    if( fChain->GetBranch("mc_truth_tWtau2_pt") ) fChain->SetBranchAddress("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, &b_mc_truth_tWtau2_pt);
    if( fChain->GetBranch("mc_truth_tWtaul1_pt") ) fChain->SetBranchAddress("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, &b_mc_truth_tWtaul1_pt);
    if( fChain->GetBranch("mc_truth_tWtaul2_pt") ) fChain->SetBranchAddress("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, &b_mc_truth_tWtaul2_pt);
    if( fChain->GetBranch("mc_truth_tWq11_pt") ) fChain->SetBranchAddress("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, &b_mc_truth_tWq11_pt);
    if( fChain->GetBranch("mc_truth_tWq21_pt") ) fChain->SetBranchAddress("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, &b_mc_truth_tWq21_pt);
    if( fChain->GetBranch("mc_truth_tWq12_pt") ) fChain->SetBranchAddress("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, &b_mc_truth_tWq12_pt);
    if( fChain->GetBranch("mc_truth_tWq22_pt") ) fChain->SetBranchAddress("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, &b_mc_truth_tWq22_pt);   

    if( fChain->GetBranch("mc_truth_t_pt") ) fChain->SetBranchAddress("mc_truth_t_pt", &mc_truth_t_pt, &b_mc_truth_t_pt);
    if( fChain->GetBranch("mc_truth_tb_pt") ) fChain->SetBranchAddress("mc_truth_tb_pt", &mc_truth_tb_pt, &b_mc_truth_tb_pt);
    if( fChain->GetBranch("mc_truth_tW_pt") ) fChain->SetBranchAddress("mc_truth_tW_pt", &mc_truth_tW_pt, &b_mc_truth_tW_pt);
    if( fChain->GetBranch("mc_truth_tWl_pt") ) fChain->SetBranchAddress("mc_truth_tWl_pt", &mc_truth_tWl_pt, &b_mc_truth_tWl_pt);
    if( fChain->GetBranch("mc_truth_tWq1_pt") ) fChain->SetBranchAddress("mc_truth_tWq1_pt", &mc_truth_tWq1_pt, &b_mc_truth_tWq1_pt);
    if( fChain->GetBranch("mc_truth_tWq2_pt") ) fChain->SetBranchAddress("mc_truth_tWq2_pt", &mc_truth_tWq2_pt, &b_mc_truth_tWq2_pt);

    if( fChain->GetBranch("mc_truth_W_pt") ) fChain->SetBranchAddress("mc_truth_W_pt", &mc_truth_W_pt, &b_mc_truth_W_pt);
    if( fChain->GetBranch("mc_truth_Wl_pt") ) fChain->SetBranchAddress("mc_truth_Wl_pt", &mc_truth_Wl_pt, &b_mc_truth_Wl_pt);
    if( fChain->GetBranch("mc_truth_Wtau_pt") ) fChain->SetBranchAddress("mc_truth_Wtau_pt", &mc_truth_Wtau_pt, &b_mc_truth_Wtau_pt);
    if( fChain->GetBranch("mc_truth_Wtaul_pt") ) fChain->SetBranchAddress("mc_truth_Wtaul_pt", &mc_truth_Wtaul_pt, &b_mc_truth_Wtaul_pt);
    if( fChain->GetBranch("mc_truth_Wq1_pt") ) fChain->SetBranchAddress("mc_truth_Wq1_pt", &mc_truth_Wq1_pt, &b_mc_truth_Wq1_pt);
    if( fChain->GetBranch("mc_truth_Wq2_pt") ) fChain->SetBranchAddress("mc_truth_Wq2_pt", &mc_truth_Wq2_pt, &b_mc_truth_Wq2_pt);
    if( fChain->GetBranch("mc_truth_Z_pt") ) fChain->SetBranchAddress("mc_truth_Z_pt", &mc_truth_Z_pt, &b_mc_truth_Z_pt);
    if( fChain->GetBranch("mc_truth_Zl1_pt") ) fChain->SetBranchAddress("mc_truth_Zl1_pt", &mc_truth_Zl1_pt, &b_mc_truth_Zl1_pt);
    if( fChain->GetBranch("mc_truth_Zl2_pt") ) fChain->SetBranchAddress("mc_truth_Zl2_pt", &mc_truth_Zl2_pt, &b_mc_truth_Zl2_pt);
    if( fChain->GetBranch("mc_truth_Ztau1_pt") ) fChain->SetBranchAddress("mc_truth_Ztau1_pt", &mc_truth_Ztau1_pt, &b_mc_truth_Ztau1_pt);
    if( fChain->GetBranch("mc_truth_Ztau2_pt") ) fChain->SetBranchAddress("mc_truth_Ztau2_pt", &mc_truth_Ztau2_pt, &b_mc_truth_Ztau2_pt);
    if( fChain->GetBranch("mc_truth_Ztaul1_pt") ) fChain->SetBranchAddress("mc_truth_Ztaul1_pt", &mc_truth_Ztaul1_pt, &b_mc_truth_Ztaul1_pt);
    if( fChain->GetBranch("mc_truth_Ztaul2_pt") ) fChain->SetBranchAddress("mc_truth_Ztaul2_pt", &mc_truth_Ztaul2_pt, &b_mc_truth_Ztaul2_pt);    
    if( fChain->GetBranch("mc_truth_Zq1_pt") ) fChain->SetBranchAddress("mc_truth_Zq1_pt", &mc_truth_Zq1_pt, &b_mc_truth_Zq1_pt);
    if( fChain->GetBranch("mc_truth_Zq2_pt") ) fChain->SetBranchAddress("mc_truth_Zq2_pt", &mc_truth_Zq2_pt, &b_mc_truth_Zq2_pt);    
    if( fChain->GetBranch("mc_truth_gammal1_pt") ) fChain->SetBranchAddress("mc_truth_gammal1_pt", &mc_truth_gammal1_pt, &b_mc_truth_gammal1_pt);
    if( fChain->GetBranch("mc_truth_gammal2_pt") ) fChain->SetBranchAddress("mc_truth_gammal2_pt", &mc_truth_gammal2_pt, &b_mc_truth_gammal2_pt);
    if( fChain->GetBranch("mc_truth_gammatau1_pt") ) fChain->SetBranchAddress("mc_truth_gammatau1_pt", &mc_truth_gammatau1_pt, &b_mc_truth_gammatau1_pt);
    if( fChain->GetBranch("mc_truth_gammatau2_pt") ) fChain->SetBranchAddress("mc_truth_gammatau2_pt", &mc_truth_gammatau2_pt, &b_mc_truth_gammatau2_pt);
    if( fChain->GetBranch("mc_truth_gammataul1_pt") ) fChain->SetBranchAddress("mc_truth_gammataul1_pt", &mc_truth_gammataul1_pt, &b_mc_truth_gammataul1_pt);
    if( fChain->GetBranch("mc_truth_gammataul2_pt") ) fChain->SetBranchAddress("mc_truth_gammataul2_pt", &mc_truth_gammataul2_pt, &b_mc_truth_gammataul2_pt);
    if( fChain->GetBranch("mc_truth_gamma_pt") ) fChain->SetBranchAddress("mc_truth_gamma_pt", &mc_truth_gamma_pt, &b_mc_truth_gamma_pt);

    if( fChain->GetBranch("mc_truth_h0_eta") ) fChain->SetBranchAddress("mc_truth_h0_eta", &mc_truth_h0_eta, &b_mc_truth_h0_eta);
    if( fChain->GetBranch("mc_truth_h0W1_eta") ) fChain->SetBranchAddress("mc_truth_h0W1_eta", &mc_truth_h0W1_eta, &b_mc_truth_h0W1_eta);
    if( fChain->GetBranch("mc_truth_h0Wl1_eta") ) fChain->SetBranchAddress("mc_truth_h0Wl1_eta", &mc_truth_h0Wl1_eta, &b_mc_truth_h0Wl1_eta);
    if( fChain->GetBranch("mc_truth_h0Wtau1_eta") ) fChain->SetBranchAddress("mc_truth_h0Wtau1_eta", &mc_truth_h0Wtau1_eta, &b_mc_truth_h0Wtau1_eta);
    if( fChain->GetBranch("mc_truth_h0Wtaul1_eta") ) fChain->SetBranchAddress("mc_truth_h0Wtaul1_eta", &mc_truth_h0Wtaul1_eta, &b_mc_truth_h0Wtaul1_eta);
    if( fChain->GetBranch("mc_truth_h0Wq11_eta") ) fChain->SetBranchAddress("mc_truth_h0Wq11_eta", &mc_truth_h0Wq11_eta, &b_mc_truth_h0Wq11_eta);
    if( fChain->GetBranch("mc_truth_h0Wq21_eta") ) fChain->SetBranchAddress("mc_truth_h0Wq21_eta", &mc_truth_h0Wq21_eta, &b_mc_truth_h0Wq21_eta);
    if( fChain->GetBranch("mc_truth_h0W2_eta") ) fChain->SetBranchAddress("mc_truth_h0W2_eta", &mc_truth_h0W2_eta, &b_mc_truth_h0W2_eta);
    if( fChain->GetBranch("mc_truth_h0Wl2_eta") ) fChain->SetBranchAddress("mc_truth_h0Wl2_eta", &mc_truth_h0Wl2_eta, &b_mc_truth_h0Wl2_eta);
    if( fChain->GetBranch("mc_truth_h0Wtau2_eta") ) fChain->SetBranchAddress("mc_truth_h0Wtau2_eta", &mc_truth_h0Wtau2_eta, &b_mc_truth_h0Wtau2_eta);
    if( fChain->GetBranch("mc_truth_h0Wtaul2_eta") ) fChain->SetBranchAddress("mc_truth_h0Wtaul2_eta", &mc_truth_h0Wtaul2_eta, &b_mc_truth_h0Wtaul2_eta);
    if( fChain->GetBranch("mc_truth_h0Wq12_eta") ) fChain->SetBranchAddress("mc_truth_h0Wq12_eta", &mc_truth_h0Wq12_eta, &b_mc_truth_h0Wq12_eta);
    if( fChain->GetBranch("mc_truth_h0Wq22_eta") ) fChain->SetBranchAddress("mc_truth_h0Wq22_eta", &mc_truth_h0Wq22_eta, &b_mc_truth_h0Wq22_eta);
    if( fChain->GetBranch("mc_truth_h0Z1_eta") ) fChain->SetBranchAddress("mc_truth_h0Z1_eta", &mc_truth_h0Z1_eta, &b_mc_truth_h0Z1_eta);
    if( fChain->GetBranch("mc_truth_h0Zl11_eta") ) fChain->SetBranchAddress("mc_truth_h0Zl11_eta", &mc_truth_h0Zl11_eta, &b_mc_truth_h0Zl11_eta);
    if( fChain->GetBranch("mc_truth_h0Zl21_eta") ) fChain->SetBranchAddress("mc_truth_h0Zl21_eta", &mc_truth_h0Zl21_eta, &b_mc_truth_h0Zl21_eta);
    if( fChain->GetBranch("mc_truth_h0Zq11_eta") ) fChain->SetBranchAddress("mc_truth_h0Zq11_eta", &mc_truth_h0Zq11_eta, &b_mc_truth_h0Zq11_eta);
    if( fChain->GetBranch("mc_truth_h0Zq21_eta") ) fChain->SetBranchAddress("mc_truth_h0Zq21_eta", &mc_truth_h0Zq21_eta, &b_mc_truth_h0Zq21_eta);
    if( fChain->GetBranch("mc_truth_h0Z2_eta") ) fChain->SetBranchAddress("mc_truth_h0Z2_eta", &mc_truth_h0Z2_eta, &b_mc_truth_h0Z2_eta);
    if( fChain->GetBranch("mc_truth_h0Zl12_eta") ) fChain->SetBranchAddress("mc_truth_h0Zl12_eta", &mc_truth_h0Zl12_eta, &b_mc_truth_h0Zl12_eta);
    if( fChain->GetBranch("mc_truth_h0Zl22_eta") ) fChain->SetBranchAddress("mc_truth_h0Zl22_eta", &mc_truth_h0Zl22_eta, &b_mc_truth_h0Zl22_eta);
    if( fChain->GetBranch("mc_truth_h0Zq12_eta") ) fChain->SetBranchAddress("mc_truth_h0Zq12_eta", &mc_truth_h0Zq12_eta, &b_mc_truth_h0Zq12_eta);
    if( fChain->GetBranch("mc_truth_h0Zq22_eta") ) fChain->SetBranchAddress("mc_truth_h0Zq22_eta", &mc_truth_h0Zq22_eta, &b_mc_truth_h0Zq22_eta);
    if( fChain->GetBranch("mc_truth_h0tau1_eta") ) fChain->SetBranchAddress("mc_truth_h0tau1_eta", &mc_truth_h0tau1_eta, &b_mc_truth_h0tau1_eta);
    if( fChain->GetBranch("mc_truth_h0tau2_eta") ) fChain->SetBranchAddress("mc_truth_h0tau2_eta", &mc_truth_h0tau2_eta, &b_mc_truth_h0tau2_eta);
    if( fChain->GetBranch("mc_truth_h0taul1_eta") ) fChain->SetBranchAddress("mc_truth_h0taul1_eta", &mc_truth_h0taul1_eta, &b_mc_truth_h0taul1_eta);
    if( fChain->GetBranch("mc_truth_h0taul2_eta") ) fChain->SetBranchAddress("mc_truth_h0taul2_eta", &mc_truth_h0taul2_eta, &b_mc_truth_h0taul2_eta);
    if( fChain->GetBranch("mc_truth_h0b1_eta") ) fChain->SetBranchAddress("mc_truth_h0b1_eta", &mc_truth_h0b1_eta, &b_mc_truth_h0b1_eta);
    if( fChain->GetBranch("mc_truth_h0b2_eta") ) fChain->SetBranchAddress("mc_truth_h0b2_eta", &mc_truth_h0b2_eta, &b_mc_truth_h0b2_eta);

    if( fChain->GetBranch("mc_truth_t1_eta") ) fChain->SetBranchAddress("mc_truth_t1_eta", &mc_truth_t1_eta, &b_mc_truth_t1_eta);
    if( fChain->GetBranch("mc_truth_t2_eta") ) fChain->SetBranchAddress("mc_truth_t2_eta", &mc_truth_t2_eta, &b_mc_truth_t2_eta);
    if( fChain->GetBranch("mc_truth_tb1_eta") ) fChain->SetBranchAddress("mc_truth_tb1_eta", &mc_truth_tb1_eta, &b_mc_truth_tb1_eta);
    if( fChain->GetBranch("mc_truth_tb2_eta") ) fChain->SetBranchAddress("mc_truth_tb2_eta", &mc_truth_tb2_eta, &b_mc_truth_tb2_eta);
    if( fChain->GetBranch("mc_truth_tW1_eta") ) fChain->SetBranchAddress("mc_truth_tW1_eta", &mc_truth_tW1_eta, &b_mc_truth_tW1_eta);
    if( fChain->GetBranch("mc_truth_tW2_eta") ) fChain->SetBranchAddress("mc_truth_tW2_eta", &mc_truth_tW2_eta, &b_mc_truth_tW2_eta);
    if( fChain->GetBranch("mc_truth_tWl1_eta") ) fChain->SetBranchAddress("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, &b_mc_truth_tWl1_eta);
    if( fChain->GetBranch("mc_truth_tWl2_eta") ) fChain->SetBranchAddress("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, &b_mc_truth_tWl2_eta);
    if( fChain->GetBranch("mc_truth_tWtau1_eta") ) fChain->SetBranchAddress("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, &b_mc_truth_tWtau1_eta);
    if( fChain->GetBranch("mc_truth_tWtau2_eta") ) fChain->SetBranchAddress("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, &b_mc_truth_tWtau2_eta);
    if( fChain->GetBranch("mc_truth_tWtaul1_eta") ) fChain->SetBranchAddress("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, &b_mc_truth_tWtaul1_eta);
    if( fChain->GetBranch("mc_truth_tWtaul2_eta") ) fChain->SetBranchAddress("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, &b_mc_truth_tWtaul2_eta);
    if( fChain->GetBranch("mc_truth_tWq11_eta") ) fChain->SetBranchAddress("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, &b_mc_truth_tWq11_eta);
    if( fChain->GetBranch("mc_truth_tWq21_eta") ) fChain->SetBranchAddress("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, &b_mc_truth_tWq21_eta);
    if( fChain->GetBranch("mc_truth_tWq12_eta") ) fChain->SetBranchAddress("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, &b_mc_truth_tWq12_eta);
    if( fChain->GetBranch("mc_truth_tWq22_eta") ) fChain->SetBranchAddress("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, &b_mc_truth_tWq22_eta);   

    if( fChain->GetBranch("mc_truth_t_eta") ) fChain->SetBranchAddress("mc_truth_t_eta", &mc_truth_t_eta, &b_mc_truth_t_eta);
    if( fChain->GetBranch("mc_truth_tb_eta") ) fChain->SetBranchAddress("mc_truth_tb_eta", &mc_truth_tb_eta, &b_mc_truth_tb_eta);
    if( fChain->GetBranch("mc_truth_tW_eta") ) fChain->SetBranchAddress("mc_truth_tW_eta", &mc_truth_tW_eta, &b_mc_truth_tW_eta);
    if( fChain->GetBranch("mc_truth_tWl_eta") ) fChain->SetBranchAddress("mc_truth_tWl_eta", &mc_truth_tWl_eta, &b_mc_truth_tWl_eta);
    if( fChain->GetBranch("mc_truth_tWq1_eta") ) fChain->SetBranchAddress("mc_truth_tWq1_eta", &mc_truth_tWq1_eta, &b_mc_truth_tWq1_eta);
    if( fChain->GetBranch("mc_truth_tWq2_eta") ) fChain->SetBranchAddress("mc_truth_tWq2_eta", &mc_truth_tWq2_eta, &b_mc_truth_tWq2_eta);

    if( fChain->GetBranch("mc_truth_W_eta") ) fChain->SetBranchAddress("mc_truth_W_eta", &mc_truth_W_eta, &b_mc_truth_W_eta);
    if( fChain->GetBranch("mc_truth_Wl_eta") ) fChain->SetBranchAddress("mc_truth_Wl_eta", &mc_truth_Wl_eta, &b_mc_truth_Wl_eta);
    if( fChain->GetBranch("mc_truth_Wtau_eta") ) fChain->SetBranchAddress("mc_truth_Wtau_eta", &mc_truth_Wtau_eta, &b_mc_truth_Wtau_eta);
    if( fChain->GetBranch("mc_truth_Wtaul_eta") ) fChain->SetBranchAddress("mc_truth_Wtaul_eta", &mc_truth_Wtaul_eta, &b_mc_truth_Wtaul_eta);
    if( fChain->GetBranch("mc_truth_Wq1_eta") ) fChain->SetBranchAddress("mc_truth_Wq1_eta", &mc_truth_Wq1_eta, &b_mc_truth_Wq1_eta);
    if( fChain->GetBranch("mc_truth_Wq2_eta") ) fChain->SetBranchAddress("mc_truth_Wq2_eta", &mc_truth_Wq2_eta, &b_mc_truth_Wq2_eta);
    if( fChain->GetBranch("mc_truth_Z_eta") ) fChain->SetBranchAddress("mc_truth_Z_eta", &mc_truth_Z_eta, &b_mc_truth_Z_eta);
    if( fChain->GetBranch("mc_truth_Zl1_eta") ) fChain->SetBranchAddress("mc_truth_Zl1_eta", &mc_truth_Zl1_eta, &b_mc_truth_Zl1_eta);
    if( fChain->GetBranch("mc_truth_Zl2_eta") ) fChain->SetBranchAddress("mc_truth_Zl2_eta", &mc_truth_Zl2_eta, &b_mc_truth_Zl2_eta);
    if( fChain->GetBranch("mc_truth_Ztau1_eta") ) fChain->SetBranchAddress("mc_truth_Ztau1_eta", &mc_truth_Ztau1_eta, &b_mc_truth_Ztau1_eta);
    if( fChain->GetBranch("mc_truth_Ztau2_eta") ) fChain->SetBranchAddress("mc_truth_Ztau2_eta", &mc_truth_Ztau2_eta, &b_mc_truth_Ztau2_eta);
    if( fChain->GetBranch("mc_truth_Ztaul1_eta") ) fChain->SetBranchAddress("mc_truth_Ztaul1_eta", &mc_truth_Ztaul1_eta, &b_mc_truth_Ztaul1_eta);
    if( fChain->GetBranch("mc_truth_Ztaul2_eta") ) fChain->SetBranchAddress("mc_truth_Ztaul2_eta", &mc_truth_Ztaul2_eta, &b_mc_truth_Ztaul2_eta);    
    if( fChain->GetBranch("mc_truth_Zq1_eta") ) fChain->SetBranchAddress("mc_truth_Zq1_eta", &mc_truth_Zq1_eta, &b_mc_truth_Zq1_eta);
    if( fChain->GetBranch("mc_truth_Zq2_eta") ) fChain->SetBranchAddress("mc_truth_Zq2_eta", &mc_truth_Zq2_eta, &b_mc_truth_Zq2_eta);
    if( fChain->GetBranch("mc_truth_gammal1_eta") ) fChain->SetBranchAddress("mc_truth_gammal1_eta", &mc_truth_gammal1_eta, &b_mc_truth_gammal1_eta);
    if( fChain->GetBranch("mc_truth_gammal2_eta") ) fChain->SetBranchAddress("mc_truth_gammal2_eta", &mc_truth_gammal2_eta, &b_mc_truth_gammal2_eta);
    if( fChain->GetBranch("mc_truth_gammatau1_eta") ) fChain->SetBranchAddress("mc_truth_gammatau1_eta", &mc_truth_gammatau1_eta, &b_mc_truth_gammatau1_eta);
    if( fChain->GetBranch("mc_truth_gammatau2_eta") ) fChain->SetBranchAddress("mc_truth_gammatau2_eta", &mc_truth_gammatau2_eta, &b_mc_truth_gammatau2_eta);
    if( fChain->GetBranch("mc_truth_gammataul1_eta") ) fChain->SetBranchAddress("mc_truth_gammataul1_eta", &mc_truth_gammataul1_eta, &b_mc_truth_gammataul1_eta);
    if( fChain->GetBranch("mc_truth_gammataul2_eta") ) fChain->SetBranchAddress("mc_truth_gammataul2_eta", &mc_truth_gammataul2_eta, &b_mc_truth_gammataul2_eta);
    if( fChain->GetBranch("mc_truth_gamma_eta") ) fChain->SetBranchAddress("mc_truth_gamma_eta", &mc_truth_gamma_eta, &b_mc_truth_gamma_eta);

    if( fChain->GetBranch("mc_truth_h0_phi") ) fChain->SetBranchAddress("mc_truth_h0_phi", &mc_truth_h0_phi, &b_mc_truth_h0_phi);
    if( fChain->GetBranch("mc_truth_h0W1_phi") ) fChain->SetBranchAddress("mc_truth_h0W1_phi", &mc_truth_h0W1_phi, &b_mc_truth_h0W1_phi);
    if( fChain->GetBranch("mc_truth_h0Wl1_phi") ) fChain->SetBranchAddress("mc_truth_h0Wl1_phi", &mc_truth_h0Wl1_phi, &b_mc_truth_h0Wl1_phi);
    if( fChain->GetBranch("mc_truth_h0Wtau1_phi") ) fChain->SetBranchAddress("mc_truth_h0Wtau1_phi", &mc_truth_h0Wtau1_phi, &b_mc_truth_h0Wtau1_phi);
    if( fChain->GetBranch("mc_truth_h0Wq11_phi") ) fChain->SetBranchAddress("mc_truth_h0Wq11_phi", &mc_truth_h0Wq11_phi, &b_mc_truth_h0Wq11_phi);
    if( fChain->GetBranch("mc_truth_h0Wq21_phi") ) fChain->SetBranchAddress("mc_truth_h0Wq21_phi", &mc_truth_h0Wq21_phi, &b_mc_truth_h0Wq21_phi);
    if( fChain->GetBranch("mc_truth_h0W2_phi") ) fChain->SetBranchAddress("mc_truth_h0W2_phi", &mc_truth_h0W2_phi, &b_mc_truth_h0W2_phi);
    if( fChain->GetBranch("mc_truth_h0Wl2_phi") ) fChain->SetBranchAddress("mc_truth_h0Wl2_phi", &mc_truth_h0Wl2_phi, &b_mc_truth_h0Wl2_phi);
    if( fChain->GetBranch("mc_truth_h0Wq12_phi") ) fChain->SetBranchAddress("mc_truth_h0Wq12_phi", &mc_truth_h0Wq12_phi, &b_mc_truth_h0Wq12_phi);
    if( fChain->GetBranch("mc_truth_h0Wq22_phi") ) fChain->SetBranchAddress("mc_truth_h0Wq22_phi", &mc_truth_h0Wq22_phi, &b_mc_truth_h0Wq22_phi);
    if( fChain->GetBranch("mc_truth_h0Z1_phi") ) fChain->SetBranchAddress("mc_truth_h0Z1_phi", &mc_truth_h0Z1_phi, &b_mc_truth_h0Z1_phi);
    if( fChain->GetBranch("mc_truth_h0Zl11_phi") ) fChain->SetBranchAddress("mc_truth_h0Zl11_phi", &mc_truth_h0Zl11_phi, &b_mc_truth_h0Zl11_phi);
    if( fChain->GetBranch("mc_truth_h0Zl21_phi") ) fChain->SetBranchAddress("mc_truth_h0Zl21_phi", &mc_truth_h0Zl21_phi, &b_mc_truth_h0Zl21_phi);
    if( fChain->GetBranch("mc_truth_h0Zq11_phi") ) fChain->SetBranchAddress("mc_truth_h0Zq11_phi", &mc_truth_h0Zq11_phi, &b_mc_truth_h0Zq11_phi);
    if( fChain->GetBranch("mc_truth_h0Zq21_phi") ) fChain->SetBranchAddress("mc_truth_h0Zq21_phi", &mc_truth_h0Zq21_phi, &b_mc_truth_h0Zq21_phi);
    if( fChain->GetBranch("mc_truth_h0Z2_phi") ) fChain->SetBranchAddress("mc_truth_h0Z2_phi", &mc_truth_h0Z2_phi, &b_mc_truth_h0Z2_phi);
    if( fChain->GetBranch("mc_truth_h0Zl12_phi") ) fChain->SetBranchAddress("mc_truth_h0Zl12_phi", &mc_truth_h0Zl12_phi, &b_mc_truth_h0Zl12_phi);
    if( fChain->GetBranch("mc_truth_h0Zl22_phi") ) fChain->SetBranchAddress("mc_truth_h0Zl22_phi", &mc_truth_h0Zl22_phi, &b_mc_truth_h0Zl22_phi);
    if( fChain->GetBranch("mc_truth_h0Zq12_phi") ) fChain->SetBranchAddress("mc_truth_h0Zq12_phi", &mc_truth_h0Zq12_phi, &b_mc_truth_h0Zq12_phi);
    if( fChain->GetBranch("mc_truth_h0Zq22_phi") ) fChain->SetBranchAddress("mc_truth_h0Zq22_phi", &mc_truth_h0Zq22_phi, &b_mc_truth_h0Zq22_phi);
    if( fChain->GetBranch("mc_truth_h0tau1_phi") ) fChain->SetBranchAddress("mc_truth_h0tau1_phi", &mc_truth_h0tau1_phi, &b_mc_truth_h0tau1_phi);
    if( fChain->GetBranch("mc_truth_h0tau2_phi") ) fChain->SetBranchAddress("mc_truth_h0tau2_phi", &mc_truth_h0tau2_phi, &b_mc_truth_h0tau2_phi);
    if( fChain->GetBranch("mc_truth_h0taul1_phi") ) fChain->SetBranchAddress("mc_truth_h0taul1_phi", &mc_truth_h0taul1_phi, &b_mc_truth_h0taul1_phi);
    if( fChain->GetBranch("mc_truth_h0taul2_phi") ) fChain->SetBranchAddress("mc_truth_h0taul2_phi", &mc_truth_h0taul2_phi, &b_mc_truth_h0taul2_phi);
    if( fChain->GetBranch("mc_truth_h0b1_phi") ) fChain->SetBranchAddress("mc_truth_h0b1_phi", &mc_truth_h0b1_phi, &b_mc_truth_h0b1_phi);
    if( fChain->GetBranch("mc_truth_h0b2_phi") ) fChain->SetBranchAddress("mc_truth_h0b2_phi", &mc_truth_h0b2_phi, &b_mc_truth_h0b2_phi);

    if( fChain->GetBranch("mc_truth_t1_phi") ) fChain->SetBranchAddress("mc_truth_t1_phi", &mc_truth_t1_phi, &b_mc_truth_t1_phi);
    if( fChain->GetBranch("mc_truth_t2_phi") ) fChain->SetBranchAddress("mc_truth_t2_phi", &mc_truth_t2_phi, &b_mc_truth_t2_phi);
    if( fChain->GetBranch("mc_truth_tb1_phi") ) fChain->SetBranchAddress("mc_truth_tb1_phi", &mc_truth_tb1_phi, &b_mc_truth_tb1_phi);
    if( fChain->GetBranch("mc_truth_tb2_phi") ) fChain->SetBranchAddress("mc_truth_tb2_phi", &mc_truth_tb2_phi, &b_mc_truth_tb2_phi);
    if( fChain->GetBranch("mc_truth_tW1_phi") ) fChain->SetBranchAddress("mc_truth_tW1_phi", &mc_truth_tW1_phi, &b_mc_truth_tW1_phi);
    if( fChain->GetBranch("mc_truth_tW2_phi") ) fChain->SetBranchAddress("mc_truth_tW2_phi", &mc_truth_tW2_phi, &b_mc_truth_tW2_phi);
    if( fChain->GetBranch("mc_truth_tWl1_phi") ) fChain->SetBranchAddress("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, &b_mc_truth_tWl1_phi);
    if( fChain->GetBranch("mc_truth_tWl2_phi") ) fChain->SetBranchAddress("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, &b_mc_truth_tWl2_phi);
    if( fChain->GetBranch("mc_truth_tWtau1_phi") ) fChain->SetBranchAddress("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, &b_mc_truth_tWtau1_phi);
    if( fChain->GetBranch("mc_truth_tWtau2_phi") ) fChain->SetBranchAddress("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, &b_mc_truth_tWtau2_phi);
    if( fChain->GetBranch("mc_truth_tWtaul1_phi") ) fChain->SetBranchAddress("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, &b_mc_truth_tWtaul1_phi);
    if( fChain->GetBranch("mc_truth_tWtaul2_phi") ) fChain->SetBranchAddress("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, &b_mc_truth_tWtaul2_phi);
    if( fChain->GetBranch("mc_truth_tWq11_phi") ) fChain->SetBranchAddress("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, &b_mc_truth_tWq11_phi);
    if( fChain->GetBranch("mc_truth_tWq21_phi") ) fChain->SetBranchAddress("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, &b_mc_truth_tWq21_phi);
    if( fChain->GetBranch("mc_truth_tWq12_phi") ) fChain->SetBranchAddress("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, &b_mc_truth_tWq12_phi);
    if( fChain->GetBranch("mc_truth_tWq22_phi") ) fChain->SetBranchAddress("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, &b_mc_truth_tWq22_phi);   

    if( fChain->GetBranch("mc_truth_t_phi") ) fChain->SetBranchAddress("mc_truth_t_phi", &mc_truth_t_phi, &b_mc_truth_t_phi);
    if( fChain->GetBranch("mc_truth_tb_phi") ) fChain->SetBranchAddress("mc_truth_tb_phi", &mc_truth_tb_phi, &b_mc_truth_tb_phi);
    if( fChain->GetBranch("mc_truth_tW_phi") ) fChain->SetBranchAddress("mc_truth_tW_phi", &mc_truth_tW_phi, &b_mc_truth_tW_phi);
    if( fChain->GetBranch("mc_truth_tWl_phi") ) fChain->SetBranchAddress("mc_truth_tWl_phi", &mc_truth_tWl_phi, &b_mc_truth_tWl_phi);
    if( fChain->GetBranch("mc_truth_tWq1_phi") ) fChain->SetBranchAddress("mc_truth_tWq1_phi", &mc_truth_tWq1_phi, &b_mc_truth_tWq1_phi);
    if( fChain->GetBranch("mc_truth_tWq2_phi") ) fChain->SetBranchAddress("mc_truth_tWq2_phi", &mc_truth_tWq2_phi, &b_mc_truth_tWq2_phi);

    if( fChain->GetBranch("mc_truth_W_phi") ) fChain->SetBranchAddress("mc_truth_W_phi", &mc_truth_W_phi, &b_mc_truth_W_phi);
    if( fChain->GetBranch("mc_truth_Wl_phi") ) fChain->SetBranchAddress("mc_truth_Wl_phi", &mc_truth_Wl_phi, &b_mc_truth_Wl_phi);
    if( fChain->GetBranch("mc_truth_Wtau_phi") ) fChain->SetBranchAddress("mc_truth_Wtau_phi", &mc_truth_Wtau_phi, &b_mc_truth_Wtau_phi);
    if( fChain->GetBranch("mc_truth_Wtaul_phi") ) fChain->SetBranchAddress("mc_truth_Wtaul_phi", &mc_truth_Wtaul_phi, &b_mc_truth_Wtaul_phi);
    if( fChain->GetBranch("mc_truth_Wq1_phi") ) fChain->SetBranchAddress("mc_truth_Wq1_phi", &mc_truth_Wq1_phi, &b_mc_truth_Wq1_phi);
    if( fChain->GetBranch("mc_truth_Wq2_phi") ) fChain->SetBranchAddress("mc_truth_Wq2_phi", &mc_truth_Wq2_phi, &b_mc_truth_Wq2_phi);
    if( fChain->GetBranch("mc_truth_Z_phi") ) fChain->SetBranchAddress("mc_truth_Z_phi", &mc_truth_Z_phi, &b_mc_truth_Z_phi);
    if( fChain->GetBranch("mc_truth_Zl1_phi") ) fChain->SetBranchAddress("mc_truth_Zl1_phi", &mc_truth_Zl1_phi, &b_mc_truth_Zl1_phi);
    if( fChain->GetBranch("mc_truth_Zl2_phi") ) fChain->SetBranchAddress("mc_truth_Zl2_phi", &mc_truth_Zl2_phi, &b_mc_truth_Zl2_phi);
    if( fChain->GetBranch("mc_truth_Ztau1_phi") ) fChain->SetBranchAddress("mc_truth_Ztau1_phi", &mc_truth_Ztau1_phi, &b_mc_truth_Ztau1_phi);
    if( fChain->GetBranch("mc_truth_Ztau2_phi") ) fChain->SetBranchAddress("mc_truth_Ztau2_phi", &mc_truth_Ztau2_phi, &b_mc_truth_Ztau2_phi);
    if( fChain->GetBranch("mc_truth_Ztaul1_phi") ) fChain->SetBranchAddress("mc_truth_Ztaul1_phi", &mc_truth_Ztaul1_phi, &b_mc_truth_Ztaul1_phi);
    if( fChain->GetBranch("mc_truth_Ztaul2_phi") ) fChain->SetBranchAddress("mc_truth_Ztaul2_phi", &mc_truth_Ztaul2_phi, &b_mc_truth_Ztaul2_phi);    
    if( fChain->GetBranch("mc_truth_Zq1_phi") ) fChain->SetBranchAddress("mc_truth_Zq1_phi", &mc_truth_Zq1_phi, &b_mc_truth_Zq1_phi);
    if( fChain->GetBranch("mc_truth_Zq2_phi") ) fChain->SetBranchAddress("mc_truth_Zq2_phi", &mc_truth_Zq2_phi, &b_mc_truth_Zq2_phi);
    if( fChain->GetBranch("mc_truth_gammal1_phi") ) fChain->SetBranchAddress("mc_truth_gammal1_phi", &mc_truth_gammal1_phi, &b_mc_truth_gammal1_phi);
    if( fChain->GetBranch("mc_truth_gammal2_phi") ) fChain->SetBranchAddress("mc_truth_gammal2_phi", &mc_truth_gammal2_phi, &b_mc_truth_gammal2_phi);
    if( fChain->GetBranch("mc_truth_gammatau1_phi") ) fChain->SetBranchAddress("mc_truth_gammatau1_phi", &mc_truth_gammatau1_phi, &b_mc_truth_gammatau1_phi);
    if( fChain->GetBranch("mc_truth_gammatau2_phi") ) fChain->SetBranchAddress("mc_truth_gammatau2_phi", &mc_truth_gammatau2_phi, &b_mc_truth_gammatau2_phi);
    if( fChain->GetBranch("mc_truth_gammataul1_phi") ) fChain->SetBranchAddress("mc_truth_gammataul1_phi", &mc_truth_gammataul1_phi, &b_mc_truth_gammataul1_phi);
    if( fChain->GetBranch("mc_truth_gammataul2_phi") ) fChain->SetBranchAddress("mc_truth_gammataul2_phi", &mc_truth_gammataul2_phi, &b_mc_truth_gammataul2_phi);
    if( fChain->GetBranch("mc_truth_gamma_phi") ) fChain->SetBranchAddress("mc_truth_gamma_phi", &mc_truth_gamma_phi, &b_mc_truth_gamma_phi);

    if( fChain->GetBranch("mc_truth_h0_E") ) fChain->SetBranchAddress("mc_truth_h0_E", &mc_truth_h0_E, &b_mc_truth_h0_E);
    if( fChain->GetBranch("mc_truth_h0W1_E") ) fChain->SetBranchAddress("mc_truth_h0W1_E", &mc_truth_h0W1_E, &b_mc_truth_h0W1_E);
    if( fChain->GetBranch("mc_truth_h0Wl1_E") ) fChain->SetBranchAddress("mc_truth_h0Wl1_E", &mc_truth_h0Wl1_E, &b_mc_truth_h0Wl1_E);
    if( fChain->GetBranch("mc_truth_h0Wtau1_E") ) fChain->SetBranchAddress("mc_truth_h0Wtau1_E", &mc_truth_h0Wtau1_E, &b_mc_truth_h0Wtau1_E);
    if( fChain->GetBranch("mc_truth_h0Wtaul1_E") ) fChain->SetBranchAddress("mc_truth_h0Wtaul1_E", &mc_truth_h0Wtaul1_E, &b_mc_truth_h0Wtaul1_E);
    if( fChain->GetBranch("mc_truth_h0Wq11_E") ) fChain->SetBranchAddress("mc_truth_h0Wq11_E", &mc_truth_h0Wq11_E, &b_mc_truth_h0Wq11_E);
    if( fChain->GetBranch("mc_truth_h0Wq21_E") ) fChain->SetBranchAddress("mc_truth_h0Wq21_E", &mc_truth_h0Wq21_E, &b_mc_truth_h0Wq21_E);
    if( fChain->GetBranch("mc_truth_h0W2_E") ) fChain->SetBranchAddress("mc_truth_h0W2_E", &mc_truth_h0W2_E, &b_mc_truth_h0W2_E);
    if( fChain->GetBranch("mc_truth_h0Wl2_E") ) fChain->SetBranchAddress("mc_truth_h0Wl2_E", &mc_truth_h0Wl2_E, &b_mc_truth_h0Wl2_E);
    if( fChain->GetBranch("mc_truth_h0Wtau2_E") ) fChain->SetBranchAddress("mc_truth_h0Wtau2_E", &mc_truth_h0Wtau2_E, &b_mc_truth_h0Wtau2_E);
    if( fChain->GetBranch("mc_truth_h0Wtaul2_E") ) fChain->SetBranchAddress("mc_truth_h0Wtaul2_E", &mc_truth_h0Wtaul2_E, &b_mc_truth_h0Wtaul2_E);   
    if( fChain->GetBranch("mc_truth_h0Wq12_E") ) fChain->SetBranchAddress("mc_truth_h0Wq12_E", &mc_truth_h0Wq12_E, &b_mc_truth_h0Wq12_E);
    if( fChain->GetBranch("mc_truth_h0Wq22_E") ) fChain->SetBranchAddress("mc_truth_h0Wq22_E", &mc_truth_h0Wq22_E, &b_mc_truth_h0Wq22_E);
    if( fChain->GetBranch("mc_truth_h0Z1_E") ) fChain->SetBranchAddress("mc_truth_h0Z1_E", &mc_truth_h0Z1_E, &b_mc_truth_h0Z1_E);
    if( fChain->GetBranch("mc_truth_h0Zl11_E") ) fChain->SetBranchAddress("mc_truth_h0Zl11_E", &mc_truth_h0Zl11_E, &b_mc_truth_h0Zl11_E);
    if( fChain->GetBranch("mc_truth_h0Zl21_E") ) fChain->SetBranchAddress("mc_truth_h0Zl21_E", &mc_truth_h0Zl21_E, &b_mc_truth_h0Zl21_E);
    if( fChain->GetBranch("mc_truth_h0Zq11_E") ) fChain->SetBranchAddress("mc_truth_h0Zq11_E", &mc_truth_h0Zq11_E, &b_mc_truth_h0Zq11_E);
    if( fChain->GetBranch("mc_truth_h0Zq21_E") ) fChain->SetBranchAddress("mc_truth_h0Zq21_E", &mc_truth_h0Zq21_E, &b_mc_truth_h0Zq21_E);
    if( fChain->GetBranch("mc_truth_h0Z2_E") ) fChain->SetBranchAddress("mc_truth_h0Z2_E", &mc_truth_h0Z2_E, &b_mc_truth_h0Z2_E);
    if( fChain->GetBranch("mc_truth_h0Zl12_E") ) fChain->SetBranchAddress("mc_truth_h0Zl12_E", &mc_truth_h0Zl12_E, &b_mc_truth_h0Zl12_E);
    if( fChain->GetBranch("mc_truth_h0Zl22_E") ) fChain->SetBranchAddress("mc_truth_h0Zl22_E", &mc_truth_h0Zl22_E, &b_mc_truth_h0Zl22_E);
    if( fChain->GetBranch("mc_truth_h0Zq12_E") ) fChain->SetBranchAddress("mc_truth_h0Zq12_E", &mc_truth_h0Zq12_E, &b_mc_truth_h0Zq12_E);
    if( fChain->GetBranch("mc_truth_h0Zq22_E") ) fChain->SetBranchAddress("mc_truth_h0Zq22_E", &mc_truth_h0Zq22_E, &b_mc_truth_h0Zq22_E);
    if( fChain->GetBranch("mc_truth_h0tau1_E") ) fChain->SetBranchAddress("mc_truth_h0tau1_E", &mc_truth_h0tau1_E, &b_mc_truth_h0tau1_E);
    if( fChain->GetBranch("mc_truth_h0tau2_E") ) fChain->SetBranchAddress("mc_truth_h0tau2_E", &mc_truth_h0tau2_E, &b_mc_truth_h0tau2_E);
    if( fChain->GetBranch("mc_truth_h0taul1_E") ) fChain->SetBranchAddress("mc_truth_h0taul1_E", &mc_truth_h0taul1_E, &b_mc_truth_h0taul1_E);
    if( fChain->GetBranch("mc_truth_h0taul2_E") ) fChain->SetBranchAddress("mc_truth_h0taul2_E", &mc_truth_h0taul2_E, &b_mc_truth_h0taul2_E);
    if( fChain->GetBranch("mc_truth_h0b1_E") ) fChain->SetBranchAddress("mc_truth_h0b1_E", &mc_truth_h0b1_E, &b_mc_truth_h0b1_E);
    if( fChain->GetBranch("mc_truth_h0b2_E") ) fChain->SetBranchAddress("mc_truth_h0b2_E", &mc_truth_h0b2_E, &b_mc_truth_h0b2_E);

    if( fChain->GetBranch("mc_truth_t1_E") ) fChain->SetBranchAddress("mc_truth_t1_E", &mc_truth_t1_E, &b_mc_truth_t1_E);
    if( fChain->GetBranch("mc_truth_t2_E") ) fChain->SetBranchAddress("mc_truth_t2_E", &mc_truth_t2_E, &b_mc_truth_t2_E);
    if( fChain->GetBranch("mc_truth_tb1_E") ) fChain->SetBranchAddress("mc_truth_tb1_E", &mc_truth_tb1_E, &b_mc_truth_tb1_E);
    if( fChain->GetBranch("mc_truth_tb2_E") ) fChain->SetBranchAddress("mc_truth_tb2_E", &mc_truth_tb2_E, &b_mc_truth_tb2_E);
    if( fChain->GetBranch("mc_truth_tW1_E") ) fChain->SetBranchAddress("mc_truth_tW1_E", &mc_truth_tW1_E, &b_mc_truth_tW1_E);
    if( fChain->GetBranch("mc_truth_tW2_E") ) fChain->SetBranchAddress("mc_truth_tW2_E", &mc_truth_tW2_E, &b_mc_truth_tW2_E);
    if( fChain->GetBranch("mc_truth_tWl1_E") ) fChain->SetBranchAddress("mc_truth_tWl1_E", &mc_truth_tWl1_E, &b_mc_truth_tWl1_E);
    if( fChain->GetBranch("mc_truth_tWl2_E") ) fChain->SetBranchAddress("mc_truth_tWl2_E", &mc_truth_tWl2_E, &b_mc_truth_tWl2_E);
    if( fChain->GetBranch("mc_truth_tWtau1_E") ) fChain->SetBranchAddress("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, &b_mc_truth_tWtau1_E);
    if( fChain->GetBranch("mc_truth_tWtau2_E") ) fChain->SetBranchAddress("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, &b_mc_truth_tWtau2_E);
    if( fChain->GetBranch("mc_truth_tWtaul1_E") ) fChain->SetBranchAddress("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, &b_mc_truth_tWtaul1_E);
    if( fChain->GetBranch("mc_truth_tWtaul2_E") ) fChain->SetBranchAddress("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, &b_mc_truth_tWtaul2_E);
    if( fChain->GetBranch("mc_truth_tWq11_E") ) fChain->SetBranchAddress("mc_truth_tWq11_E", &mc_truth_tWq11_E, &b_mc_truth_tWq11_E);
    if( fChain->GetBranch("mc_truth_tWq21_E") ) fChain->SetBranchAddress("mc_truth_tWq21_E", &mc_truth_tWq21_E, &b_mc_truth_tWq21_E);
    if( fChain->GetBranch("mc_truth_tWq12_E") ) fChain->SetBranchAddress("mc_truth_tWq12_E", &mc_truth_tWq12_E, &b_mc_truth_tWq12_E);
    if( fChain->GetBranch("mc_truth_tWq22_E") ) fChain->SetBranchAddress("mc_truth_tWq22_E", &mc_truth_tWq22_E, &b_mc_truth_tWq22_E);   

    if( fChain->GetBranch("mc_truth_t_E") ) fChain->SetBranchAddress("mc_truth_t_E", &mc_truth_t_E, &b_mc_truth_t_E);
    if( fChain->GetBranch("mc_truth_tb_E") ) fChain->SetBranchAddress("mc_truth_tb_E", &mc_truth_tb_E, &b_mc_truth_tb_E);
    if( fChain->GetBranch("mc_truth_tW_E") ) fChain->SetBranchAddress("mc_truth_tW_E", &mc_truth_tW_E, &b_mc_truth_tW_E);
    if( fChain->GetBranch("mc_truth_tWl_E") ) fChain->SetBranchAddress("mc_truth_tWl_E", &mc_truth_tWl_E, &b_mc_truth_tWl_E);
    if( fChain->GetBranch("mc_truth_tWq1_E") ) fChain->SetBranchAddress("mc_truth_tWq1_E", &mc_truth_tWq1_E, &b_mc_truth_tWq1_E);
    if( fChain->GetBranch("mc_truth_tWq2_E") ) fChain->SetBranchAddress("mc_truth_tWq2_E", &mc_truth_tWq2_E, &b_mc_truth_tWq2_E);   

    if( fChain->GetBranch("mc_truth_W_E") ) fChain->SetBranchAddress("mc_truth_W_E", &mc_truth_W_E, &b_mc_truth_W_E);
    if( fChain->GetBranch("mc_truth_Wl_E") ) fChain->SetBranchAddress("mc_truth_Wl_E", &mc_truth_Wl_E, &b_mc_truth_Wl_E);
    if( fChain->GetBranch("mc_truth_Wtau_E") ) fChain->SetBranchAddress("mc_truth_Wtau_E", &mc_truth_Wtau_E, &b_mc_truth_Wtau_E);
    if( fChain->GetBranch("mc_truth_Wtaul_E") ) fChain->SetBranchAddress("mc_truth_Wtaul_E", &mc_truth_Wtaul_E, &b_mc_truth_Wtaul_E);
    if( fChain->GetBranch("mc_truth_Wq1_E") ) fChain->SetBranchAddress("mc_truth_Wq1_E", &mc_truth_Wq1_E, &b_mc_truth_Wq1_E);
    if( fChain->GetBranch("mc_truth_Wq2_E") ) fChain->SetBranchAddress("mc_truth_Wq2_E", &mc_truth_Wq2_E, &b_mc_truth_Wq2_E);
    if( fChain->GetBranch("mc_truth_Z_E") ) fChain->SetBranchAddress("mc_truth_Z_E", &mc_truth_Z_E, &b_mc_truth_Z_E);
    if( fChain->GetBranch("mc_truth_Zl1_E") ) fChain->SetBranchAddress("mc_truth_Zl1_E", &mc_truth_Zl1_E, &b_mc_truth_Zl1_E);
    if( fChain->GetBranch("mc_truth_Zl2_E") ) fChain->SetBranchAddress("mc_truth_Zl2_E", &mc_truth_Zl2_E, &b_mc_truth_Zl2_E);
    if( fChain->GetBranch("mc_truth_Ztau1_E") ) fChain->SetBranchAddress("mc_truth_Ztau1_E", &mc_truth_Ztau1_E, &b_mc_truth_Ztau1_E);
    if( fChain->GetBranch("mc_truth_Ztau2_E") ) fChain->SetBranchAddress("mc_truth_Ztau2_E", &mc_truth_Ztau2_E, &b_mc_truth_Ztau2_E);
    if( fChain->GetBranch("mc_truth_Ztaul1_E") ) fChain->SetBranchAddress("mc_truth_Ztaul1_E", &mc_truth_Ztaul1_E, &b_mc_truth_Ztaul1_E);
    if( fChain->GetBranch("mc_truth_Ztaul2_E") ) fChain->SetBranchAddress("mc_truth_Ztaul2_E", &mc_truth_Ztaul2_E, &b_mc_truth_Ztaul2_E);    
    if( fChain->GetBranch("mc_truth_Zq1_E") ) fChain->SetBranchAddress("mc_truth_Zq1_E", &mc_truth_Zq1_E, &b_mc_truth_Zq1_E);
    if( fChain->GetBranch("mc_truth_Zq2_E") ) fChain->SetBranchAddress("mc_truth_Zq2_E", &mc_truth_Zq2_E, &b_mc_truth_Zq2_E);
    if( fChain->GetBranch("mc_truth_gammal1_E") ) fChain->SetBranchAddress("mc_truth_gammal1_E", &mc_truth_gammal1_E, &b_mc_truth_gammal1_E);
    if( fChain->GetBranch("mc_truth_gammal2_E") ) fChain->SetBranchAddress("mc_truth_gammal2_E", &mc_truth_gammal2_E, &b_mc_truth_gammal2_E);
    if( fChain->GetBranch("mc_truth_gammatau1_E") ) fChain->SetBranchAddress("mc_truth_gammatau1_E", &mc_truth_gammatau1_E, &b_mc_truth_gammatau1_E);
    if( fChain->GetBranch("mc_truth_gammatau2_E") ) fChain->SetBranchAddress("mc_truth_gammatau2_E", &mc_truth_gammatau2_E, &b_mc_truth_gammatau2_E);
    if( fChain->GetBranch("mc_truth_gammataul1_E") ) fChain->SetBranchAddress("mc_truth_gammataul1_E", &mc_truth_gammataul1_E, &b_mc_truth_gammataul1_E);
    if( fChain->GetBranch("mc_truth_gammataul2_E") ) fChain->SetBranchAddress("mc_truth_gammataul2_E", &mc_truth_gammataul2_E, &b_mc_truth_gammataul2_E);
    if( fChain->GetBranch("mc_truth_gamma_E") ) fChain->SetBranchAddress("mc_truth_gamma_E", &mc_truth_gamma_E, &b_mc_truth_gamma_E);
}
