//ThreeLeptonSelection_THQ3l

#include "../include/TTbarHiggsMultileptonAnalysis.h"
#include "TSystem.h"
#include "SignalExtractionMVA.cxx"
#include "Helper.cxx"
#include "BTagging.cxx"
#include "FakeRate.cxx"
#include "ChargeFlip.cxx"

#define kCat_3l_2b_2j   0
#define kCat_3l_1b_2j   1
#define kCat_3l_2b_1j   2
#define kCat_3l_1b_1j   3
#define kCat_3l_2b_0j   4

#define DEBUG           false

using namespace std;

//Ccutflow vector -- new
vector<double> v_cutflow(15);


TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis()
{
    cout<<BOLD(FRED("DO NOT USE THE DEFAULT CONSTRUCTOR !"))<<endl<<endl;
}

TTbarHiggsMultileptonAnalysis::~TTbarHiggsMultileptonAnalysis()
{
    delete vEvent; delete vElectron; delete vMuon; delete vTau; delete vJet; delete vTruth;

    delete tOutput;
    delete outputfile;

    delete f_CSVwgt_HF; delete f_CSVwgt_LF; delete f_BTag_eff; delete f_QFwgt; delete f_FRwgt;
}

TTbarHiggsMultileptonAnalysis::TTbarHiggsMultileptonAnalysis(TString inputFileName, TString sampleName, TString treeName, TString outputFileName, bool isdata, bool doSystCombine, float xsec, float lumi, int nowe, int nmax)
{
    //
    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;
    _process = "toto";


    fChain = new TChain(treeName);

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        fChain->Add(fnameStr.c_str());
    }
    infile.close();
    InitTree();

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot.Data(), "RECREATE");

    Init();
    initializeOutputTree();
}



void TTbarHiggsMultileptonAnalysis::InitTree()
{
    // Set branch addresses and branch pointers
    if (!fChain) return;

    vEvent    = new std::vector<Event>();
    vElectron = new std::vector<Electron>();
    vMuon     = new std::vector<Muon>();
    vTau      = new std::vector<Tau>();
    vJet      = new std::vector<Jet>();
    vTruth    = new std::vector<Truth>();

    fChain->SetBranchAddress("Event",    &vEvent   );
    fChain->SetBranchAddress("Electron", &vElectron);
    fChain->SetBranchAddress("Muon",     &vMuon    );
    fChain->SetBranchAddress("Tau",      &vTau     );
    fChain->SetBranchAddress("Jet",      &vJet     );
    fChain->SetBranchAddress("Truth",    &vTruth   );
}

//NEW
void TTbarHiggsMultileptonAnalysis::Init()
{
    //--- Load MVA weight files
    Load_MVA();

    // Loading weight files and creating corrsponding histograms

    // b-tagging
    std::string inputFileHF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_hf_76x_2016_02_08.root"; //FIXME -- needs update ??
    std::string inputFileLF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_lf_76x_2016_02_08.root";
    f_CSVwgt_HF = new TFile ((inputFileHF).c_str());
    f_CSVwgt_LF = new TFile ((inputFileLF).c_str());
    fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);

    // new b-tagging (using BTagCalibrationXStandaloneWhatever)
    // setup calibration + reader
    // BTagCalibrationX *
    calib = BTagCalibrationX("csvv2", "/opt/sbg/scratch1/cms/TTH/weight/CSVv2_Moriond17_B_H.csv");

    // BTagCalibrationXReader *
    reader = BTagCalibrationXReader(  BTagEntryX::OP_LOOSE,  // operating point
                                                            "central",            // central sys type
                                                            {"up", "down"});      // other sys types

    reader.load(   calib,                // calibration instance
                    BTagEntryX::FLAV_B,    // btag flavour
                    "comb");               // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_C,    // btag flavour
                    "comb");              // measurement type

    reader.load(    calib,                // calibration instance
                    BTagEntryX::FLAV_UDSG, // btag flavour
                    "comb");              // measurement type


    // BTag Efficiencies
    std::string inputFileBTagEff = "/opt/sbg/scratch1/cms/TTH/weight/output_tt_effBtag.root";
    f_BTag_eff = new TFile ((inputFileBTagEff).c_str());
    fill_eff_btagging_histos(f_BTag_eff);

    // charge flip
    std::string inputFileQF = "/opt/sbg/scratch1/cms/TTH/weightMoriond2017/FakeFlip/QF_data_el.root";
    //std::string inputFileQF = "/opt/sbg/scratch1/cms/TTH/weight/QF_data_el.root";
    f_QFwgt    = new TFile ((inputFileQF).c_str());
    fillQFhistos(f_QFwgt);

    // fake rate
    std::string inputFileFR = "/opt/sbg/scratch1/cms/TTH/weightMoriond2017/FakeFlip/FR_data_ttH_mva.root";
    //std::string inputFileFR = "/opt/sbg/scratch1/cms/TTH/weight/FR_data_ttH_mva.root";
    f_FRwgt    = new TFile ((inputFileFR).c_str());
    fillFRhistos(f_FRwgt);
}

void TTbarHiggsMultileptonAnalysis::InitEvent()
{
    vLeptons.clear();
    vSelectedMuons.clear();
    vSelectedElectrons.clear();
    vSelectedTaus.clear();
    vSelectedLeptons.clear();
    vFakeMuons.clear();
    vFakeElectrons.clear();
    vFakeLeptons.clear();
    vFakeableLeptons.clear();
    vLooseLeptons.clear();
    vTightLeptons.clear();

    vSelectedJets.clear();
    vLooseBTagJets.clear();
    vLightJets.clear();

    weightfake = 0;
    weightflip = 0;

    is_emu_TT_CR        = false;

    is_3l_THQ_SR        = 0;
    is_3l_THQ_Training  = 0;
    is_3l_AppFakes_SR   = 0;

    signal_3l_TT_MVA = 0; signal_3l_TTV_MVA = 0;

    cat_HtoWW       = 0;   cat_HtoZZ       = 0;   cat_Htott     = 0;

    n_tight       = 0;

    return;
}


void TTbarHiggsMultileptonAnalysis::Loop()
{
    if (!fChain) return;

    for(int i=0; i<15; i++)
    {
        v_cutflow[i] = 0;
    }

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Will process " << nentries_max << " events" << std::endl << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++)
    {
        // cout<<endl<<jentry<< " / "<<nentries_max<<endl;

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;

        if(DEBUG) std::cout << "New event ===" << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //
        int pvn = vEvent->at(0).pv_n();

        mc_event = vEvent->at(0).id();

        weights_pdf.clear();
        ids_pdf.clear();

        if(!_isdata )
        {
            weight = _lumi*_xsec/_nowe;

            mc_weight = vEvent->at(0).mc_weight();
            //weight_PV = _h_PV->GetBinContent(pvn);
            weight = weight * mc_weight; //*weight_PV;

            weight_scale_muF0p5 = vEvent->at(0).weight_scale_muF0p5();
            weight_scale_muF2 = vEvent->at(0).weight_scale_muF2();
            weight_scale_muR0p5 = vEvent->at(0).weight_scale_muR0p5();
            weight_scale_muR2 = vEvent->at(0).weight_scale_muR2();

            weights_pdf = vEvent->at(0).pdf_weights();
            ids_pdf = vEvent->at(0).pdf_ids();
        }
        else
        {
            weight    = 1.;
            mc_weight = 1.;
            weight_PV = 1.;
        }


        // ###########################################################
        // #  _       _ _   _       _ _           _   _              #
        // # (_)_ __ (_) |_(_) __ _| (_)___  __ _| |_(_) ___  _ __   #
        // # | | '_ \| | __| |/ _` | | / __|/ _` | __| |/ _ \| '_ \  #
        // # | | | | | | |_| | (_| | | \__ \ (_| | |_| | (_) | | | | #
        // # |_|_| |_|_|\__|_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_| #
        // #                                                         #
        // ###########################################################

        InitEvent();


        // ######################################
        // #  _        _                        #
        // # | |_ _ __(_) __ _  __ _  ___ _ __  #
        // # | __| '__| |/ _` |/ _` |/ _ \ '__| #
        // # | |_| |  | | (_| | (_| |  __/ |    #
        // #  \__|_|  |_|\__, |\__, |\___|_|    #
        // #             |___/ |___/            #
        // #                                    #
        // ######################################

        if(_isdata)
        {
            //
            bool TRIGm   = vEvent->at(0).is_TRIGm()  ;
            bool TRIGe   = vEvent->at(0).is_TRIGe()  ;
            bool TRIGmTk = vEvent->at(0).is_TRIGmTk();
            bool TRIGee  = vEvent->at(0).is_TRIGee() ;
            bool TRIGmm  = vEvent->at(0).is_TRIGmm() ;
            bool TRIGme  = vEvent->at(0).is_TRIGme() ;
            bool TRIGem  = vEvent->at(0).is_TRIGem() ;
            bool TRIGmmTk= vEvent->at(0).is_TRIGmmTk();
            bool TRIGeee = vEvent->at(0).is_TRIGeee();
            bool TRIGmme = vEvent->at(0).is_TRIGmme();
            bool TRIGeem = vEvent->at(0).is_TRIGeem();
            bool TRIGmmm = vEvent->at(0).is_TRIGmmm();

            bool E = false, M = false, EE = false, MM = false, EM = false;
            if ( TRIGme || TRIGem || TRIGeem || TRIGmme ) EM = true;
            if ( TRIGmm || TRIGmmTk || TRIGmmm )          MM = true;
            if ( TRIGee || TRIGeee )	                  EE = true;
            if ( TRIGm  || TRIGmTk )                      M  = true;
            if ( TRIGe  )                                 E  = true;

            bool emdataset = _sampleName.Contains("MuonE");
            bool mmdataset = _sampleName.Contains("DoubleM");
            bool eedataset = _sampleName.Contains("DoubleE");
            bool mdataset  = _sampleName.Contains("SingleM");
            bool edataset  = _sampleName.Contains("SingleE");

            is_trigger = false;
            if ( EM  &&                               (emdataset) ) is_trigger = true;
            if ( !EM && MM  &&                        (mmdataset) ) is_trigger = true;
            if ( !EM && !MM && EE  &&                 (eedataset) ) is_trigger = true;
            if ( !EM && !MM && !EE && M  &&           (mdataset ) ) is_trigger = true;
            if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) is_trigger = true;

			// cout<<"eedataset : "<<eedataset<<endl;
			// cout<<" / EE : "<<EE;
			// cout<<" / EM : "<<EM;
			// cout<<" / MM : "<<MM;


            if(is_trigger)
            {
                weight *= 1;
            }
            else
            {
                weight *= 0;
            }
        }

        // #####################################
        // #  _ __ ___  _   _  ___  _ __  ___  #
        // # | '_ ` _ \| | | |/ _ \| '_ \/ __| #
        // # | | | | | | |_| | (_) | | | \__ \ #
        // # |_| |_| |_|\__,_|\___/|_| |_|___/ #
        // #                                   #
        // #####################################

        for(unsigned int imuon=0; imuon < vMuon->size() ; imuon++)
        {
            Lepton l; l.setLepton(&vMuon->at(imuon),imuon,0,1);

            vLeptons.push_back(l);

            if(vMuon->at(imuon).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }
            if(vMuon->at(imuon).isLooseTTH() )
            {
                vLooseLeptons.push_back(l);
            }
            if(vMuon->at(imuon).isTightTTH() )
            {
                vTightLeptons.push_back(l);
                n_tight++;
            }
        }


        // ##############################################
        // #       _           _                        #
        // #   ___| | ___  ___| |_ _ __ ___  _ __  ___  #
        // #  / _ \ |/ _ \/ __| __| '__/ _ \| '_ \/ __| #
        // # |  __/ |  __/ (__| |_| | | (_) | | | \__ \ #
        // #  \___|_|\___|\___|\__|_|  \___/|_| |_|___/ #
        // #                                            #
        // ##############################################

        for(unsigned int ielectron=0; ielectron < vElectron->size() ; ielectron++)
        {
            Lepton l; l.setLepton(&vElectron->at(ielectron),ielectron,1,0);

            vLeptons.push_back(l);

            if(vElectron->at(ielectron).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }
            if(vElectron->at(ielectron).isLooseTTH() )
            {
                vLooseLeptons.push_back(l);
            }
            if(vElectron->at(ielectron).isTightTTH() )
            {
                if(!vElectron->at(ielectron).cutEventSel() || !vElectron->at(ielectron).noLostHits() ) {continue;} //conversion veto
                vTightLeptons.push_back(l);
                n_tight++;
            }
        }

        // ########################
        // #  _                   #
        // # | |_ __ _ _   _ ___  #
        // # | __/ _` | | | / __| #
        // # | || (_| | |_| \__ \ #
        // #  \__\__,_|\__,_|___/ #
        // #                      #
        // ########################

        for(unsigned int itau=0; itau < vTau->size() ; itau++)
        {
            Lepton l; l.setLepton(&vTau->at(itau),itau,0,0);

            if( vTau->at(itau).isTightTTH() )
            {
                vSelectedTaus.push_back(vTau->at(itau));
            }
        }

        // #############################################
        // #                _           _              #
        // #   ___  _ __ __| | ___ _ __(_)_ __   __ _  #
        // #  / _ \| '__/ _` |/ _ \ '__| | '_ \ / _` | #
        // # | (_) | | | (_| |  __/ |  | | | | | (_| | #
        // #  \___/|_|  \__,_|\___|_|  |_|_| |_|\__, | #
        // #                                    |___/  #
        // #                                           #
        // #############################################

        std::sort(vLeptons.begin(), vLeptons.end(), SortingLeptonPt);
        std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);
        std::sort(vFakeableLeptons.begin(), vFakeableLeptons.end(), SortingLeptonPt);
        std::sort(vLooseLeptons.begin(), vLooseLeptons.end(), SortingLeptonPt);
        std::sort(vTightLeptons.begin(), vTightLeptons.end(), SortingLeptonPt);


        // ################################
        // #                              #
        // #  _           _      _        #
        // # | |__       (_) ___| |_ ___  #
        // # | '_ \ _____| |/ _ \ __/ __| #
        // # | |_) |_____| |  __/ |_\__ \ #
        // # |_.__/     _/ |\___|\__|___/ #
        // #           |__/               #
        // #                              #
        // ################################

        nLooseBJets  = 0;
        nMediumBJets = 0;
        nForwardJets = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {
            if(vJet->at(ijet).pt() < 25) {continue;}

            //-- should cleanup jet within R<0.4 of fakeable leptons here ?

            //NB : should also consider jets with CSV>0.54 & eta > 2.4 ? (transition region)

            //--- B-tagged jets
            // if(vJet->at(ijet).CSVv2() > 0.5426 && fabs(vJet->at(ijet).eta()) < 2.4) //CHANGED
            if(vJet->at(ijet).CSVv2() > 0.5426)
            {
                vLooseBTagJets.push_back(vJet->at(ijet));
                nLooseBJets++;

                if(vJet->at(ijet).CSVv2() > 0.8484) {nMediumBJets++;}
            }
            //--- Light jets
            else if(vJet->at(ijet).CSVv2() < 0.5426)
            {
                vLightJets.push_back(vJet->at(ijet));
                nForwardJets++;

                //FIXME : cf. AN, should apply this cut. But then yield < yield from AN... !! Coincidence ? ASK
                // if(fabs(vJet->at(ijet).eta() ) < 2.4)
                // {
                //     nForwardJets++;
                // }
                // else if(vJet->at(ijet).pt() > 40 && fabs(vJet->at(ijet).eta() ) > 2.4)
                // {
                //     nForwardJets++;
                // }
            }

            //---------
			//All jets with pT>25 in vJets are "selected jets"
            vSelectedJets.push_back(vJet->at(ijet));
        }


        // ###################################
        // #  _____ ____  _   _ _____ _   _  #
        // # |_   _|  _ \| | | |_   _| | | | #
        // #   | | | |_) | | | | | | | |_| | #
        // #   | | |  _ <| |_| | | | |  _  | #
        // #   |_| |_| \_\\___/  |_| |_| |_| #
        // #                                 #
        // ###################################

        if( !_isdata )
        {
            for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
            {

                //std::cout << vTruth->at(0).mc_truth_label().at(itruth) << std::endl;
                if( vTruth->at(0).mc_truth_label().at(itruth) == 12 )
                {   cat_HtoWW = true;
                    //std::cout << "H to WW" << std::endl;
                }
                else if( vTruth->at(0).mc_truth_label().at(itruth) == 14 )
                {   cat_HtoZZ = true;
                    //std::cout << "H to ZZ" << std::endl;
                }
                else if( vTruth->at(0).mc_truth_label().at(itruth) == 16 )
                {   cat_Htott = true;
                    //std::cout << "H to tau tau" << std::endl;
                }
            }
        }

        // ################################################################################
        // #  ____  ____    ____  ____ _____                   _       _     _            #
        // # |___ \|  _ \  | __ )|  _ \_   _| __   ____ _ _ __(_) __ _| |__ | | ___  ___  #
        // #   __) | | | | |  _ \| | | || |   \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __| #
        // #  / __/| |_| | | |_) | |_| || |    \ V / (_| | |  | | (_| | |_) | |  __/\__ \ #
        // # |_____|____/  |____/|____/ |_|     \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/ #
        // #                                                                              #
        // ################################################################################

        //-- Old ttH variables
        max_Lep_eta     = 0. ;
        nJet25_Recl     = 0  ;
        mindr_lep1_jet  = 0. ;
        mindr_lep2_jet  = 0. ;
        met             = 0. ;
        avg_dr_jet      = 0. ;
        MT_met_lep1     = 0. ;
        LepGood_conePt0 = 0. ;
        LepGood_conePt1 = 0. ;
        mhtJet25_Recl   = 0. ;

        //--- NEW -- tHq2016 vars
        nJet25 = 0.;
        maxEtaJet25 = 0.;
        lepCharge = 0.;
        nJetEta1  = 0.;
        dEtaFwdJetBJet = 0.;
        dEtaFwdJet2BJet = 0.;
        dEtaFwdJetClosestLep = 0.;
        dPhiHighestPtSSPair = 0.;
        minDRll = 0.;
        Lep3Pt = 0.;


        // ############################################
        // #           _           _   _              #
        // #  ___  ___| | ___  ___| |_(_) ___  _ __   #
        // # / __|/ _ \ |/ _ \/ __| __| |/ _ \| '_ \  #
        // # \__ \  __/ |  __/ (__| |_| | (_) | | | | #
        // # |___/\___|_|\___|\___|\__|_|\___/|_| |_| #
        // #                                          #
        // ############################################

        ThreeLeptonSelection_THQ3l_SR(jentry);
        ThreeLeptonSelection_THQ3l_Training(jentry);

        //ThreeLeptonSelection_TTH3l(jentry);
        //ThreeLeptonSelection_ApplicationFakes(jentry);
        //ThreeLeptonSelection_CR_WZ(jentry);
        //ThreeLeptonSelection_CR_WZ_ApplicationFakes(jentry);
        //ThreeLeptonSelection_CR_WZrelaxed(jentry);
        //ThreeLeptonSelection_TTZ(jentry);
        //ThreeLeptonSelection_CR_Zl(jentry);
    }

	ofstream file_out("cutflow.txt");
    for(int i=0; i<v_cutflow.size(); i++)
    {
		// cout<<"v_cutflow "<<i<<" = "<<v_cutflow[i]<<endl;
		file_out<<"v_cutflow "<<i<<" = "<<v_cutflow[i]<<endl;
    }

	outputfile->cd();
    tOutput->Write();
}

// ##############################################################################
// #        _ _       _                   _                  _                  #
// #   __ _| | |  ___(_) __ _ _ __   __ _| |  _ __ ___  __ _(_) ___  _ __  ___  #
// #  / _` | | | / __| |/ _` | '_ \ / _` | | | '__/ _ \/ _` | |/ _ \| '_ \/ __| #
// # | (_| | | | \__ \ | (_| | | | | (_| | | | | |  __/ (_| | | (_) | | | \__ \ #
// #  \__,_|_|_| |___/_|\__, |_| |_|\__,_|_| |_|  \___|\__, |_|\___/|_| |_|___/ #
// #                    |___/                          |___/                    #
// #                                                                            #
// ##############################################################################






/*
████████ ██   ██  ██████  ██████  ██
   ██    ██   ██ ██    ██      ██ ██
   ██    ███████ ██    ██  █████  ██
   ██    ██   ██ ██ ▄▄ ██      ██ ██
   ██    ██   ██  ██████  ██████  ███████
                     ▀▀
*/



/*

//FIXME -- TESTING BACK PREVIOUS VERSION OF FUNCTION --> FS COUNT CORRECT ???

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_SR(int evt)
{
	v_cutflow[0]++;

    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)
    //if(vSelectedTaus.size()>0) return;

    v_cutflow[1]++;


    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss SR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    v_cutflow[2]++;


    if(DEBUG) std::cout << "nLep Ok... ";

    bool nTight             = ( n_tight                     == 3 ); //CHANGED
    if(!nTight)             return;

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              //&& vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge() && vSelectedLeptons.at(2).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    // if(!are_fullytight)     return; //CHANGED
    if(DEBUG) std::cout << "nTight Ok... ";


    v_cutflow[3]++;


    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 25 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 15) && (vSelectedLeptons.at(2).pt() > 15) );
    if(!following_lep_pt)   return;


    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    v_cutflow[4]++;


    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nMediumBtag)      return; //CHANGED

    v_cutflow[5]++;

    bool nForward        = ( nForwardJets                >= 1 );
    if(!nForward)      return; //CHANGED

    v_cutflow[6]++;


    if(DEBUG) std::cout << "Btag Ok... ";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;

    v_cutflow[7]++;

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id()     == -vSelectedLeptons.at(j).id()                             )
               && ( vSelectedLeptons.at(i).charge() == -vSelectedLeptons.at(j).charge()                         )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 15 ) )
            { pass_Zveto = false ;} //CHANGED //15 GeV window
        }
    }
    if(!pass_Zveto)       return;

    v_cutflow[8]++;


    if(DEBUG) std::cout << "Zveto Ok... ";

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;


    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }


    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    // if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return; //CHANGED

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    v_cutflow[9]++;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    is_3l_THQ_SR = 1;

    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py

    // ======================================================================================================
    // variables against ttbar

    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    nJet25_Recl     = vSelectedJets.size() ;
    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ); }
    }

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");

    if( DEBUG ) std::cout << "Signal 3l TT MVA ok" << std::endl;

    // ======================================================================================================
    // variables against ttV

    met             = vEvent->at(0).metpt() ;

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(2).pt() ;

    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");

    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    // ======================================================================================================
    //variables from tHq2016
    //NEW //FIXME
    //NB : "forward jet" = most forward non-CSV loose jet

    // double nJet25;
    // double MaxEtaJet25;
    // double totCharge;
    // double nJetEta1 ;
    // double detaFwdJetBJet;
    // double detaFwdJet2BJet;
    // double detaFwdJetClosestLep;
    // double dphiHighestPtSPPair;
    // double minDRll;
    // double Lep3Pt;

    nJet25 = 0.;
    maxEtaJet25 = 0.;
    lepCharge = 0.;
    nJetEta1  = 0.;
    dEtaFwdJetBJet = 0.;
    dEtaFwdJet2BJet = 0.;
    dEtaFwdJetClosestLep = 0.;
    dPhiHighestPtSSPair = 0.;
    minDRll = 0.;
    Lep3Pt = 0.;

    int ijet_forward=-1, ijet_hardest_btag=-1, ijet_2nd_hardest_btag=-1;
    double tmp = -999;

	//--- Find "forward jet"
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp)
        {
            tmp = fabs(vLightJets.at(ijet).eta());
            ijet_forward = ijet;
        }
    }

	//--- Find "hardest" and "second hardest" tagged jets
    double pt1=0, pt2=0;
    for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
    {
        if(fabs(vLooseBTagJets.at(ijet).pt()) > pt1)
        {
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            ijet_hardest_btag = ijet;
        }
        else if(fabs(vLooseBTagJets.at(ijet).pt()) > pt2)
        {
            pt2 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet;
        }
    }



	//--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt() > 25 && fabs(vSelectedJets.at(ijet).eta() ) < 2.4) {nJet25++;}
    }

    //--- Var2: max eta of any 'non-CSV-loose' jet
    tmp = -999;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp) {tmp = fabs(vLightJets.at(ijet).eta());}
    }
    maxEtaJet25 = tmp;

	//--- Var3 : sum of leptons charges
    lepCharge = vSelectedLeptons.at(0).charge() + vSelectedLeptons.at(1).charge() + vSelectedLeptons.at(2).charge();

	//--- Var4 : nof 'non-csv-loose' jets with eta>1.0
    nJetEta1 = 0;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if( fabs(vLightJets.at(ijet).eta() ) > 1.0) {nJetEta1++;}
    }
    nJetEta1 = tmp;

	//--- Var5 : dEta between forward light jet and hardest tagged jet
    dEtaFwdJetBJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_hardest_btag).eta() );

	//--- Var6 : dEta between forward and 2nd hardest tagged jet
    if(ijet_2nd_hardest_btag < 0) {dEtaFwdJet2BJet = -1;}
    else {dEtaFwdJet2BJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_2nd_hardest_btag).eta() );}

	//--- Var7 : dEta between forward light jet and closet lepton (angular dist.) //FIXME -- only distance in eta ?
    dEtaFwdJetClosestLep = 0; tmp = 999;
	TLorentzVector jet_tmp;
	jet_tmp.SetPtEtaPhiE(vLightJets.at(ijet_forward).pt(), vLightJets.at(ijet_forward).eta(), vLightJets.at(ijet_forward).phi(), vLightJets.at(ijet_forward).E());
    for(int ilep=0; ilep<vSelectedLeptons.size(); ilep++)
    {
		TLorentzVector lep_tmp;
		lep_tmp.SetPtEtaPhiE(vSelectedLeptons.at(ilep).pt(), vSelectedLeptons.at(ilep).eta(), vSelectedLeptons.at(ilep).phi(), vSelectedLeptons.at(ilep).E());

        if(jet_tmp.DeltaR(lep_tmp) < tmp) {dEtaFwdJetClosestLep = fabs( jet_tmp.Eta() - lep_tmp.Eta() ); tmp = jet_tmp.DeltaR(lep_tmp);}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair (?)
    if(vSelectedLeptons.at(0).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs(vSelectedLeptons.at(1).phi() - vSelectedLeptons.at(2).phi() );}
    else if(vSelectedLeptons.at(1).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs(vSelectedLeptons.at(0).phi() - vSelectedLeptons.at(2).phi() );}
    else {dPhiHighestPtSSPair = fabs(vSelectedLeptons.at(0).phi() - vSelectedLeptons.at(1).phi() );}

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999; tmp = 999;
	TLorentzVector lep1, lep2, lep3;
	lep1.SetPtEtaPhiE(vSelectedLeptons.at(0).pt(), vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi(), vSelectedLeptons.at(0).E() );
	lep2.SetPtEtaPhiE(vSelectedLeptons.at(1).pt(), vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi(), vSelectedLeptons.at(1).E() );
	lep3.SetPtEtaPhiE(vSelectedLeptons.at(2).pt(), vSelectedLeptons.at(2).eta(), vSelectedLeptons.at(2).phi(), vSelectedLeptons.at(2).E() );
	if(lep1.DeltaR(lep2) < tmp) {minDRll = lep1.DeltaR(lep2); tmp = minDRll;}
	if(lep1.DeltaR(lep3) < tmp) {minDRll = lep1.DeltaR(lep3); tmp = minDRll;}
	if(lep2.DeltaR(lep3) < tmp) {minDRll = lep2.DeltaR(lep3); tmp = minDRll;}


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vSelectedLeptons.at(2).pt();






    // ======================================================================================================

    if( DEBUG ) std::cout << "Signal 3l TTV MVA ok" << std::endl;

    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;

    // if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}



*/

















//NEW -------------------------------------------
void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_SR(int evt)
{
//---------------------------------------------------------------------------
	v_cutflow[0]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}
    v_cutflow[1]++;

    // ####################
    // # Common selection #
    // ####################

    //At least 3 leptons
    if(vFakeableLeptons.size() < 3) {return;}
    v_cutflow[2]++;

    //==3 tight leptons, pT cuts
	if(vTightLeptons.size() != 3) {return;}
    v_cutflow[3]++;

    if(vTightLeptons.at(0).pt() < 25 || vTightLeptons.at(1).pt() < 15 || vTightLeptons.at(2).pt() < 15) {return;}

    //cf. AN : electrons must have 0 lost hit & pass Conversion veto (always true for muons)
    // bool are_fullytight = true;
    // if(!vTightLeptons.at(0).isTightTTH() || !vTightLeptons.at(0).cutEventSel() || !vTightLeptons.at(0).noLostHits() ) {are_fullytight = false;}
    // if(!vTightLeptons.at(1).isTightTTH() || !vTightLeptons.at(1).cutEventSel() || !vTightLeptons.at(1).noLostHits() ) {are_fullytight = false;}
    // if(!vTightLeptons.at(2).isTightTTH() || !vTightLeptons.at(2).cutEventSel() || !vTightLeptons.at(2).noLostHits() ) {are_fullytight = false;}
    // if(!are_fullytight) {return;} //done at selection level

    if(vSelectedJets.size() < 2) {return;}
    v_cutflow[4]++;

    if(nMediumBJets == 0) {return;}
    v_cutflow[5]++;

    if(nForwardJets == 0)  {return;}
    v_cutflow[6]++;

    //Cleanup -- Veto if loose lepton pair with mll < 12
    bool pass_cleanup = true;
    for(int i=0; i<vLooseLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLooseLeptons.size(); j++)
        {
            if( fabs( (vLooseLeptons.at(i).p4() + vLooseLeptons.at(j).p4()).M() ) < 12 ) {pass_cleanup = false;}
        }
    }
    if(!pass_cleanup) {return;}
    v_cutflow[7]++;


    // ##########
    // # OSSF pair - Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vTightLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vTightLeptons.size(); j++)
        {
            if ( fabs(vTightLeptons.at(i).id()) == fabs(vTightLeptons.at(j).id()) && vTightLeptons.at(i).charge() == -vTightLeptons.at(j).charge() && (fabs( ( vTightLeptons.at(i).p4() + vTightLeptons.at(j).p4()).M() - 91.188) < 15 ) ) {pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;
    v_cutflow[8]++;



//---------------------------------------------------------------------------


    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vTightLeptons.size(); i++)
    {
        lepton_px = lepton_px + vTightLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vTightLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    if(DEBUG) std::cout << " MHT =  " << MHT << "MET = " << vEvent->at(0).metpt() << " met_ld = " << met_ld ;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vTightLeptons.size(); i++)
    {
        for(int j=0; j<vTightLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vTightLeptons.at(i).id() == -vTightLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    // if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return; //CHANGED

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vTightLeptons.at(i).charge();
    }
    // if( fabs(sum_charges_3l) != 1 ) return; //Cf. AN : triple charge consistency not required in 3l ?


    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;
    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    is_3l_THQ_SR = 1;





    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py


    //variables from tHq2016
    //NEW //FIXME
    //NB : "forward jet" = most forward non-CSV loose jet

    // double nJet25;
    // double maxEtaJet25;
    // double lepCharge;
    // double nJetEta1 ;
    // double dEtaFwdJetBJet;
    // double dEtaFwdJet2BJet;
    // double dEtaFwdJetClosestLep;
    // double dPhiHighestPtSSPair;
    // double minDRll;
    // double Lep3Pt;

    int ijet_forward=-1, ijet_hardest_btag=-1, ijet_2nd_hardest_btag=-1;
    double tmp = -999;

	//--- Find "forward jet"
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp)
        {
            tmp = fabs(vLightJets.at(ijet).eta());
            ijet_forward = ijet;
        }
    }

	//--- Find "hardest" and "second hardest" tagged jets
    double pt1=0, pt2=0;
    for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
    {
        if(fabs(vLooseBTagJets.at(ijet).pt()) > pt1)
        {
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            ijet_hardest_btag = ijet;
        }
        else if(fabs(vLooseBTagJets.at(ijet).pt()) > pt2)
        {
            pt2 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet;
        }
    }



	//--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt() > 25 && fabs(vSelectedJets.at(ijet).eta() ) < 2.4) {nJet25++;}
    }

    //--- Var2: max eta of any 'non-CSV-loose' jet
    tmp = -999;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp) {tmp = fabs(vLightJets.at(ijet).eta());}
    }
    maxEtaJet25 = tmp;

	//--- Var3 : sum of leptons charges
    lepCharge = vTightLeptons.at(0).charge() + vTightLeptons.at(1).charge() + vTightLeptons.at(2).charge();

	//--- Var4 : nof 'non-csv-loose' jets with eta>1.0
    nJetEta1 = 0;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if( fabs(vLightJets.at(ijet).eta() ) > 1.0) {nJetEta1++;}
    }

	//--- Var5 : dEta between forward light jet and hardest tagged jet
    dEtaFwdJetBJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_hardest_btag).eta() );

	//--- Var6 : dEta between forward and 2nd hardest tagged jet
    if(ijet_2nd_hardest_btag < 0) {dEtaFwdJet2BJet = -1;}
    else {dEtaFwdJet2BJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_2nd_hardest_btag).eta() );}

	//--- Var7 : dEta between forward light jet and closet lepton (angular dist.)
    dEtaFwdJetClosestLep = 0; tmp = 999;
    for(int ilep=0; ilep<vTightLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vTightLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vTightLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair (?)
    if(vTightLeptons.at(0).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vTightLeptons.at(1).phi() - vTightLeptons.at(2).phi() ) );}
    else if(vTightLeptons.at(1).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs( Phi_MPi_Pi(vTightLeptons.at(0).phi() - vTightLeptons.at(2).phi() ) );}
    else {dPhiHighestPtSSPair = fabs( Phi_MPi_Pi(vTightLeptons.at(0).phi() - vTightLeptons.at(1).phi() ) );}

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999; tmp = 999;
	TLorentzVector lep1, lep2, lep3;
	lep1.SetPtEtaPhiE(vTightLeptons.at(0).pt(), vTightLeptons.at(0).eta(), vTightLeptons.at(0).phi(), vTightLeptons.at(0).E() );
	lep2.SetPtEtaPhiE(vTightLeptons.at(1).pt(), vTightLeptons.at(1).eta(), vTightLeptons.at(1).phi(), vTightLeptons.at(1).E() );
	lep3.SetPtEtaPhiE(vTightLeptons.at(2).pt(), vTightLeptons.at(2).eta(), vTightLeptons.at(2).phi(), vTightLeptons.at(2).E() );
	if(lep1.DeltaR(lep2) < tmp) {minDRll = lep1.DeltaR(lep2); tmp = minDRll;}
	if(lep1.DeltaR(lep3) < tmp) {minDRll = lep1.DeltaR(lep3); tmp = minDRll;}
	if(lep2.DeltaR(lep3) < tmp) {minDRll = lep2.DeltaR(lep3); tmp = minDRll;}


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vTightLeptons.at(2).pt();




    //--------------------------
    // PRINTOUT OF INFOS & INPUT VARS
    //--------------------------

    bool do_printout = false;

    if(do_printout)
    {
        cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

        cout<<FYEL("--- Tagged Jets : ")<<endl;
        for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt()<<" / eta = "<<vLooseBTagJets.at(ijet).eta()<<" / phi = "<<vLooseBTagJets.at(ijet).phi()<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2()<<endl;
        }
        cout<<"Hardest & 2nd hardest jets are "<<ijet_hardest_btag<<", "<<ijet_2nd_hardest_btag<<endl<<endl;

        cout<<FYEL("--- Forward Jets : ")<<endl;
        for(int ijet=0; ijet<vLightJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt()<<" / eta = "<<vLightJets.at(ijet).eta()<<" / phi = "<<vLightJets.at(ijet).phi()<<" / CSV = "<<vLightJets.at(ijet).CSVv2()<<endl;
        }
        cout<<"Forwardest jet is "<<ijet_forward<<endl<<endl;

        cout<<FYEL("--- Selected Leptons : ")<<endl;
        for(int ilep=0; ilep<vTightLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vTightLeptons.at(ilep).pt()<<" / eta = "<<vTightLeptons.at(ilep).eta()<<" / phi = "<<vTightLeptons.at(ilep).phi()<<" : Charge = "<<vTightLeptons.at(ilep).charge()<<endl;
        }

        cout<<FYEL("--- Input variables : ")<<endl;

        cout<<"nJet25 = "<<nJet25<<endl;
        cout<<"maxEtaJet25 = "<<maxEtaJet25<<endl;
        cout<<"lepCharge = "<<lepCharge<<endl;
        cout<<"nJetEta1 = "<<nJetEta1 <<endl;
        cout<<"dEtaFwdJetBJet = "<<dEtaFwdJetBJet<<endl;
        cout<<"dEtaFwdJet2BJet = "<<dEtaFwdJet2BJet<<endl;
        cout<<"dEtaFwdJetClosestLep = "<<dEtaFwdJetClosestLep<<endl;
        cout<<"dPhiHighestPtSSPair = "<<dPhiHighestPtSSPair<<endl;
        cout<<"minDRll = "<<minDRll<<endl;
        cout<<"Lep3Pt = "<<Lep3Pt<<endl;

        cout<<"------------------"<<endl<<endl;
    }


    //--- COMPUTE BDT OUTPUT from weight files
    signal_3l_TT_MVA   = mva_3l_tt->EvaluateMVA("BDTG method");
    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}










//NEW -------------------------------------------
void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_Training(int evt)
{
    if(weight==0) {return;}

    if(vFakeableLeptons.size() < 3) {return;}

    // ##########
    // # Z veto #
    // ##########

    //-- Only consider the 3 hardest 'FakeableObject' leptons

    if(vFakeableLeptons.at(0).pt() < 20 || vFakeableLeptons.at(1).pt() < 10 || vFakeableLeptons.at(2).pt() < 10)
    {
        return;
    }

    if(vSelectedJets.size() < 2) {return;}

    if(nLooseBJets == 0) {return;}

    if(nForwardJets == 0) {return;}

    bool pass_Zveto = true;
    for(int i=0; i<3; i++)
    {
        for(int j=i+1; j<3; j++)
        {
            if(vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() && vFakeableLeptons.at(i).charge() == -vFakeableLeptons.at(j).charge() && fabs((vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4()).M() - 91.188) < 10 ) {pass_Zveto = false;}
        }
    }
    if(!pass_Zveto) {return;}

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vLooseLeptons.size(); i++)
    {
        lepton_px = lepton_px + vLooseLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vLooseLeptons.at(i).p4().Py();
    }

    TLorentzVector taup4;
    for(int i=0; i<vSelectedTaus.size(); i++)
    {
        taup4.SetPtEtaPhiE(vSelectedTaus.at(i).pt(), vSelectedTaus.at(i).eta(), vSelectedTaus.at(i).phi(), vSelectedTaus.at(i).E());
        tau_px = tau_px + taup4.Px();
        tau_py = tau_py + taup4.Py();
    }

    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        for(int j=0; j<vSelectedLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i).id() == -vSelectedLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;
    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    is_3l_THQ_Training = 1;



    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py


    //variables from tHq2016
    //NEW //FIXME
    //NB : "forward jet" = most forward non-CSV loose jet

    // double nJet25;
    // double maxEtaJet25;
    // double lepCharge;
    // double nJetEta1 ;
    // double dEtaFwdJetBJet;
    // double dEtaFwdJet2BJet;
    // double dEtaFwdJetClosestLep;
    // double dPhiHighestPtSSPair;
    // double minDRll;
    // double Lep3Pt;

    int ijet_forward=-1, ijet_hardest_btag=-1, ijet_2nd_hardest_btag=-1;
    double tmp = -999;

	//--- Find "forward jet"
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp)
        {
            tmp = fabs(vLightJets.at(ijet).eta());
            ijet_forward = ijet;
        }
    }

	//--- Find "hardest" and "second hardest" tagged jets
    double pt1=0, pt2=0;
    for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
    {
        if(fabs(vLooseBTagJets.at(ijet).pt()) > pt1)
        {
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            ijet_hardest_btag = ijet;
        }
        else if(fabs(vLooseBTagJets.at(ijet).pt()) > pt2)
        {
            pt2 = fabs(vLooseBTagJets.at(ijet).pt());
            ijet_2nd_hardest_btag = ijet;
        }
    }



	//--- Var1 : nof jets with pT>25 and |eta|<2.4
    nJet25 = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if(vSelectedJets.at(ijet).pt() > 25 && fabs(vSelectedJets.at(ijet).eta() ) < 2.4) {nJet25++;}
    }

    //--- Var2: max eta of any 'non-CSV-loose' jet
    tmp = -999;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if(fabs(vLightJets.at(ijet).eta()) > tmp) {tmp = fabs(vLightJets.at(ijet).eta());}
    }
    maxEtaJet25 = tmp;

	//--- Var3 : sum of leptons charges
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge() + vFakeableLeptons.at(2).charge();

	//--- Var4 : nof 'non-csv-loose' jets with eta>1.0
    nJetEta1 = 0;
    for(int ijet=0; ijet<vLightJets.size(); ijet++)
    {
        if( fabs(vLightJets.at(ijet).eta() ) > 1.0) {nJetEta1++;}
    }

	//--- Var5 : dEta between forward light jet and hardest tagged jet
    dEtaFwdJetBJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_hardest_btag).eta() );

	//--- Var6 : dEta between forward and 2nd hardest tagged jet
    if(ijet_2nd_hardest_btag < 0) {dEtaFwdJet2BJet = -1;}
    else {dEtaFwdJet2BJet = fabs( vLightJets.at(ijet_forward).eta() - vLooseBTagJets.at(ijet_2nd_hardest_btag).eta() );}

	//--- Var7 : dEta between forward light jet and closet lepton (angular dist.)
    dEtaFwdJetClosestLep = 0; tmp = 999;
    for(int ilep=0; ilep<3; ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair (?)
    if(vFakeableLeptons.at(0).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vSelectedLeptons.at(1).phi() - vFakeableLeptons.at(2).phi() ) );}
    else if(vFakeableLeptons.at(1).charge()!=lepCharge) {dPhiHighestPtSSPair = fabs( Phi_MPi_Pi(vFakeableLeptons.at(0).phi() - vFakeableLeptons.at(2).phi() ) );}
    else {dPhiHighestPtSSPair = fabs( Phi_MPi_Pi(vFakeableLeptons.at(0).phi() - vFakeableLeptons.at(1).phi() ) );}

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999; tmp = 999;
	TLorentzVector lep1, lep2, lep3;
	lep1.SetPtEtaPhiE(vFakeableLeptons.at(0).pt(), vFakeableLeptons.at(0).eta(), vFakeableLeptons.at(0).phi(), vFakeableLeptons.at(0).E() );
	lep2.SetPtEtaPhiE(vFakeableLeptons.at(1).pt(), vFakeableLeptons.at(1).eta(), vFakeableLeptons.at(1).phi(), vFakeableLeptons.at(1).E() );
	lep3.SetPtEtaPhiE(vFakeableLeptons.at(2).pt(), vFakeableLeptons.at(2).eta(), vFakeableLeptons.at(2).phi(), vFakeableLeptons.at(2).E() );
	if(lep1.DeltaR(lep2) < tmp) {minDRll = lep1.DeltaR(lep2); tmp = minDRll;}
	if(lep1.DeltaR(lep3) < tmp) {minDRll = lep1.DeltaR(lep3); tmp = minDRll;}
	if(lep2.DeltaR(lep3) < tmp) {minDRll = lep2.DeltaR(lep3); tmp = minDRll;}


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(2).pt();

    //--------------------------
    // PRINTOUT OF INFOS & INPUT VARS
    //--------------------------

    bool do_printout = false;

    if(do_printout)
    {
        cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

        cout<<FYEL("--- Tagged Jets : ")<<endl;
        for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt()<<" / eta = "<<vLooseBTagJets.at(ijet).eta()<<" / phi = "<<vLooseBTagJets.at(ijet).phi()<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2()<<endl;
        }
        cout<<"Hardest & 2nd hardest jets are "<<ijet_hardest_btag<<", "<<ijet_2nd_hardest_btag<<endl<<endl;

        cout<<FYEL("--- Forward Jets : ")<<endl;
        for(int ijet=0; ijet<vLightJets.size(); ijet++)
        {
            cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt()<<" / eta = "<<vLightJets.at(ijet).eta()<<" / phi = "<<vLightJets.at(ijet).phi()<<" / CSV = "<<vLightJets.at(ijet).CSVv2()<<endl;
        }
        cout<<"Forwardest jet is "<<ijet_forward<<endl<<endl;

        cout<<FYEL("--- Selected Leptons : ")<<endl;
        for(int ilep=0; ilep<vSelectedLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vSelectedLeptons.at(ilep).pt()<<" / eta = "<<vSelectedLeptons.at(ilep).eta()<<" / phi = "<<vSelectedLeptons.at(ilep).phi()<<" : Charge = "<<vSelectedLeptons.at(ilep).charge()<<endl;
        }

        cout<<FYEL("--- Input variables : ")<<endl;

        cout<<"nJet25 = "<<nJet25<<endl;
        cout<<"maxEtaJet25 = "<<maxEtaJet25<<endl;
        cout<<"lepCharge = "<<lepCharge<<endl;
        cout<<"nJetEta1 = "<<nJetEta1 <<endl;
        cout<<"dEtaFwdJetBJet = "<<dEtaFwdJetBJet<<endl;
        cout<<"dEtaFwdJet2BJet = "<<dEtaFwdJet2BJet<<endl;
        cout<<"dEtaFwdJetClosestLep = "<<dEtaFwdJetClosestLep<<endl;
        cout<<"dPhiHighestPtSSPair = "<<dPhiHighestPtSSPair<<endl;
        cout<<"minDRll = "<<minDRll<<endl;
        cout<<"Lep3Pt = "<<Lep3Pt<<endl;

        cout<<"------------------"<<endl<<endl;
    }


    //--- COMPUTE BDT OUTPUT from weight files
    signal_3l_TT_MVA   = mva_3l_tt->EvaluateMVA("BDTG method");
    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================

    fillOutputTree();

    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}














void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_ApplicationFakes(int evt)
{
    if(weight==0)         return; // For data not passing the relevant trigger, clean up the histograms from events with weight 0)
    //if(vSelectedTaus.size()>0) return;

    // ####################
    // # Common selection #
    // ####################

    if(DEBUG)             std::cout << std::endl << " 3l ss FR ==================================================================" << std::endl;

    bool nLep               = ( vSelectedLeptons.size()     >= 3 );
    if(!nLep)               return;

    if(DEBUG) std::cout << "nLep Ok... ";

    bool are_fullytight     =  ( vSelectedLeptons.at(0).isTightTTH()        && vSelectedLeptons.at(1).isTightTTH()      && vSelectedLeptons.at(2).isTightTTH()
                              //&& vSelectedLeptons.at(0).passTightCharge()   && vSelectedLeptons.at(1).passTightCharge() && vSelectedLeptons.at(2).passTightCharge()
                              && vSelectedLeptons.at(0).cutEventSel()       && vSelectedLeptons.at(1).cutEventSel()     && vSelectedLeptons.at(2).cutEventSel()
                              && vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()      && vSelectedLeptons.at(2).noLostHits()      );
    if(are_fullytight)     return;
    if(DEBUG) std::cout << "nNonTight Ok... ";

    bool pass_nolosthits    =  ( vSelectedLeptons.at(0).noLostHits()        && vSelectedLeptons.at(1).noLostHits()       && vSelectedLeptons.at(2).noLostHits()     );
    //if(!pass_nolosthits)     return;
    if(DEBUG) std::cout << "pass_nolosthits Ok... ";

    bool leading_lep_pt     = ( vSelectedLeptons.at(0).pt() > 25 );
    if(!leading_lep_pt)     return;

    if(DEBUG) std::cout << "leading_lep_pt Ok... ";

    bool following_lep_pt   = ( (vSelectedLeptons.at(1).pt() > 15) && (vSelectedLeptons.at(2).pt() > 15) );
    if(!following_lep_pt)   return;

    if(DEBUG) std::cout << "following_lep_pt Ok... ";

    bool nJets              = ( vSelectedJets.size()        >= 2 );
    if(!nJets)              return;

    if(DEBUG) std::cout << "nJets Ok... ";

    bool nLooseBtag         = ( nLooseBJets                 >= 2 );
    bool nMediumBtag        = ( nMediumBJets                >= 1 );
    if(!nLooseBtag && !nMediumBtag)      return;

    if(DEBUG) std::cout << "Btag Ok...";

    // Adding invariant mass cut on loose leptons pairs
    bool pass_invariantemasscut = true;
    for(int i=0; i<vLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vLeptons.size(); j++)
        {
            if ( fabs( ( vLeptons.at(i).p4() + vLeptons.at(j).p4() ).M() ) < 12 )
            { pass_invariantemasscut = false ;}
        }
    }
    if(!pass_invariantemasscut) return;

    if(DEBUG) std::cout << "invariantmasscut Ok... ";

    // ##########
    // # Z veto #
    // ##########

    bool pass_Zveto = true;
    for(int i=0; i<vSelectedLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            if (  ( vSelectedLeptons.at(i).id()     == -vSelectedLeptons.at(j).id()                             )
               && ( vSelectedLeptons.at(i).charge() == -vSelectedLeptons.at(j).charge()                         )
                    && ( fabs( ( vSelectedLeptons.at(i).p4() + vSelectedLeptons.at(j).p4() ).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
    }
    if(!pass_Zveto)       return;

    if(DEBUG) std::cout << "Zveto Ok... ";

    // ##########
    // # MET LD #
    // ##########

    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        lepton_px = lepton_px + vSelectedLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vSelectedLeptons.at(i).p4().Py();
    }

    MHT = sqrt( (jet_px + lepton_px) * (jet_px + lepton_px) + (jet_py + lepton_py) * (jet_py + lepton_py) );

    met_ld = 0.00397 * vEvent->at(0).metpt() + 0.00265 * MHT;

    // ########
    // # SFOS #
    // ########

    bool isSFOS = false;
    for(int i=0; i<vLeptons.size(); i++)
    {
        for(int j=0; j<vLeptons.size(); j++)
        {
            if ( i != j && ( vLeptons.at(i).id() == -vLeptons.at(j).id() ) ) {isSFOS = true;}
        }
    }

    if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges = 0;
    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        sum_charges = sum_charges + vSelectedLeptons.at(i).charge();
    }

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vSelectedLeptons.at(i).charge();
    }
    if( fabs(sum_charges_3l) != 1 ) return;

    if(DEBUG) std::cout << "sumOfCharges Ok... ";

    // ##################################################################################################################################

    // #################################
    // # b-tagging nominal reweighting #
    // #################################

    std::vector<double> jetPts;
    std::vector<double> jetEtas;
    std::vector<double> jetCSVs;
    std::vector<int>    jetFlavors;
    int iSys = 0;
    double wgt_csv, wgt_csv_def, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf, new_weight;

    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetPts.push_back(     vSelectedJets.at(i).pt()                );
        jetEtas.push_back(    vSelectedJets.at(i).eta()               );
        jetCSVs.push_back(    vSelectedJets.at(i).CSVv2()             );
        jetFlavors.push_back( vSelectedJets.at(i).jet_hadronFlavour() );
    }

    wgt_csv_def = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, iSys, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf);
    new_weight = weight * wgt_csv_def; // weight = weight * wgt_csv_def;

    // ##################################################################################################################################

    // ##################################
    // # b-tagging deriving systematics #
    // ##################################

    std::vector<double> weights_csv;
    double wgt_csv_def_sys = 0;

    for(int i=7; i<25; i++)
    {
        wgt_csv_def_sys = get_csv_wgt(jetPts, jetEtas, jetCSVs, jetFlavors, i, wgt_csv_hf, wgt_csv_lf, wgt_csv_cf)/wgt_csv_def;
        weights_csv.push_back(wgt_csv_def_sys);
    }

    double min_weight_csv = *min_element(weights_csv.begin(),weights_csv.end());
    double max_weight_csv = *max_element(weights_csv.begin(),weights_csv.end());

    weight_csv_down = min_weight_csv;
    weight_csv_up   = max_weight_csv;

    // ##################################################################################################################################

    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    for(int i=0; i<vSelectedLeptons.size(); i++)
    {
        if( !vSelectedLeptons.at(i).isTightTTH())
        {
            //std::cout << "Applying fake rate" << std::endl;
            //std::cout << "i: " << i << std::endl;
            //std::cout << "pt: " << vSelectedLeptons.at(i).pt() << std::endl;
            //std::cout << "fake? " << vSelectedLeptons.at(i).isFakeableTTH() << std::endl;
            leptonsPts.push_back(     vSelectedLeptons.at(i).pt()                );
            leptonsEtas.push_back(    vSelectedLeptons.at(i).eta()               );
            leptonsIds.push_back(     vSelectedLeptons.at(i).id()                );
        }
    }

    weightfake = get_FR_wgt_3l(leptonsPts, leptonsEtas, leptonsIds);

    // ##################################################################################################################################

    is_3l_AppFakes_SR = true;


    // ####################################
    // #  ____  ____    ____  ____ _____  #
    // # |___ \|  _ \  | __ )|  _ \_   _| #
    // #   __) | | | | |  _ \| | | || |   #
    // #  / __/| |_| | | |_) | |_| || |   #
    // # |_____|____/  |____/|____/ |_|   #
    // #                                  #
    // ####################################

    // based on https://github.com/CERN-PH-CMG/cmgtools-lite/blob/0b47d4d1c50ea0e24ef0d9cf1c24c763e78c1bf0/TTHAnalysis/python/tools/kinMVA_2D_2lss_3l.py

    // ======================================================================================================
    // variables against ttbar

    max_Lep_eta     = std::max( fabs(vSelectedLeptons.at(0).eta()), fabs(vSelectedLeptons.at(1).eta()) ) ;

    MT_met_lep1     = sqrt( 2 * vSelectedLeptons.at(0).p4().Pt() * vEvent->at(0).metpt() * (1 - cos( vSelectedLeptons.at(0).phi() - vEvent->at(0).metphi() )));

    nJet25_Recl     = vSelectedJets.size() ;
    mhtJet25_Recl   = sqrt( (jet_px*jet_px) + (jet_py*jet_py) );

    int njj = 0;
    avg_dr_jet = 0.;
    for (int ijet=0; ijet < vSelectedJets.size()-1 ; ijet++)
    {
        for (int kjet=ijet+1; kjet < vSelectedJets.size() ; kjet++)
        {
            avg_dr_jet += GetDeltaR( vSelectedJets.at(ijet).eta(), vSelectedJets.at(ijet).phi(), vSelectedJets.at(kjet).eta(), vSelectedJets.at(kjet).phi() );
            njj++;
        }
    }
    if ( njj > 0 ) avg_dr_jet = avg_dr_jet / njj;

    mindr_lep1_jet = 1000.;
    mindr_lep2_jet = 1000.;

    for (int i=0; i<vSelectedJets.size(); i++)
    {
        if( DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ) < mindr_lep1_jet )
        { mindr_lep1_jet = DeltaRLeptonJet( vSelectedLeptons.at(0), vSelectedJets.at(i) ); }
        if( DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ) < mindr_lep2_jet )
        { mindr_lep2_jet = DeltaRLeptonJet( vSelectedLeptons.at(2), vSelectedJets.at(i) ); }
    }

    signal_3l_TT_MVA    = mva_3l_tt->EvaluateMVA("BDTG method");


    // ======================================================================================================
    // variables against ttV

    met             = vEvent->at(0).metpt() ;

    LepGood_conePt0 = vSelectedLeptons.at(0).pt() ;
    LepGood_conePt1 = vSelectedLeptons.at(2).pt() ;

    signal_3l_TTV_MVA   = mva_3l_ttV->EvaluateMVA("BDTG method");


    //std::cout << " signal 3l   TT MVA: "  << signal_3l_TT_MVA
    //          << " signal 3l   TTV MVA: " << signal_3l_TTV_MVA << std::endl;

    // ======================================================================================================

    fillOutputTree();

    // if (_printLHCO_RECO) PrintLHCOforMadweight_RECO(evt);
}

// # ######################################################################
// #                  _             _                  _                  #
// #   ___ ___  _ __ | |_ _ __ ___ | |  _ __ ___  __ _(_) ___  _ __  ___  #
// #  / __/ _ \| '_ \| __| '__/ _ \| | | '__/ _ \/ _` | |/ _ \| '_ \/ __| #
// # | (_| (_) | | | | |_| | | (_) | | | | |  __/ (_| | | (_) | | | \__ \ #
// #  \___\___/|_| |_|\__|_|  \___/|_| |_|  \___|\__, |_|\___/|_| |_|___/ #
// #                                             |___/                    #
// #                                                                      #
// ########################################################################















void TTbarHiggsMultileptonAnalysis::initializeOutputTree()
{
    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

    tOutput->Branch("mc_event",&mc_event,"mc_event/I");
    tOutput->Branch("weight",&weight,"weight/F");
    tOutput->Branch("weightfake",&weightfake,"weightfake/F");
    tOutput->Branch("weightflip",&weightflip,"weightflip/F");
    tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");
    tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
    tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
    tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
    tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
    tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
    tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");
    tOutput->Branch("weights_pdf","std::vector<float>",&weights_pdf);
    tOutput->Branch("ids_pdf","std::vector<std::string>",&ids_pdf);

    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");
    tOutput->Branch("mc_3l_category",&mc_3l_category,"mc_3l_category/I");
    tOutput->Branch("mc_ttbar_decay",&mc_ttbar_decay,"mc_ttbar_decay/I");
    tOutput->Branch("mc_boson_decay",&mc_boson_decay,"mc_boson_decay/I");
    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
    tOutput->Branch("mc_nJets25",&mc_nJets25,"mc_nJets25/I");
    tOutput->Branch("mc_nBtagJets25",&mc_nBtagJets25,"mc_nBtagJets25/I");

    tOutput->Branch("catJets",&catJets,"catJets/I");

    //-- NB : coded as floats bc treated as any other variable in analysis code (floats) !
    tOutput->Branch("is_3l_THQ_SR",&is_3l_THQ_SR,"is_3l_THQ_SR/F");
    tOutput->Branch("is_3l_THQ_Training",&is_3l_THQ_Training,"is_3l_THQ_Training/F");

    tOutput->Branch("cat_HtoWW",&cat_HtoWW,"cat_HtoWW/B");
    tOutput->Branch("cat_HtoZZ",&cat_HtoZZ,"cat_HtoZZ/B");
    tOutput->Branch("cat_Htott",&cat_Htott,"cat_Htott/B");

    tOutput->Branch("is_trigger",&is_trigger,"is_trigger/B");

    tOutput->Branch("max_Lep_eta", &max_Lep_eta, "max_Lep_eta/F");
    tOutput->Branch("MT_met_lep1",&MT_met_lep1,"MT_met_lep1/F");
    tOutput->Branch("nJet25_Recl",&nJet25_Recl,"nJet25_Recl/F");
    tOutput->Branch("mindr_lep1_jet",&mindr_lep1_jet,"mindr_lep1_jet/F");
    tOutput->Branch("mindr_lep2_jet",&mindr_lep2_jet,"mindr_lep2_jet/F");
    tOutput->Branch("LepGood_conePt0",&LepGood_conePt0,"LepGood_conePt0/F");
    tOutput->Branch("LepGood_conePt1",&LepGood_conePt1,"LepGood_conePt1/F");
    tOutput->Branch("met",&met,"met/F");
    tOutput->Branch("mhtJet25_Recl",&mhtJet25_Recl,"mhtJet25_Recl/F");
    tOutput->Branch("avg_dr_jet",&avg_dr_jet,"avg_dr_jet/F");

    // tOutput->Branch("signal_2lss_TT_MVA",&signal_2lss_TT_MVA,"signal_2lss_TT_MVA/F");
    // tOutput->Branch("signal_2lss_TTV_MVA",&signal_2lss_TTV_MVA,"signal_2lss_TTV_MVA/F");

    tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
    tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");

    tOutput->Branch("multilepton_Lepton1_Id",               &multilepton_Lepton1_Id,                "multilepton_Lepton1_Id/I");
    tOutput->Branch("multilepton_Lepton1_P4",               "TLorentzVector",                       &multilepton_Lepton1_P4);
    tOutput->Branch("multilepton_Lepton1_DeltaR_Matched",   &multilepton_Lepton1_DeltaR_Matched,    "multilepton_Lepton1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton1_Label_Matched",    &multilepton_Lepton1_Label_Matched,     "multilepton_Lepton1_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton1_Id_Matched",       &multilepton_Lepton1_Id_Matched,        "multilepton_Lepton1_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton1_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton1_P4_Matched);
    tOutput->Branch("multilepton_Lepton2_Id",               &multilepton_Lepton2_Id,                "multilepton_Lepton2_Id/I");
    tOutput->Branch("multilepton_Lepton2_P4",               "TLorentzVector",                       &multilepton_Lepton2_P4);
    tOutput->Branch("multilepton_Lepton2_DeltaR_Matched",   &multilepton_Lepton2_DeltaR_Matched,    "multilepton_Lepton2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton2_Label_Matched",    &multilepton_Lepton2_Label_Matched,     "multilepton_Lepton2_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton2_Id_Matched",       &multilepton_Lepton2_Id_Matched,        "multilepton_Lepton2_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton2_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton2_P4_Matched);
    tOutput->Branch("multilepton_Lepton3_Id",               &multilepton_Lepton3_Id,                "multilepton_Lepton3_Id/I");
    tOutput->Branch("multilepton_Lepton3_P4",               "TLorentzVector",                       &multilepton_Lepton3_P4);
    tOutput->Branch("multilepton_Lepton3_DeltaR_Matched",   &multilepton_Lepton3_DeltaR_Matched,    "multilepton_Lepton3_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton3_Label_Matched",    &multilepton_Lepton3_Label_Matched,     "multilepton_Lepton3_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton3_Id_Matched",       &multilepton_Lepton3_Id_Matched,        "multilepton_Lepton3_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton3_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton3_P4_Matched);
    tOutput->Branch("multilepton_Lepton4_Id",               &multilepton_Lepton4_Id,                "multilepton_Lepton4_Id/I");
    tOutput->Branch("multilepton_Lepton4_P4",               "TLorentzVector",                       &multilepton_Lepton4_P4);
    tOutput->Branch("multilepton_Lepton4_DeltaR_Matched",   &multilepton_Lepton4_DeltaR_Matched,    "multilepton_Lepton4_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Lepton4_Label_Matched",    &multilepton_Lepton4_Label_Matched,     "multilepton_Lepton4_Label_Matched/I");
    tOutput->Branch("multilepton_Lepton4_Id_Matched",       &multilepton_Lepton4_Id_Matched,        "multilepton_Lepton4_Id_Matched/I");
    tOutput->Branch("multilepton_Lepton4_P4_Matched",       "TLorentzVector",                       &multilepton_Lepton4_P4_Matched);

    tOutput->Branch("multilepton_Bjet1_Id",                 &multilepton_Bjet1_Id,                  "multilepton_Bjet1_Id/I");
    tOutput->Branch("multilepton_Bjet1_P4",                 "TLorentzVector",                       &multilepton_Bjet1_P4);
    tOutput->Branch("multilepton_Bjet1_CSV",                &multilepton_Bjet1_CSV,                 "multilepton_Bjet1_CSV/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Up",             &multilepton_Bjet1_JEC_Up,              "multilepton_Bjet1_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet1_JEC_Down",           &multilepton_Bjet1_JEC_Down,            "multilepton_Bjet1_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet1_JER_Up",             &multilepton_Bjet1_JER_Up,              "multilepton_Bjet1_JER_Up/F");
    tOutput->Branch("multilepton_Bjet1_JER_Down",           &multilepton_Bjet1_JER_Down,            "multilepton_Bjet1_JER_Down/F");
    tOutput->Branch("multilepton_Bjet1_DeltaR_Matched",     &multilepton_Bjet1_DeltaR_Matched,      "multilepton_Bjet1_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet1_Label_Matched",      &multilepton_Bjet1_Label_Matched,       "multilepton_Bjet1_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet1_Id_Matched",         &multilepton_Bjet1_Id_Matched,          "multilepton_Bjet1_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet1_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet1_P4_Matched);

    tOutput->Branch("multilepton_Bjet2_Id",                 &multilepton_Bjet2_Id,                  "multilepton_Bjet2_Id/I");
    tOutput->Branch("multilepton_Bjet2_P4",                 "TLorentzVector",                       &multilepton_Bjet2_P4);
    tOutput->Branch("multilepton_Bjet2_CSV",                &multilepton_Bjet2_CSV,                 "multilepton_Bjet2_CSV/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Up",             &multilepton_Bjet2_JEC_Up,              "multilepton_Bjet2_JEC_Up/F");
    tOutput->Branch("multilepton_Bjet2_JEC_Down",           &multilepton_Bjet2_JEC_Down,            "multilepton_Bjet2_JEC_Down/F");
    tOutput->Branch("multilepton_Bjet2_JER_Up",             &multilepton_Bjet2_JER_Up,              "multilepton_Bjet2_JER_Up/F");
    tOutput->Branch("multilepton_Bjet2_JER_Down",           &multilepton_Bjet2_JER_Down,            "multilepton_Bjet2_JER_Down/F");
    tOutput->Branch("multilepton_Bjet2_DeltaR_Matched",     &multilepton_Bjet2_DeltaR_Matched,      "multilepton_Bjet2_DeltaR_Matched/F");
    tOutput->Branch("multilepton_Bjet2_Label_Matched",      &multilepton_Bjet2_Label_Matched,       "multilepton_Bjet2_Label_Matched/I");
    tOutput->Branch("multilepton_Bjet2_Id_Matched",         &multilepton_Bjet2_Id_Matched,          "multilepton_Bjet2_Id_Matched/I");
    tOutput->Branch("multilepton_Bjet2_P4_Matched",         "TLorentzVector",                       &multilepton_Bjet2_P4_Matched);

    tOutput->Branch("multilepton_JetHighestPt1_Id",&multilepton_JetHighestPt1_Id,"multilepton_JetHighestPt1_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_P4","TLorentzVector",&multilepton_JetHighestPt1_P4);
    tOutput->Branch("multilepton_JetHighestPt1_CSV",&multilepton_JetHighestPt1_CSV,"multilepton_JetHighestPt1_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Up",&multilepton_JetHighestPt1_JEC_Up,"multilepton_JetHighestPt1_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JEC_Down",&multilepton_JetHighestPt1_JEC_Down,"multilepton_JetHighestPt1_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Up",&multilepton_JetHighestPt1_JER_Up,"multilepton_JetHighestPt1_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_JER_Down",&multilepton_JetHighestPt1_JER_Down,"multilepton_JetHighestPt1_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_Id",&multilepton_JetHighestPt2_Id,"multilepton_JetHighestPt2_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_P4","TLorentzVector",&multilepton_JetHighestPt2_P4);
    tOutput->Branch("multilepton_JetHighestPt2_CSV",&multilepton_JetHighestPt2_CSV,"multilepton_JetHighestPt2_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Up",&multilepton_JetHighestPt2_JEC_Up,"multilepton_JetHighestPt2_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JEC_Down",&multilepton_JetHighestPt2_JEC_Down,"multilepton_JetHighestPt2_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Up",&multilepton_JetHighestPt2_JER_Up,"multilepton_JetHighestPt2_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_JER_Down",&multilepton_JetHighestPt2_JER_Down,"multilepton_JetHighestPt2_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_Id",&multilepton_JetClosestMw1_Id,"multilepton_JetClosestMw1_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_P4","TLorentzVector",&multilepton_JetClosestMw1_P4);
    tOutput->Branch("multilepton_JetClosestMw1_CSV",&multilepton_JetClosestMw1_CSV,"multilepton_JetClosestMw1_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Up",&multilepton_JetClosestMw1_JEC_Up,"multilepton_JetClosestMw1_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JEC_Down",&multilepton_JetClosestMw1_JEC_Down,"multilepton_JetClosestMw1_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Up",&multilepton_JetClosestMw1_JER_Up,"multilepton_JetClosestMw1_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_JER_Down",&multilepton_JetClosestMw1_JER_Down,"multilepton_JetClosestMw1_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_Id",&multilepton_JetClosestMw2_Id,"multilepton_JetClosestMw2_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_P4","TLorentzVector",&multilepton_JetClosestMw2_P4);
    tOutput->Branch("multilepton_JetClosestMw2_CSV",&multilepton_JetClosestMw2_CSV,"multilepton_JetClosestMw2_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Up",&multilepton_JetClosestMw2_JEC_Up,"multilepton_JetClosestMw2_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JEC_Down",&multilepton_JetClosestMw2_JEC_Down,"multilepton_JetClosestMw2_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Up",&multilepton_JetClosestMw2_JER_Up,"multilepton_JetClosestMw2_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_JER_Down",&multilepton_JetClosestMw2_JER_Down,"multilepton_JetClosestMw2_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_Id",&multilepton_JetLowestMjj1_Id,"multilepton_JetLowestMjj1_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_P4","TLorentzVector",&multilepton_JetLowestMjj1_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_CSV",&multilepton_JetLowestMjj1_CSV,"multilepton_JetLowestMjj1_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Up",&multilepton_JetLowestMjj1_JEC_Up,"multilepton_JetLowestMjj1_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JEC_Down",&multilepton_JetLowestMjj1_JEC_Down,"multilepton_JetLowestMjj1_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Up",&multilepton_JetLowestMjj1_JER_Up,"multilepton_JetLowestMjj1_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_JER_Down",&multilepton_JetLowestMjj1_JER_Down,"multilepton_JetLowestMjj1_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_Id",&multilepton_JetLowestMjj2_Id,"multilepton_JetLowestMjj2_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_P4","TLorentzVector",&multilepton_JetLowestMjj2_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_CSV",&multilepton_JetLowestMjj2_CSV,"multilepton_JetLowestMjj2_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Up",&multilepton_JetLowestMjj2_JEC_Up,"multilepton_JetLowestMjj2_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JEC_Down",&multilepton_JetLowestMjj2_JEC_Down,"multilepton_JetLowestMjj2_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Up",&multilepton_JetLowestMjj2_JER_Up,"multilepton_JetLowestMjj2_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_JER_Down",&multilepton_JetLowestMjj2_JER_Down,"multilepton_JetLowestMjj2_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_Id",&multilepton_JetHighestPt1_2ndPair_Id,"multilepton_JetHighestPt1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt1_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_CSV",&multilepton_JetHighestPt1_2ndPair_CSV,"multilepton_JetHighestPt1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Up",&multilepton_JetHighestPt1_2ndPair_JEC_Up,"multilepton_JetHighestPt1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JEC_Down",&multilepton_JetHighestPt1_2ndPair_JEC_Down,"multilepton_JetHighestPt1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Up",&multilepton_JetHighestPt1_2ndPair_JER_Up,"multilepton_JetHighestPt1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt1_2ndPair_JER_Down",&multilepton_JetHighestPt1_2ndPair_JER_Down,"multilepton_JetHighestPt1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_Id",&multilepton_JetHighestPt2_2ndPair_Id,"multilepton_JetHighestPt2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_P4","TLorentzVector",&multilepton_JetHighestPt2_2ndPair_P4);
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_CSV",&multilepton_JetHighestPt2_2ndPair_CSV,"multilepton_JetHighestPt2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Up",&multilepton_JetHighestPt2_2ndPair_JEC_Up,"multilepton_JetHighestPt2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JEC_Down",&multilepton_JetHighestPt2_2ndPair_JEC_Down,"multilepton_JetHighestPt2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Up",&multilepton_JetHighestPt2_2ndPair_JER_Up,"multilepton_JetHighestPt2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetHighestPt2_2ndPair_JER_Down",&multilepton_JetHighestPt2_2ndPair_JER_Down,"multilepton_JetHighestPt2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_Id",&multilepton_JetClosestMw1_2ndPair_Id,"multilepton_JetClosestMw1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw1_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_CSV",&multilepton_JetClosestMw1_2ndPair_CSV,"multilepton_JetClosestMw1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Up",&multilepton_JetClosestMw1_2ndPair_JEC_Up,"multilepton_JetClosestMw1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JEC_Down",&multilepton_JetClosestMw1_2ndPair_JEC_Down,"multilepton_JetClosestMw1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Up",&multilepton_JetClosestMw1_2ndPair_JER_Up,"multilepton_JetClosestMw1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw1_2ndPair_JER_Down",&multilepton_JetClosestMw1_2ndPair_JER_Down,"multilepton_JetClosestMw1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_Id",&multilepton_JetClosestMw2_2ndPair_Id,"multilepton_JetClosestMw2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_P4","TLorentzVector",&multilepton_JetClosestMw2_2ndPair_P4);
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_CSV",&multilepton_JetClosestMw2_2ndPair_CSV,"multilepton_JetClosestMw2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Up",&multilepton_JetClosestMw2_2ndPair_JEC_Up,"multilepton_JetClosestMw2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JEC_Down",&multilepton_JetClosestMw2_2ndPair_JEC_Down,"multilepton_JetClosestMw2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Up",&multilepton_JetClosestMw2_2ndPair_JER_Up,"multilepton_JetClosestMw2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetClosestMw2_2ndPair_JER_Down",&multilepton_JetClosestMw2_2ndPair_JER_Down,"multilepton_JetClosestMw2_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_Id",&multilepton_JetLowestMjj1_2ndPair_Id,"multilepton_JetLowestMjj1_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj1_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_CSV",&multilepton_JetLowestMjj1_2ndPair_CSV,"multilepton_JetLowestMjj1_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Up",&multilepton_JetLowestMjj1_2ndPair_JEC_Up,"multilepton_JetLowestMjj1_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JEC_Down",&multilepton_JetLowestMjj1_2ndPair_JEC_Down,"multilepton_JetLowestMjj1_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Up",&multilepton_JetLowestMjj1_2ndPair_JER_Up,"multilepton_JetLowestMjj1_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj1_2ndPair_JER_Down",&multilepton_JetLowestMjj1_2ndPair_JER_Down,"multilepton_JetLowestMjj1_2ndPair_JER_Down/F");

    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_Id",&multilepton_JetLowestMjj2_2ndPair_Id,"multilepton_JetLowestMjj2_2ndPair_Id/I");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_P4","TLorentzVector",&multilepton_JetLowestMjj2_2ndPair_P4);
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_CSV",&multilepton_JetLowestMjj2_2ndPair_CSV,"multilepton_JetLowestMjj2_2ndPair_CSV/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Up",&multilepton_JetLowestMjj2_2ndPair_JEC_Up,"multilepton_JetLowestMjj2_2ndPair_JEC_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JEC_Down",&multilepton_JetLowestMjj2_2ndPair_JEC_Down,"multilepton_JetLowestMjj2_2ndPair_JEC_Down/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Up",&multilepton_JetLowestMjj2_2ndPair_JER_Up,"multilepton_JetLowestMjj2_2ndPair_JER_Up/F");
    tOutput->Branch("multilepton_JetLowestMjj2_2ndPair_JER_Down",&multilepton_JetLowestMjj2_2ndPair_JER_Down,"multilepton_JetLowestMjj2_2ndPair_JER_Down/F");

    // Test adding truth information

    tOutput->Branch("multilepton_h0_Id",                    &multilepton_h0_Id,                 "multilepton_h0_Id/I");
    tOutput->Branch("multilepton_h0_P4",                    "TLorentzVector",                   &multilepton_h0_P4);
    tOutput->Branch("multilepton_h0_Label",                    &multilepton_h0_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t1_Id",                    &multilepton_t1_Id,                 "multilepton_t1_Id/I");
    tOutput->Branch("multilepton_t1_P4",                    "TLorentzVector",                   &multilepton_t1_P4);
    tOutput->Branch("multilepton_t1_Label",                    &multilepton_t1_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t2_Id",                    &multilepton_t2_Id,                 "multilepton_t2_Id/I");
    tOutput->Branch("multilepton_t2_P4",                    "TLorentzVector",                   &multilepton_t2_P4);
    tOutput->Branch("multilepton_t2_Label",                    &multilepton_t2_Label,                 "multilepton_h0_Label/I");

    // End test adding truth information

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
    tOutput->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
    tOutput->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
    tOutput->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");
    tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    //-- Input variables from tHq2016 analysis
    tOutput->Branch("nJet25",&nJet25,"nJet25/F");
    tOutput->Branch("maxEtaJet25",&maxEtaJet25,"maxEtaJet25/F");
    tOutput->Branch("lepCharge",&lepCharge,"lepCharge/F");
    tOutput->Branch("nJetEta1",&nJetEta1,"nJetEta1/F");
    tOutput->Branch("dEtaFwdJetBJet",&dEtaFwdJetBJet,"dEtaFwdJetBJet/F");
    tOutput->Branch("dEtaFwdJet2BJet",&dEtaFwdJet2BJet,"dEtaFwdJet2BJet/F");
    tOutput->Branch("dEtaFwdJetClosestLep",&dEtaFwdJetClosestLep,"dEtaFwdJetClosestLep/F");
    tOutput->Branch("dPhiHighestPtSSPair",&dPhiHighestPtSSPair,"dPhiHighestPtSSPair/F");
    tOutput->Branch("minDRll",&minDRll,"minDRll/F");
    tOutput->Branch("Lep3Pt",&Lep3Pt,"Lep3Pt/F");

    return;
}

void TTbarHiggsMultileptonAnalysis::fillOutputTree()
{
    int tot_charge = 0;
    int tot_id = 0;
    if (vSelectedLeptons.size()>=4) {
        for (unsigned int i=0; i<4; i++) {
            tot_charge += vSelectedLeptons.at(i).charge();
            tot_id += vSelectedLeptons.at(i).id();
        }
    }

    //Choosing 2 b-jets
    bool doSelectOnlyBjets = true;

    //std::cout << "doSelectOnlyBjets = true" << std::endl;

    TLorentzVector Bjet1, Bjet2;
    int ib1=-1, ib2=-1;
    selectBjets("HighestBtagDiscrim", &ib1, &ib2, doSelectOnlyBjets);
    if (ib1!=-1) Bjet1.SetPtEtaPhiE(vSelectedJets.at(ib1).pt(), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi(), vSelectedJets.at(ib1).E());
    if (ib2!=-1) Bjet2.SetPtEtaPhiE(vSelectedJets.at(ib2).pt(), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi(), vSelectedJets.at(ib2).E());

    //std::cout << "Setting bjets ok" << std::endl;

    //3l
    if (ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=2) catJets = kCat_3l_2b_2j;
    else if (ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=2) catJets = kCat_3l_1b_2j;
    else if (ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) catJets = kCat_3l_2b_1j;
    else if (ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) catJets = kCat_3l_1b_1j;
    else if (ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==0) catJets = kCat_3l_2b_0j;
    else catJets = -1;

    //std::cout << "catJets="<<catJets<<std::endl;

    multilepton_Lepton1_Id = -999;
    multilepton_Lepton2_Id = -999;
    multilepton_Lepton3_Id = -999;
    multilepton_Lepton4_Id = -999;

    if (vSelectedLeptons.size()>=2){
        multilepton_Lepton1_P4 = vSelectedLeptons.at(0).p4();
        multilepton_Lepton1_Id = vSelectedLeptons.at(0).id();
        multilepton_Lepton2_P4 = vSelectedLeptons.at(1).p4();
        multilepton_Lepton2_Id = vSelectedLeptons.at(1).id();
    }

    if (vSelectedLeptons.size()>=3)
    {
        multilepton_Lepton3_P4 = vSelectedLeptons.at(2).p4();
        multilepton_Lepton3_Id = vSelectedLeptons.at(2).id();
    }


    //std::cout << "Setting leptons ok" << std::endl;

    multilepton_Bjet1_Id = -999;
    if (ib1!=-1){
        FillJetInfoOutputTree(&multilepton_Bjet1_Id, 5, &multilepton_Bjet1_P4, Bjet1, &multilepton_Bjet1_CSV, vSelectedJets.at(ib1).CSVv2(), &multilepton_Bjet1_JEC_Up, &multilepton_Bjet1_JEC_Down, vSelectedJets.at(ib1).JES_uncert(), &multilepton_Bjet1_JER_Up, &multilepton_Bjet1_JER_Down, vSelectedJets.at(ib1).pt_JER(), vSelectedJets.at(ib1).pt_JER_up(), vSelectedJets.at(ib1).pt_JER_down());
        //multilepton_Bjet1_P4 = Bjet1;
        //multilepton_Bjet1_Id = 5;
    }
    multilepton_Bjet2_Id = -999;
    if (ib2!=-1){
        FillJetInfoOutputTree(&multilepton_Bjet2_Id, 5, &multilepton_Bjet2_P4, Bjet2, &multilepton_Bjet2_CSV, vSelectedJets.at(ib2).CSVv2(), &multilepton_Bjet2_JEC_Up, &multilepton_Bjet2_JEC_Down, vSelectedJets.at(ib2).JES_uncert(), &multilepton_Bjet2_JER_Up, &multilepton_Bjet2_JER_Down, vSelectedJets.at(ib2).pt_JER(), vSelectedJets.at(ib2).pt_JER_up(), vSelectedJets.at(ib2).pt_JER_down());
        //multilepton_Bjet2_P4 = Bjet2;
        //multilepton_Bjet2_Id = 5;
    }


    // ###############################################################################
    // #                  _       _     _               _                            #
    // #  _ __ ___   __ _| |_ ___| |__ (_)_ __   __ _  | |_ ___     __ _  ___ _ __   #
    // # | '_ ` _ \ / _` | __/ __| '_ \| | '_ \ / _` | | __/ _ \   / _` |/ _ \ '_ \  #
    // # | | | | | | (_| | || (__| | | | | | | | (_| | | || (_) | | (_| |  __/ | | | #
    // # |_| |_| |_|\__,_|\__\___|_| |_|_|_| |_|\__, |  \__\___/   \__, |\___|_| |_| #
    // #                                        |___/              |___/             #
    // #                                                                             #
    // ###############################################################################

    if ( !_isdata )
    {
        float lep1_dr_gen       = 100.,     lep2_dr_gen     = 100.,     lep3_dr_gen     = 100.,     lep4_dr_gen     = 100. ;
        float jet1_dr_gen       = 100.,     jet2_dr_gen     = 100.;
        float lep1_dr_gen_min   = 100.,     lep2_dr_gen_min = 100.,     lep3_dr_gen_min = 100.,     lep4_dr_gen_min = 100. ;
        float jet1_dr_gen_min   = 100.,     jet2_dr_gen_min = 100.;
        int   lep1_matched      = -1,       lep2_matched    = -1,       lep3_matched   = -1,       lep4_matched    = -1;
        int   jet1_matched      = -1,       jet2_matched    = -1;

        TLorentzVector LeptonX;

        for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
        {
            if( abs(vTruth->at(0).mc_truth_id().at(itruth)) < 18 )
            {
                // cout<<"vTruth->at(0).mc_truth_eta().at(itruth) "<<vTruth->at(0).mc_truth_eta().at(itruth)<<endl;
                // cout<<"vTruth->at(0).mc_truth_phi().at(itruth) "<<vTruth->at(0).mc_truth_phi().at(itruth)<<endl;

                //FIXME -- some values of vTruth->at(0).mc_truth_eta() above 10^24 --> freeze program !
                // if(vTruth->at(0).mc_truth_eta().at(itruth) > pow(10, 6) || vTruth->at(0).mc_truth_phi().at(itruth) > pow(10, 6))
                // {
                //     cout<<"=== cf. line : "<<__LINE__<<endl;
                //     continue;
                // }

                lep1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(0).eta(), vSelectedLeptons.at(0).phi() );
                if( lep1_dr_gen < lep1_dr_gen_min)
                {
                    lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                }

                lep2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(1).eta(), vSelectedLeptons.at(1).phi() );
                if( lep2_dr_gen < lep2_dr_gen_min)
                {   lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;  }

                if(vSelectedLeptons.size()>=3)
                {
                    lep3_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(2).eta(), vSelectedLeptons.at(2).phi() );
                    if( lep3_dr_gen < lep3_dr_gen_min)
                    {   lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;  }
                }

                if(vSelectedLeptons.size()>=4)
                {
                    lep4_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedLeptons.at(3).eta(), vSelectedLeptons.at(3).phi() );
                    if( lep4_dr_gen < lep4_dr_gen_min)
                    {   lep4_dr_gen_min = lep4_dr_gen;  lep4_matched = itruth;  }
                }

                jet1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedJets.at(ib1).eta(), vSelectedJets.at(ib1).phi() );
                if( jet1_dr_gen < jet1_dr_gen_min)
                {   jet1_dr_gen_min = jet1_dr_gen;  jet1_matched = itruth;  }


                if(ib2!=-1)
                {
                    jet2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vSelectedJets.at(ib2).eta(), vSelectedJets.at(ib2).phi() );
                    if( jet2_dr_gen < jet2_dr_gen_min)
                    {   jet2_dr_gen_min = jet2_dr_gen;  jet2_matched = itruth;  }
                }
            }
        }


        if(false)
        {
            std::cout << "lep1_matched: "   << lep1_matched                                     << std::endl;
            std::cout << "pt: "             << vTruth->at(0).mc_truth_pt().at(lep1_matched)     << std::endl;
            std::cout << "eta: "            << vTruth->at(0).mc_truth_eta().at(lep1_matched)    << std::endl;
            std::cout << "phi: "            << vTruth->at(0).mc_truth_phi().at(lep1_matched)    << std::endl;
            std::cout << "E: "              << vTruth->at(0).mc_truth_E().at(lep1_matched)      << std::endl;
            std::cout << "Id: "             << vTruth->at(0).mc_truth_id().at(lep1_matched)     << std::endl;
            std::cout << "Label: "          << vTruth->at(0).mc_truth_label().at(lep1_matched)  << std::endl;
        }


        if(lep1_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep1_matched),       vTruth->at(0).mc_truth_eta().at(lep1_matched),
                                    vTruth->at(0).mc_truth_phi().at(lep1_matched),      vTruth->at(0).mc_truth_E().at(lep1_matched)     );
            multilepton_Lepton1_P4_Matched      = LeptonX;
            multilepton_Lepton1_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep1_matched);
            multilepton_Lepton1_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep1_matched);
            multilepton_Lepton1_DeltaR_Matched  = lep1_dr_gen_min;
        }

        if(lep2_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep2_matched),       vTruth->at(0).mc_truth_eta().at(lep2_matched),
                                    vTruth->at(0).mc_truth_phi().at(lep2_matched),      vTruth->at(0).mc_truth_E().at(lep2_matched)     );
            multilepton_Lepton2_P4_Matched      = LeptonX;
            multilepton_Lepton2_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep2_matched);
            multilepton_Lepton2_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep2_matched);
            multilepton_Lepton2_DeltaR_Matched  = lep2_dr_gen_min;
        }

        if(lep3_matched >= 0)
        {
            if(vSelectedLeptons.size()>=3)
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep3_matched),       vTruth->at(0).mc_truth_eta().at(lep3_matched),
                                        vTruth->at(0).mc_truth_phi().at(lep3_matched),      vTruth->at(0).mc_truth_E().at(lep3_matched)     );
                multilepton_Lepton3_P4_Matched      = LeptonX;
                multilepton_Lepton3_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep3_matched);
                multilepton_Lepton3_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep3_matched);
                multilepton_Lepton3_DeltaR_Matched  = lep3_dr_gen_min;
            }
        }

        if(lep4_matched >= 0)
        {
            if(vSelectedLeptons.size()>=4)
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(lep4_matched),       vTruth->at(0).mc_truth_eta().at(lep4_matched),
                                        vTruth->at(0).mc_truth_phi().at(lep4_matched),      vTruth->at(0).mc_truth_E().at(lep4_matched)     );
                multilepton_Lepton4_P4_Matched      = LeptonX;
                multilepton_Lepton4_Id_Matched      = vTruth->at(0).mc_truth_id().at(lep4_matched);
                multilepton_Lepton4_Label_Matched   = vTruth->at(0).mc_truth_label().at(lep4_matched);
                multilepton_Lepton4_DeltaR_Matched  = lep4_dr_gen_min;
            }
        }

        if( DEBUG ) std::cout << "Matching lepton to gen ok" << std::endl;

        if(false)
        {
            std::cout << " ============ "   << std::endl;
            std::cout << "jet1_matched: "   << jet1_matched                                     << std::endl;
            std::cout << "pt: "             << vTruth->at(0).mc_truth_pt().at(jet1_matched)     << std::endl;
            std::cout << "eta: "            << vTruth->at(0).mc_truth_eta().at(jet1_matched)    << std::endl;
            std::cout << "phi: "            << vTruth->at(0).mc_truth_phi().at(jet1_matched)    << std::endl;
            std::cout << "E: "              << vTruth->at(0).mc_truth_E().at(jet1_matched)      << std::endl;
            std::cout << "Id: "             << vTruth->at(0).mc_truth_id().at(jet1_matched)     << std::endl;
            std::cout << "Label: "          << vTruth->at(0).mc_truth_label().at(jet1_matched)  << std::endl;
        }

        if(jet1_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(jet1_matched),       vTruth->at(0).mc_truth_eta().at(jet1_matched),
                                    vTruth->at(0).mc_truth_phi().at(jet1_matched),      vTruth->at(0).mc_truth_E().at(jet1_matched)     );
            multilepton_Bjet1_P4_Matched        = LeptonX;
            multilepton_Bjet1_Id_Matched        = vTruth->at(0).mc_truth_id().at(jet1_matched);
            multilepton_Bjet1_Label_Matched     = vTruth->at(0).mc_truth_label().at(jet1_matched);
            multilepton_Bjet1_DeltaR_Matched    = jet1_dr_gen_min;
        }

        if(jet2_matched >= 0)
        {
            LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(jet2_matched),       vTruth->at(0).mc_truth_eta().at(jet2_matched),
                                    vTruth->at(0).mc_truth_phi().at(jet2_matched),      vTruth->at(0).mc_truth_E().at(jet2_matched)     );
            multilepton_Bjet2_P4_Matched        = LeptonX;
            multilepton_Bjet2_Id_Matched        = vTruth->at(0).mc_truth_id().at(jet2_matched);
            multilepton_Bjet2_Label_Matched     = vTruth->at(0).mc_truth_label().at(jet2_matched);
            multilepton_Bjet2_DeltaR_Matched    = jet2_dr_gen_min;
        }

    }

    // ========================

    //Choose 2 jets
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
    float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
    for (unsigned int ij=0; ij<vSelectedJets.size(); ij++){
        if (ij==ib1 || ij==ib2) continue;
        if (vSelectedJets.at(ij).pt() > pt_max ) {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vSelectedJets.at(ij).pt();
            ij1 = ij;
        }
        if (vSelectedJets.at(ij).pt() < pt_max && vSelectedJets.at(ij).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(ij).pt();
            ij2 = ij;
        }
        for (unsigned int ik=0; ik<vSelectedJets.size(); ik++){
            if (ik==ij) continue;
            if (ik==ib1 || ik==ib2) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(ij).pt(), vSelectedJets.at(ij).eta(), vSelectedJets.at(ij).phi(), vSelectedJets.at(ij).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(ik).pt(), vSelectedJets.at(ik).eta(), vSelectedJets.at(ik).phi(), vSelectedJets.at(ik).E());
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                ik1=ij;
                ik2=ik;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if ((Pjet1+Pjet2).M()<mass_min){
                il1=ij;
                il2=ik;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }

    //std::cout << "Choosing 2 jets ok" << std::endl;

    //Choose 2 more jets
    int io1=-1, io2=-1, ip1=-1, ip2=-1, im1=-1, im2=-1;
    diffmass_min = 10000, mass_min = 10000, pt_max2 = 0, pt_max = 0;
    for (unsigned int im=0; im<vSelectedJets.size(); im++){
        if (im==ib1 || im==ib2 || im==ik1 || im==ik2) continue;
        if (vSelectedJets.at(im).pt() > pt_max ) {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vSelectedJets.at(im).pt();
            im1 = im;
        }
        if (vSelectedJets.at(im).pt() < pt_max && vSelectedJets.at(im).pt() > pt_max2){
            pt_max2 = vSelectedJets.at(im).pt();
            im2 = im;
        }
        for (unsigned int in=0; in<vSelectedJets.size(); in++){
            if (in==ib1 || in==ib2 || in==ik1 || in==ik2 || in==im) continue;
            Pjet1.SetPtEtaPhiE(vSelectedJets.at(im).pt(), vSelectedJets.at(im).eta(), vSelectedJets.at(im).phi(), vSelectedJets.at(im).E());
            Pjet2.SetPtEtaPhiE(vSelectedJets.at(in).pt(), vSelectedJets.at(in).eta(), vSelectedJets.at(in).phi(), vSelectedJets.at(in).E());
            if (TMath::Abs((Pjet1+Pjet2).M()-80.419)<diffmass_min){
                io1=im;
                io2=in;
                diffmass_min = TMath::Abs((Pjet1+Pjet2).M()-80.419);
            }
            if ((Pjet1+Pjet2).M()<mass_min){
                ip1=im;
                ip2=in;
                mass_min = (Pjet1+Pjet2).M();
            }
        }
    }

    multilepton_JetHighestPt1_Id = -999;
    multilepton_JetHighestPt2_Id = -999;
    multilepton_JetClosestMw1_Id = -999;
    multilepton_JetClosestMw2_Id = -999;
    multilepton_JetLowestMjj1_Id = -999;
    multilepton_JetLowestMjj2_Id = -999;

    multilepton_JetHighestPt1_2ndPair_Id = -999;
    multilepton_JetHighestPt2_2ndPair_Id = -999;
    multilepton_JetClosestMw1_2ndPair_Id = -999;
    multilepton_JetClosestMw2_2ndPair_Id = -999;
    multilepton_JetLowestMjj1_2ndPair_Id = -999;
    multilepton_JetLowestMjj2_2ndPair_Id = -999;

    TLorentzVector Jet1, Jet2;

    if (ij1!=-1 && ij2==-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
    }
    if (ij1!=-1 && ij2!=-1) {
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
        FillJetInfoOutputTree(&multilepton_JetHighestPt1_Id, 1, &multilepton_JetHighestPt1_P4, Jet1, &multilepton_JetHighestPt1_CSV, vSelectedJets.at(ij1).CSVv2(), &multilepton_JetHighestPt1_JEC_Up, &multilepton_JetHighestPt1_JEC_Down, vSelectedJets.at(ij1).JES_uncert(), &multilepton_JetHighestPt1_JER_Up, &multilepton_JetHighestPt1_JER_Down, vSelectedJets.at(ij1).pt_JER(), vSelectedJets.at(ij1).pt_JER_up(), vSelectedJets.at(ij1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetHighestPt2_Id, 1, &multilepton_JetHighestPt2_P4, Jet2, &multilepton_JetHighestPt2_CSV, vSelectedJets.at(ij2).CSVv2(), &multilepton_JetHighestPt2_JEC_Up, &multilepton_JetHighestPt2_JEC_Down, vSelectedJets.at(ij2).JES_uncert(), &multilepton_JetHighestPt2_JER_Up, &multilepton_JetHighestPt2_JER_Down, vSelectedJets.at(ij2).pt_JER(), vSelectedJets.at(ij2).pt_JER_up(), vSelectedJets.at(ij2).pt_JER_down());
        //multilepton_JetHighestPt1_Id = 1;
        //multilepton_JetHighestPt2_Id = 1;
        //multilepton_JetHighestPt1_P4.SetPtEtaPhiE(vSelectedJets.at(ij1).pt(), vSelectedJets.at(ij1).eta(), vSelectedJets.at(ij1).phi(), vSelectedJets.at(ij1).E());
        //multilepton_JetHighestPt2_P4.SetPtEtaPhiE(vSelectedJets.at(ij2).pt(), vSelectedJets.at(ij2).eta(), vSelectedJets.at(ij2).phi(), vSelectedJets.at(ij2).E());
    }
    if (ik1!=-1 && ik2!=-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
        FillJetInfoOutputTree(&multilepton_JetClosestMw1_Id, 2, &multilepton_JetClosestMw1_P4, Jet1, &multilepton_JetClosestMw1_CSV, vSelectedJets.at(ik1).CSVv2(), &multilepton_JetClosestMw1_JEC_Up, &multilepton_JetClosestMw1_JEC_Down, vSelectedJets.at(ik1).JES_uncert(), &multilepton_JetClosestMw1_JER_Up, &multilepton_JetClosestMw1_JER_Down, vSelectedJets.at(ik1).pt_JER(), vSelectedJets.at(ik1).pt_JER_up(), vSelectedJets.at(ik1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetClosestMw2_Id, 2, &multilepton_JetClosestMw2_P4, Jet2, &multilepton_JetClosestMw2_CSV, vSelectedJets.at(ik2).CSVv2(), &multilepton_JetClosestMw2_JEC_Up, &multilepton_JetClosestMw2_JEC_Down, vSelectedJets.at(ik2).JES_uncert(), &multilepton_JetClosestMw2_JER_Up, &multilepton_JetClosestMw2_JER_Down, vSelectedJets.at(ik2).pt_JER(), vSelectedJets.at(ik2).pt_JER_up(), vSelectedJets.at(ik2).pt_JER_down());
        //multilepton_JetClosestMw1_Id = 2;
        //multilepton_JetClosestMw2_Id = 2;
        //multilepton_JetClosestMw1_P4.SetPtEtaPhiE(vSelectedJets.at(ik1).pt(), vSelectedJets.at(ik1).eta(), vSelectedJets.at(ik1).phi(), vSelectedJets.at(ik1).E());
        //multilepton_JetClosestMw2_P4.SetPtEtaPhiE(vSelectedJets.at(ik2).pt(), vSelectedJets.at(ik2).eta(), vSelectedJets.at(ik2).phi(), vSelectedJets.at(ik2).E());
    }
    if (il1!=-1 && il2!=-1){
        Jet1.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        Jet2.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
        FillJetInfoOutputTree(&multilepton_JetLowestMjj1_Id, 3, &multilepton_JetLowestMjj1_P4, Jet1, &multilepton_JetLowestMjj1_CSV, vSelectedJets.at(il1).CSVv2(), &multilepton_JetLowestMjj1_JEC_Up, &multilepton_JetLowestMjj1_JEC_Down, vSelectedJets.at(il1).JES_uncert(), &multilepton_JetLowestMjj1_JER_Up, &multilepton_JetLowestMjj1_JER_Down, vSelectedJets.at(il1).pt_JER(), vSelectedJets.at(il1).pt_JER_up(), vSelectedJets.at(il1).pt_JER_down());
        FillJetInfoOutputTree(&multilepton_JetLowestMjj2_Id, 3, &multilepton_JetLowestMjj2_P4, Jet2, &multilepton_JetLowestMjj2_CSV, vSelectedJets.at(il2).CSVv2(), &multilepton_JetLowestMjj2_JEC_Up, &multilepton_JetLowestMjj2_JEC_Down, vSelectedJets.at(il2).JES_uncert(), &multilepton_JetLowestMjj2_JER_Up, &multilepton_JetLowestMjj2_JER_Down, vSelectedJets.at(il2).pt_JER(), vSelectedJets.at(il2).pt_JER_up(), vSelectedJets.at(il2).pt_JER_down());
        //multilepton_JetLowestMjj1_Id = 3;
        //multilepton_JetLowestMjj2_Id = 3;
        //multilepton_JetLowestMjj1_P4.SetPtEtaPhiE(vSelectedJets.at(il1).pt(), vSelectedJets.at(il1).eta(), vSelectedJets.at(il1).phi(), vSelectedJets.at(il1).E());
        //multilepton_JetLowestMjj2_P4.SetPtEtaPhiE(vSelectedJets.at(il2).pt(), vSelectedJets.at(il2).eta(), vSelectedJets.at(il2).phi(), vSelectedJets.at(il2).E());
    }


    // ##########################################################################
    // #      _                  _       _                                      #
    // #  ___| |_ __ _ _ __   __| | __ _| | ___  _ __   ___    __ _  ___ _ __   #
    // # / __| __/ _` | '_ \ / _` |/ _` | |/ _ \| '_ \ / _ \  / _` |/ _ \ '_ \  #
    // # \__ \ || (_| | | | | (_| | (_| | | (_) | | | |  __/ | (_| |  __/ | | | #
    // # |___/\__\__,_|_| |_|\__,_|\__,_|_|\___/|_| |_|\___|  \__, |\___|_| |_| #
    // #                                                      |___/             #
    // #                                                                        #
    // ##########################################################################

    if( !_isdata )
    {
        for(unsigned int itruth = 0; itruth < vTruth->at(0).mc_truth_label().size() ; itruth++)
        {
            TLorentzVector LeptonX;

            if( vTruth->at(0).mc_truth_label().at(itruth) == 1 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth),
                                        vTruth->at(0).mc_truth_eta().at(itruth),
                                        vTruth->at(0).mc_truth_phi().at(itruth),
                                        vTruth->at(0).mc_truth_E().at(itruth) );

                multilepton_h0_P4 = LeptonX;
                multilepton_h0_Id = vTruth->at(0).mc_truth_id().at(itruth);
            }

            if( vTruth->at(0).mc_truth_label().at(itruth) == 2 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth),
                                        vTruth->at(0).mc_truth_eta().at(itruth),
                                        vTruth->at(0).mc_truth_phi().at(itruth),
                                        vTruth->at(0).mc_truth_E().at(itruth) );

                multilepton_t1_P4 = LeptonX;
                multilepton_t1_Id = vTruth->at(0).mc_truth_id().at(itruth);
            }

            if( vTruth->at(0).mc_truth_label().at(itruth) == 3 )
            {
                LeptonX.SetPtEtaPhiE(   vTruth->at(0).mc_truth_pt().at(itruth),
                                        vTruth->at(0).mc_truth_eta().at(itruth),
                                        vTruth->at(0).mc_truth_phi().at(itruth),
                                        vTruth->at(0).mc_truth_E().at(itruth) );

                multilepton_t2_P4 = LeptonX;
                multilepton_t2_Id = vTruth->at(0).mc_truth_id().at(itruth);
            }
        }
    }

    multilepton_mET.SetPtEtaPhiE(vEvent->at(0).metpt(), 0, vEvent->at(0).metphi(), vEvent->at(0).metpt());
    multilepton_mETcov00 = vEvent->at(0).metcov00();
    multilepton_mETcov01 = vEvent->at(0).metcov01();
    multilepton_mETcov10 = vEvent->at(0).metcov10();
    multilepton_mETcov11 = vEvent->at(0).metcov11();
    multilepton_mHT = vEvent->at(0).metsumet();

    mc_ttZhypAllowed = 0;
    /*
        if(vSelectedLeptons.size()==3)
        {
            if ( vSelectedLeptons.at(0).charge()==vSelectedLeptons.at(1).charge() && vSelectedLeptons.at(1).charge()==vSelectedLeptons.at(2).charge() ) mc_ttZhypAllowed =-1;
            else if (  ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(1).id() )
            || ( vSelectedLeptons.at(0).id() == -vSelectedLeptons.at(2).id() )
            || ( vSelectedLeptons.at(1).id() == -vSelectedLeptons.at(2).id() ))
            mc_ttZhypAllowed = 1;
        }
    */

    if (multilepton_Lepton1_Id!=-999 && multilepton_Lepton2_Id!=-999 && multilepton_Lepton3_Id!=-999){
        if (multilepton_Lepton1_Id*multilepton_Lepton2_Id>0 && multilepton_Lepton2_Id*multilepton_Lepton3_Id>0) mc_ttZhypAllowed =-1;
        else if ( (multilepton_Lepton1_Id==-multilepton_Lepton2_Id)
                || (multilepton_Lepton1_Id==-multilepton_Lepton3_Id)
                || (multilepton_Lepton2_Id==-multilepton_Lepton3_Id))
            mc_ttZhypAllowed = 1; }

    mc_nJets25 = vSelectedJets.size();

    tOutput->Fill();
}

void TTbarHiggsMultileptonAnalysis::FillJetInfoOutputTree(int* tree_Id, int Id, TLorentzVector* tree_P4, TLorentzVector P4, float* tree_CSV, float CSV, float* tree_JEC_Up, float* tree_JEC_Down, float JEC_value, float* tree_JER_Up, float* tree_JER_Down, float JER, float JER_Up, float JER_Down){

    *tree_Id = Id;
    *tree_P4 = P4;

    *tree_CSV = CSV;

    *tree_JEC_Up = P4.E()*(1.+JEC_value);
    *tree_JEC_Down = P4.E()*(1.-JEC_value);

    *tree_JER_Up = P4.E()*JER_Up/JER;
    *tree_JER_Down = P4.E()*JER_Down/JER;

    return;
}



void TTbarHiggsMultileptonAnalysis::selectBjets(std::string BjetSel, int* ibsel1, int* ibsel2, bool doSelectOnlyBjets){

    //Selects the two highest b-tag jets. If only one b-tag select just this one.
    int ib1=-1, ib2=-1;

    if (BjetSel=="HighestBtagDiscrim"){
        float btag_max=-999, btag_max2=-999;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            if (doSelectOnlyBjets && (vSelectedJets.at(ib).CSVv2()<0.5426)) continue;
            if (vSelectedJets.at(ib).CSVv2()>btag_max){
                btag_max2 = btag_max;
                ib2 = ib1;
                btag_max = vSelectedJets.at(ib).CSVv2();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).CSVv2()<btag_max && vSelectedJets.at(ib).CSVv2()>btag_max2){
                btag_max2 = vSelectedJets.at(ib).CSVv2();
                ib2 = ib;
            }
        }
    }
    if (BjetSel=="BtagHighestPt"){
        float pt_max=0, pt_max2=0;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            //            if (vSelectedJets.at(ib).CSVv2()<0.423) continue;
            if (vSelectedJets.at(ib).pt()>pt_max){
                pt_max2 = pt_max;
                ib2 = ib1;
                pt_max = vSelectedJets.at(ib).pt();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).pt()<pt_max && vSelectedJets.at(ib).pt()>pt_max2){
                pt_max2 = vSelectedJets.at(ib).pt();
                ib2 = ib;
            }
        }
    }

    *ibsel1 = ib1;
    *ibsel2 = ib2;

}

float TTbarHiggsMultileptonAnalysis::Phi_0_2Pi(float phi)
{
    float phi_0_2pi = phi;
    if (phi>= TMath::Pi()) phi_0_2pi  -= 2.*TMath::Pi();
    if (phi<0.)            phi_0_2pi  += 2.*TMath::Pi();
    return phi_0_2pi;
}

float TTbarHiggsMultileptonAnalysis::Phi_MPi_Pi(float phi)
{
    float phi_MPi_Pi = phi;
    while(phi_MPi_Pi > TMath::Pi()) {phi_MPi_Pi -= 2.*TMath::Pi();}
    while(phi_MPi_Pi < - TMath::Pi()) {phi_MPi_Pi += 2.*TMath::Pi();}

	return phi_MPi_Pi;
}

float TTbarHiggsMultileptonAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    if(fabs(phi2 - phi1) > pow(10, 3) ) {return -999;} //Some phi values are huge in vTruth (bad init ?)

    float DeltaPhi = Phi_MPi_Pi( phi2 - phi1 );
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}


/*
//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString TTbarHiggsMultileptonAnalysis::Convert_Number_To_TString(double number, int precision=10)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a double
double TTbarHiggsMultileptonAnalysis::Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}


//Use stat function (from library sys/stat) to check if a file exists
bool TTbarHiggsMultileptonAnalysis::Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}
*/
