//ThreeLeptonSelection_THQ3l
//TwoLeptonSelection_THQ2l

//FIXME -- modified leptons pT cuts !!!

#include "../include/tHqMultileponAnalysis.h"

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

//Cutflow vector --> Count events for THQ_3l_SR selection, after each cut
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

    InitFiles();
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
void TTbarHiggsMultileptonAnalysis::InitFiles()
{
    //--- Load MVA weight files
    Load_MVA();

    // Loading weight files and creating corrsponding histograms

    // b-tagging
    std::string inputFileHF = "/opt/sbg/scratch1/cms/TTH/weight/csv_rwt_fit_hf_76x_2016_02_08.root"; //-- needs update ??
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


//--- To be called for each event
void TTbarHiggsMultileptonAnalysis::InitVectors()
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
    vLightJets_FwdPtCut.clear();
}

//--- To be called in between each selection function
void TTbarHiggsMultileptonAnalysis::InitVariables()
{
	weightfake = 0;
	weightflip = 0;

	signal_3l_TT_MVA = -9; signal_3l_TTV_MVA = -9;
	signal_2l_TT_MVA = -9; signal_2l_TTV_MVA = -9;

	channel = -1;

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

    lep1Pt = 0; lep2Pt = 0; lep3Pt = 0; inv_mll = 0; hardestBjetPt = 0; hardestBjetEta = 0; fwdJetPt = 0;

    //Some categories (e.g. 3l SR & Training) overlap, whereas all the others are orthogonal
    //---> Make sure that when we fill the TTree, only 1 category is TRUE !
    is_3l_THQ_SR = false;
    is_3l_THQ_Training = false;
    is_3l_Z_CR = false;
    is_3l_AppFakes_SR = false;
    is_3l_GammaConv_SR = false;

    is_2l_THQ_SR = false;
    is_2l_THQ_Training = false;
    is_2l_EMU_CR = false;
    is_2l_AppFakes_SR = false;
    is_2l_GammaConv_SR = false;
    is_2l_QFlip_SR = false;

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
			weightfake = 1.;
			weightflip = 1.;
        }


        // ###########################################################
        // #  _       _ _   _       _ _           _   _              #
        // # (_)_ __ (_) |_(_) __ _| (_)___  __ _| |_(_) ___  _ __   #
        // # | | '_ \| | __| |/ _` | | / __|/ _` | __| |/ _ \| '_ \  #
        // # | | | | | | |_| | (_| | | \__ \ (_| | |_| | (_) | | | | #
        // # |_|_| |_|_|\__|_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_| #
        // #                                                         #
        // ###########################################################

        InitVectors();


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
            if      ( EM  &&                               (emdataset) ) is_trigger = true;
            else if ( !EM && MM  &&                        (mmdataset) ) is_trigger = true;
            else if ( !EM && !MM && EE  &&                 (eedataset) ) is_trigger = true;
            else if ( !EM && !MM && !EE && M  &&           (mdataset ) ) is_trigger = true;
            else if ( !EM && !MM && !EE && !M && E &&      (edataset ) ) is_trigger = true;

            //The event will be rejected in sel. func., but counted
            if(!is_trigger)
            {
                weight = 0;
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

            if(vMuon->at(imuon).isLooseTTH() )
            {
                vLooseLeptons.push_back(l);
            }

            //FIXME -- tmp lepton pT cut
            if(vMuon->at(imuon).pt() < 15) {continue;}

            if(vMuon->at(imuon).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }
            if(vMuon->at(imuon).isTightTTH() )
            {
                vTightLeptons.push_back(l);
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

            if(vElectron->at(ielectron).isLooseTTH() )
            {
                vLooseLeptons.push_back(l);
            }

            //FIXME -- tmp lepton pT cut
            if(vElectron->at(ielectron).pt() < 15) {continue;}

            if(vElectron->at(ielectron).isFakeableTTH() )
            {
                vSelectedLeptons.push_back(l);
                vFakeableLeptons.push_back(l);
            }
            if(vElectron->at(ielectron).isTightTTH() )
            {
                if(!vElectron->at(ielectron).cutEventSel() || !vElectron->at(ielectron).noLostHits() ) {continue;} //conversion veto
                vTightLeptons.push_back(l);
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

        nTightLep = vTightLeptons.size();
        nFakeableLep = vFakeableLeptons.size();

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
        nLightJets = 0;
        nLightJets_Fwd40 = 0;

        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {
            if(vJet->at(ijet).pt() < 25) {continue;}

			//---------
			//All jets with pT>25 in vJets are "selected jets"
			vSelectedJets.push_back(vJet->at(ijet));

            //--- B-tagged jets
            if(vJet->at(ijet).CSVv2() > 0.5426 && fabs(vJet->at(ijet).eta()) < 2.4)
            {
                vLooseBTagJets.push_back(vJet->at(ijet)); nLooseBJets++;

                if(vJet->at(ijet).CSVv2() > 0.8484) {nMediumBJets++;}
            }
            //--- Light jets
            else if(vJet->at(ijet).CSVv2() < 0.5426)
            {
                vLightJets.push_back(vJet->at(ijet)); nLightJets++;

                if(fabs(vJet->at(ijet).eta() ) < 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJet->at(ijet)); nLightJets_Fwd40++;
                }
                else if(vJet->at(ijet).pt() > 40 && fabs(vJet->at(ijet).eta() ) > 2.4)
                {
                    vLightJets_FwdPtCut.push_back(vJet->at(ijet)); nLightJets_Fwd40++;
                }
            }
        }



        //--------------------------------------------
        // #####  ###### #####  #    #  ####     #####  #####  # #    # #####
        // #    # #      #    # #    # #    #    #    # #    # # ##   #   #
        // #    # #####  #####  #    # #         #    # #    # # # #  #   #
        // #    # #      #    # #    # #  ###    #####  #####  # #  # #   #
        // #    # #      #    # #    # #    #    #      #   #  # #   ##   #
        // #####  ###### #####   ####   ####     #      #    # # #    #   #
        //--------------------------------------------

        bool debug_printout = false;
        if(debug_printout)
        {
            cout<<endl<<endl<<BOLD(FBLU("------------ EVENT -----------"))<<endl;

            cout<<FYEL("--- Loose Tagged Jets : ")<<endl;
            for(int ijet=0; ijet<vLooseBTagJets.size(); ijet++)
            {
                cout<<ijet<<" pT = "<<vLooseBTagJets.at(ijet).pt()<<" / eta = "<<vLooseBTagJets.at(ijet).eta()<<" / phi = "<<vLooseBTagJets.at(ijet).phi()<<" / CSV = "<<vLooseBTagJets.at(ijet).CSVv2()<<endl;
            }

            cout<<FYEL("--- Light Jets : ")<<endl;
            for(int ijet=0; ijet<vLightJets.size(); ijet++)
            {
                cout<<ijet<<" pT = "<<vLightJets.at(ijet).pt()<<" / eta = "<<vLightJets.at(ijet).eta()<<" / phi = "<<vLightJets.at(ijet).phi()<<" / CSV = "<<vLightJets.at(ijet).CSVv2()<<endl;
            }

            cout<<FYEL("--- Fakeable Leptons : ")<<endl;
            for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
            {
                cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
            }
        }

        // ############################################
        // #           _           _   _              #
        // #  ___  ___| | ___  ___| |_(_) ___  _ __   #
        // # / __|/ _ \ |/ _ \/ __| __| |/ _ \| '_ \  #
        // # \__ \  __/ |  __/ (__| |_| | (_) | | | | #
        // # |___/\___|_|\___|\___|\__|_|\___/|_| |_| #
        // #                                          #
        // ############################################

        // if(Sample_isUsed_forTraining() == true) //Only few samples are used for training //FIXME
        {
            ThreeLeptonSelection_THQ3l_Training(jentry);
            TwoLeptonSelection_THQ2l_Training(jentry);
            if(_sampleName.Contains("ttWJets_13TeV_madgraphMLM") || _sampleName.Contains("ttZJets_13TeV_madgraphMLM") ) {continue;} //Used *only* for training
        }

        if( _isdata ) //Apply to datasets only
        {
            ThreeLeptonSelection_ApplicationFakes(jentry);
            TwoLeptonSelection_ApplicationFakes(jentry);

            TwoLeptonSelection_ApplicationChargeFlip(jentry);
        }

        ThreeLeptonSelection_THQ3l_SR(jentry);
        ThreeLeptonSelection_THQ3l_Z_CR(jentry);

        TwoLeptonSelection_THQ2l_SR(jentry);
        TwoLeptonSelection_THQ2l_EMU_OS_CR(jentry); // -- disactivate it for now : only TTbar + ST, >1M events !
    }


	ofstream file_out("cutflow.txt");
    for(int i=0; i<v_cutflow.size(); i++)
    {
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


//---------------------------------------------------------------------------

void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_SR(int evt)
{
    InitVariables();

	v_cutflow[0]++;

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}
    v_cutflow[1]++;

    // ####################
    // # Common selection #
    // ####################
	//At least 3 FO leptons

    if(vFakeableLeptons.size() < 3) {return;}
    v_cutflow[2]++;

    //Ask exactly 3 tight leptons, which must be the 3 hardest FO leptons
	if(vTightLeptons.size() != 3) {return;}
	if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH() && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}
    v_cutflow[3]++;


    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15 || vFakeableLeptons.at(2).pt() < 15) {return;}

    if(vSelectedJets.size() < 2) {return;}
    v_cutflow[4]++;

    if(nMediumBJets == 0) {return;}
    v_cutflow[5]++;

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}
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
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && vFakeableLeptons.at(i).charge() == -vFakeableLeptons.at(j).charge() )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}
    v_cutflow[8]++;


    is_3l_THQ_SR = 1; //3l SR category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 39) {channel = 0;} //uuu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 37) {channel = 1;} //uue
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 35) {channel = 2;} //eeu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 33) {channel = 3;} //eee

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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
    }
    // if( fabs(sum_charges_3l) != 1 ) return; //Cf. AN : triple charge consistency not required in 3l


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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {

                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(2).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt(); lep3Pt = vFakeableLeptons.at(2).pt();
    inv_mll = mll_tmp;
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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






void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_Training(int evt)
{
    InitVariables();

    if(weight==0) {return;}

    if(vFakeableLeptons.size() < 3) {return;}

    //-- Only consider the 3 hardest 'FakeableObject' leptons
    if(vFakeableLeptons.at(0).pt() < 20 || vFakeableLeptons.at(1).pt() < 10 || vFakeableLeptons.at(2).pt() < 10)
    {
        return;
    }

	int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
    	if( fabs(vSelectedJets.at(ijet).eta()) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nLooseBJets == 0) {return;}

    if(vLightJets.size() == 0) {return;}

    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && vFakeableLeptons.at(i).charge() == -vFakeableLeptons.at(j).charge() )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_3l_THQ_Training = 1;


// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 39) {channel = 0;} //uuu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 37) {channel = 1;} //uue
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 35) {channel = 2;} //eeu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 33) {channel = 3;} //eee


    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;

    TLorentzVector jetp4;
    for(int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiE(vSelectedJets.at(i).pt(), vSelectedJets.at(i).eta(), vSelectedJets.at(i).phi(), vSelectedJets.at(i).E());
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(2).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt(); lep3Pt = vFakeableLeptons.at(2).pt();
    inv_mll = mll_tmp;
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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

        cout<<FYEL("--- FO Leptons : ")<<endl;
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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







void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_THQ3l_Z_CR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################

    //At least 3 leptons
    if(vFakeableLeptons.size() < 3) {return;}

    //Ask exactly 3 tight leptons, which must be the 3 hardest FO leptons
	if(vTightLeptons.size() != 3) {return;}
	if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH() && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}

    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15 || vFakeableLeptons.at(2).pt() < 15) {return;}

    if( fabs(vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge() + vFakeableLeptons.at(2).charge() ) != 1) {return;}

    if(vSelectedJets.size() < 2) {return;}

    if(nLooseBJets == 0) {return;}

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}


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


    // ##########
    // # OSSF pair - Z INVERTED veto (selection enriched in ttZ/tZ, WZ events)
    // ##########
    bool pass_inverted_Zveto = false;
    float mll_tmp = 0; float best_mll=999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && vFakeableLeptons.at(i).charge() == -vFakeableLeptons.at(j).charge() )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if(fabs(mll_tmp - 91.188) < 15)
                {
                    pass_inverted_Zveto = true ;
                    if(fabs(mll_tmp - 91.188) < best_mll) {best_mll = mll_tmp;}
                }
            }

        }
    }
    if(!pass_inverted_Zveto) {return;}


    is_3l_Z_CR = 1; //3l SR category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 39) {channel = 0;} //uuu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 37) {channel = 1;} //uue
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 35) {channel = 2;} //eeu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 33) {channel = 3;} //eee


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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;

    // if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return; //CHANGED

    if(DEBUG) std::cout << "nJets and met_ld Ok... ";

    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(2).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt(); lep3Pt = vFakeableLeptons.at(2).pt();
    inv_mll = best_mll;
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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







void TTbarHiggsMultileptonAnalysis::TwoLeptonSelection_THQ2l_SR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 leptons
    if(vFakeableLeptons.size() < 2) {return;}

    //==2 tight leptons, which must be the 2 hardest FO leptons
	if(vTightLeptons.size() != 2) {return;}
    if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) ) {return;}


    //SS lepton pair
    if(vFakeableLeptons.at(0).charge() * vFakeableLeptons.at(1).charge() < 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vFakeableLeptons.at(0).passTightCharge() ) {return;}
    if(!vFakeableLeptons.at(1).passTightCharge() ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}


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


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && fabs(vFakeableLeptons.at(i).id())==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_2l_THQ_SR = 1; //2l SR category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 26) {channel = 0;} //uu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 24) {channel = 1;} //ue+eu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 22) {channel = 2;} //ee


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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }


    int sum_charges_3l = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();
    for(int i=0; i<2; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();

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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(1).pt();

    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt();
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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
    signal_2l_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2l_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}






void TTbarHiggsMultileptonAnalysis::TwoLeptonSelection_THQ2l_Training(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}


    //At most 2 tight leptons
	if(vTightLeptons.size() > 2) {return;}

    //SS lepton pair
    if(vFakeableLeptons.at(0).charge() * vFakeableLeptons.at(1).charge() < 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    // if(!vFakeableLeptons.at(0).passTightCharge() ) {return;}
    // if(!vFakeableLeptons.at(1).passTightCharge() ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt() < 20 || vFakeableLeptons.at(1).pt() < 10) {return;}

    //Jet cuts
    int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if( fabs(vSelectedJets.at(ijet).eta()) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nLooseBJets == 0) {return;}

    if(vLightJets.size() == 0)  {return;}

    is_2l_THQ_Training = 1; //2l Training category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 26) {channel = 0;} //uu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 24) {channel = 1;} //ue+eu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 22) {channel = 2;} //ee

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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;


    int sum_charges_3l = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();
    for(int i=0; i<2; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();

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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(1).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt();
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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
    signal_2l_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2l_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}








//--- Enriched in dileptonic ttbar events
//-- E/mu 2l SS control region
void TTbarHiggsMultileptonAnalysis::TwoLeptonSelection_THQ2l_EMU_OS_CR(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 leptons
    if(vFakeableLeptons.size() < 2) {return;}


    //==2 tight leptons, which must be the 2 hardest leptons
	if(vTightLeptons.size() != 2) {return;}
    if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) ) {return;}

    //OS pair
    if(vFakeableLeptons.at(0).charge() * vFakeableLeptons.at(1).charge() > 0) {return;}

    //e-mu pair
    if( fabs(vFakeableLeptons.at(0).id() * vFakeableLeptons.at(1).id()) != 143 ) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vFakeableLeptons.at(0).passTightCharge() ) {return;}
    if(!vFakeableLeptons.at(1).passTightCharge() ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 10) {return;}
    if( fabs(vFakeableLeptons.at(1).id()) == 11 && vFakeableLeptons.at(1).pt() < 15) {return;}

    //Jet cuts
    int jets_tmp = 0;
    for(int ijet=0; ijet<vSelectedJets.size(); ijet++)
    {
        if( fabs(vSelectedJets.at(ijet).eta()) < 2.4 ) jets_tmp++;
    }
    if(jets_tmp < 2) {return;} //At least 2 jets with pT>25 & eta<2.4

    if(nMediumBJets == 0) {return;}

	// if(vLightJets_FwdPtCut.size() == 0) {return;}
    if(vLightJets.size() == 0) {return;}


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

    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && fabs(vFakeableLeptons.at(i).id())==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_2l_EMU_CR = 1; //2l SR category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 26) {channel = 0;} //uu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 24) {channel = 1;} //ue+eu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 22) {channel = 2;} //ee

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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }

    if(DEBUG) std::cout << std::endl << "nJets: " << vSelectedJets.size() << " met_ld: " << met_ld << " isSFOS " << isSFOS << std::endl;


    int sum_charges_3l = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();
    for(int i=0; i<2; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();

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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }


    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(1).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt();
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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
    signal_2l_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2l_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}







void TTbarHiggsMultileptonAnalysis::ThreeLeptonSelection_ApplicationFakes(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}

    // ####################
    // # Common selection #
    // ####################
	//At least 3 FO leptons

    if(vFakeableLeptons.size() < 3) {return;}

    if(vTightLeptons.size() > 3) {return;}

	// if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) || !(vFakeableLeptons.at(2).isTightTTH() && vFakeableLeptons.at(2).cutEventSel() && vFakeableLeptons.at(2).noLostHits()) ) {return;}

    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15 || vFakeableLeptons.at(2).pt() < 15) {return;}

    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}

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


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && vFakeableLeptons.at(i).charge() == -vFakeableLeptons.at(j).charge() )
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 15) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_3l_AppFakes_SR = 1; //3l AppFakes SR category

//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 39) {channel = 0;} //uuu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 37) {channel = 1;} //uue
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 35) {channel = 2;} //eeu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id())+fabs(vFakeableLeptons.at(2).id()) == 33) {channel = 3;} //eee


// ##################################################################################################################################

    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    //only considering the 3 hardest FO leptons
    for(int i=0; i<3; i++)
    {
        if( !vFakeableLeptons.at(i).isTightTTH() ) //For each lepton failing the tight requirements, multiply event by a weight
        {
            //std::cout << "Applying fake rate" << std::endl;
            //std::cout << "i: " << i << std::endl;
            //std::cout << "pt: " << vFakeableLeptons.at(i).pt() << std::endl;
            //std::cout << "fake? " << vFakeableLeptons.at(i).isFakeableTTH() << std::endl;
            leptonsPts.push_back(     vFakeableLeptons.at(i).pt() );
            leptonsEtas.push_back(    vFakeableLeptons.at(i).eta() );
            leptonsIds.push_back(     vFakeableLeptons.at(i).id() );
        }
    }

    weightfake = get_FR_weight(leptonsPts, leptonsEtas, leptonsIds);

    // cout<<"weightfake = "<<weightfake<<endl;


// ##################################################################################################################################


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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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

    // if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) return;


    int sum_charges_3l = 0;
    for(int i=0; i<3; i++)
    {
        sum_charges_3l = sum_charges_3l + vFakeableLeptons.at(i).charge();
    }

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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {

                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(2).pt();


    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt(); lep3Pt = vFakeableLeptons.at(2).pt();
    inv_mll = mll_tmp;
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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







void TTbarHiggsMultileptonAnalysis::TwoLeptonSelection_ApplicationFakes(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}
    if(vTightLeptons.size() >= 2) {return;}

    // if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) ) {return;}


    //SS lepton pair
    if(vFakeableLeptons.at(0).charge() * vFakeableLeptons.at(1).charge() < 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vFakeableLeptons.at(0).passTightCharge() ) {return;}
    if(!vFakeableLeptons.at(1).passTightCharge() ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}


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


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && fabs(vFakeableLeptons.at(i).id())==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_2l_AppFakes_SR = 1; //2l SR category

// ##################################################################################################################################

    // #########################
    // # fake rate reweighting #
    // #########################

    std::vector<double> leptonsPts;
    std::vector<double> leptonsEtas;
    std::vector<int>    leptonsIds;

    //Only consider the 2 hardest FO leptons
    for(int i=0; i<2; i++)
    {
        if( !vFakeableLeptons.at(i).isTightTTH() ) //For each lepton failing the tight requirement, apply a weight to event
        {
            //std::cout << "Applying fake rate" << std::endl;
            //std::cout << "i: " << i << std::endl;
            //std::cout << "pt: " << vFakeableLeptons.at(i).pt() << std::endl;
            //std::cout << "fake? " << vFakeableLeptons.at(i).isFakeableTTH() << std::endl;
            leptonsPts.push_back(     vFakeableLeptons.at(i).pt()                );
            leptonsEtas.push_back(    vFakeableLeptons.at(i).eta()               );
            leptonsIds.push_back(     vFakeableLeptons.at(i).id()                );
        }
    }

    weightfake = get_FR_weight(leptonsPts, leptonsEtas, leptonsIds);

    // cout<<"weightfake = "<<weightfake<<endl;



//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 26) {channel = 0;} //uu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 24) {channel = 1;} //ue+eu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 22) {channel = 2;} //ee


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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }


    int sum_charges_3l = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();
    for(int i=0; i<2; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();

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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(1).pt();

    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt();
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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
    signal_2l_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2l_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}







void TTbarHiggsMultileptonAnalysis::TwoLeptonSelection_ApplicationChargeFlip(int evt)
{
    InitVariables();

    //If data doesn't pass trigger, weight=0
    if(weight==0) {return;}


    // ####################
    // # Common selection #
    // ####################

    //At least 2 FO leptons
    if(vFakeableLeptons.size() < 2) {return;}
    if(vTightLeptons.size() >= 2) {return;}

    // if( !(vFakeableLeptons.at(0).isTightTTH() && vFakeableLeptons.at(0).cutEventSel() && vFakeableLeptons.at(0).noLostHits()) || !(vFakeableLeptons.at(1).isTightTTH() && vFakeableLeptons.at(1).cutEventSel() && vFakeableLeptons.at(1).noLostHits()) ) {return;}


    //OS lepton pair
    if(vFakeableLeptons.at(0).charge() * vFakeableLeptons.at(1).charge() > 0) {return;}

    //Dilepton channel  : leptons need to pass tightCharge cut (dPt/pT>0.2 for muons, tightCharge>1 for ele)
    if(!vFakeableLeptons.at(0).passTightCharge() ) {return;}
    if(!vFakeableLeptons.at(1).passTightCharge() ) {return;}

    //pT cuts
    if(vFakeableLeptons.at(0).pt() < 25 || vFakeableLeptons.at(1).pt() < 15) {return;}

    //Jet cuts
    if(vSelectedJets.size() < 2) {return;}

    if(nMediumBJets == 0) {return;}

    // if(vLightJets_FwdPtCut.size() == 0)  {return;}
    if(vLightJets.size() == 0)  {return;}


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


    // ##########
    // # OSSF pair - Z veto #
    // ##########
    bool pass_Zveto = true;
    float mll_tmp = 0;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            if( fabs(vFakeableLeptons.at(i).id()) == fabs(vFakeableLeptons.at(j).id()) && fabs(vFakeableLeptons.at(i).id())==11)
            {
                mll_tmp = (vFakeableLeptons.at(i).p4() + vFakeableLeptons.at(j).p4() ).M();
                if( fabs(mll_tmp - 91.188) < 10) {pass_Zveto = false ;}
            }
        }
    }
    if(!pass_Zveto) {return;}


    is_2l_QFlip_SR = 1; //2l SR category

// ##################################################################################################################################

    // #########################
    // # QFlip reweighting     #
    // #########################

    weightflip = get_QF_weight(vFakeableLeptons.at(0).pt(), vFakeableLeptons.at(0).eta(), vFakeableLeptons.at(0).id(), vFakeableLeptons.at(1).pt(), vFakeableLeptons.at(1).eta(), vFakeableLeptons.at(1).id());

    // cout<<"weightflip = "<<weightflip<<endl;
	// weight*= weightflip;



//---------------------------------------------------------------------------
// ##################################################################################################################################

    //--- Determine leptonic channel of event
    if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 26) {channel = 0;} //uu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 24) {channel = 1;} //ue+eu
    else if(fabs(vFakeableLeptons.at(0).id())+fabs(vFakeableLeptons.at(1).id()) == 22) {channel = 2;} //ee


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

    for(int i=0; i<vFakeableLeptons.size(); i++)
    {
        lepton_px = lepton_px + vFakeableLeptons.at(i).p4().Px();
        lepton_py = lepton_py + vFakeableLeptons.at(i).p4().Py();
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
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<vFakeableLeptons.size(); j++)
        {
            if (  ( i                           != j                            )
                    && ( vFakeableLeptons.at(i).id() == -vFakeableLeptons.at(j).id() ) )
            { isSFOS = true ;}
        }
    }


    int sum_charges_3l = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();
    for(int i=0; i<2; i++)
    {
        sum_charges_3l += vFakeableLeptons.at(i).charge();
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
            pt2 = pt1;
            ijet_2nd_hardest_btag = ijet_hardest_btag;
            pt1 = fabs(vLooseBTagJets.at(ijet).pt());
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
    lepCharge = vFakeableLeptons.at(0).charge() + vFakeableLeptons.at(1).charge();

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
    for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
    {
        if( fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ) < tmp) {dEtaFwdJetClosestLep = fabs(vLightJets.at(ijet_forward).eta() - vFakeableLeptons.at(ilep).eta() ); tmp = dEtaFwdJetClosestLep;}
    }

    //--- Var8 : dPhi of highest pT SS lepton pair
    dPhiHighestPtSSPair = 0; tmp = -999;
    TLorentzVector lepi, lepj;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(vFakeableLeptons.at(i).charge()==vFakeableLeptons.at(j).charge() && (lepi+lepj).Pt() > tmp)
            {
                tmp = (lepi+lepj).Pt();
                dPhiHighestPtSSPair = fabs(Phi_MPi_Pi(vFakeableLeptons.at(i).phi() - vFakeableLeptons.at(j).phi() ) );
            }
        }
    }

    //--- Var9 : min. dR between any 2 leptons
    minDRll = 999;
    for(int i=0; i<vFakeableLeptons.size()-1; i++)
    {
        lepi.SetPtEtaPhiE(vFakeableLeptons.at(i).pt(), vFakeableLeptons.at(i).eta(), vFakeableLeptons.at(i).phi(), vFakeableLeptons.at(i).E() );
        for(int j=i+1; j<vFakeableLeptons.size(); j++)
        {
            lepj.SetPtEtaPhiE(vFakeableLeptons.at(j).pt(), vFakeableLeptons.at(j).eta(), vFakeableLeptons.at(j).phi(), vFakeableLeptons.at(j).E() );
            if(lepi.DeltaR(lepj) < minDRll)
            {
                minDRll = lepi.DeltaR(lepj);
            }
        }
    }

    //--- Var10 : pT of 3rd hardest lepton
    Lep3Pt = vFakeableLeptons.at(1).pt();

    //--- Fill additionnal variables, used for control only
    lep1Pt = vFakeableLeptons.at(0).pt(); lep2Pt = vFakeableLeptons.at(1).pt();
    hardestBjetPt = vLooseBTagJets.at(ijet_hardest_btag).pt(); hardestBjetEta = vLooseBTagJets.at(ijet_hardest_btag).eta();
    fwdJetPt = vLightJets.at(ijet_forward).pt();

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
        for(int ilep=0; ilep<vFakeableLeptons.size(); ilep++)
        {
            cout<<ilep<<" pT = "<<vFakeableLeptons.at(ilep).pt()<<" / eta = "<<vFakeableLeptons.at(ilep).eta()<<" / phi = "<<vFakeableLeptons.at(ilep).phi()<<" : Charge = "<<vFakeableLeptons.at(ilep).charge()<<endl;
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
    signal_2l_TT_MVA   = mva_2lss_tt->EvaluateMVA("BDTG method");
    signal_2l_TTV_MVA   = mva_2lss_ttV->EvaluateMVA("BDTG method");
    // ======================================================================================================


    fillOutputTree();
    // cout<<BOLD(FYEL("TREE FILLED !"))<<endl;
}











void TTbarHiggsMultileptonAnalysis::initializeOutputTree()
{
    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

	//-- Event main infos
	tOutput->Branch("channel",&channel,"channel/F");
    tOutput->Branch("weight",&weight,"weight/F");
	tOutput->Branch("is_trigger",&is_trigger,"is_trigger/B");
	tOutput->Branch("weightfake",&weightfake,"weightfake/F");
	tOutput->Branch("weightflip",&weightflip,"weightflip/F");
	tOutput->Branch("mc_event",&mc_event,"mc_event/I");
	tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F");


	//--- Categories & MVA
	tOutput->Branch("is_3l_THQ_SR",&is_3l_THQ_SR,"is_3l_THQ_SR/B");
	tOutput->Branch("is_3l_THQ_Training",&is_3l_THQ_Training,"is_3l_THQ_Training/B");
    tOutput->Branch("is_3l_Z_CR",&is_3l_Z_CR,"is_3l_Z_CR/B");
    tOutput->Branch("is_3l_AppFakes_SR",&is_3l_AppFakes_SR,"is_3l_AppFakes_SR/B");
    tOutput->Branch("is_3l_GammaConv_SR",&is_3l_GammaConv_SR,"is_3l_GammaConv_SR/B");

    tOutput->Branch("is_2l_THQ_SR",&is_2l_THQ_SR,"is_2l_THQ_SR/B");
    tOutput->Branch("is_2l_THQ_Training",&is_2l_THQ_Training,"is_2l_THQ_Training/B");
	tOutput->Branch("is_2l_EMU_CR",&is_2l_EMU_CR,"is_2l_EMU_CR/B");
    tOutput->Branch("is_2l_AppFakes_SR",&is_2l_AppFakes_SR,"is_2l_AppFakes_SR/B");
    tOutput->Branch("is_2l_GammaConv_SR",&is_2l_GammaConv_SR,"is_2l_GammaConv_SR/B");
    tOutput->Branch("is_2l_QFlip_SR",&is_2l_QFlip_SR,"is_2l_QFlip_SR/B");

	tOutput->Branch("signal_3l_TT_MVA",&signal_3l_TT_MVA,"signal_3l_TT_MVA/F");
	tOutput->Branch("signal_3l_TTV_MVA",&signal_3l_TTV_MVA,"signal_3l_TTV_MVA/F");
    tOutput->Branch("signal_2l_TT_MVA",&signal_2l_TT_MVA,"signal_2l_TT_MVA/F");
    tOutput->Branch("signal_2l_TTV_MVA",&signal_2l_TTV_MVA,"signal_2l_TTV_MVA/F");

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

	//-- More control vars
    tOutput->Branch("inv_mll",&inv_mll,"inv_mll/F");
    tOutput->Branch("hardestBjetPt",&hardestBjetPt,"hardestBjetPt/F");
    tOutput->Branch("hardestBjetEta",&hardestBjetEta,"hardestBjetEta/F");
    tOutput->Branch("fwdJetPt",&fwdJetPt,"fwdJetPt/F");
    tOutput->Branch("lep1Pt",&lep1Pt,"lep1Pt/F");
    tOutput->Branch("lep2Pt",&lep2Pt,"lep2Pt/F");
    tOutput->Branch("lep3Pt",&lep3Pt,"lep3Pt/F");

	//-- Cut variables
	tOutput->Branch("nLooseBJets",&nLooseBJets,"nLooseBJets/F");
	tOutput->Branch("nMediumBJets",&nMediumBJets,"nMediumBJets/F");
	tOutput->Branch("nLightJets",&nLightJets,"nLightJets/F");
	tOutput->Branch("nTightLep",&nTightLep,"nTightLep/F");
	tOutput->Branch("nFakeableLep",&nFakeableLep,"nFakeableLep/F");
	tOutput->Branch("nLightJets_Fwd40",&nLightJets_Fwd40,"nLightJets_Fwd40/F");


	//-- Other weights
    tOutput->Branch("weight_scale_muF0p5",&weight_scale_muF0p5,"weight_scale_muF0p5/F");
    tOutput->Branch("weight_scale_muF2",&weight_scale_muF2,"weight_scale_muF2/F");
    tOutput->Branch("weight_scale_muR0p5",&weight_scale_muR0p5,"weight_scale_muR0p5/F");
    tOutput->Branch("weight_scale_muR2",&weight_scale_muR2,"weight_scale_muR2/F");
    tOutput->Branch("weight_csv_down",&weight_csv_down,"weight_csv_down/F");
    tOutput->Branch("weight_csv_up",&weight_csv_up,"weight_csv_up/F");
    tOutput->Branch("weights_pdf","std::vector<float>",&weights_pdf);
    tOutput->Branch("ids_pdf","std::vector<std::string>",&ids_pdf);

	//-- Other, MEM necessary vars, ...
    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");
    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");

    tOutput->Branch("catJets",&catJets,"catJets/I");

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

    return;
}



void TTbarHiggsMultileptonAnalysis::fillOutputTree()
{
    int tot_charge = 0;
    int tot_id = 0;
    if (vFakeableLeptons.size()>=4) {
        for (unsigned int i=0; i<4; i++) {
            tot_charge += vFakeableLeptons.at(i).charge();
            tot_id += vFakeableLeptons.at(i).id();
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

    if (vFakeableLeptons.size()>=2){
        multilepton_Lepton1_P4 = vFakeableLeptons.at(0).p4();
        multilepton_Lepton1_Id = vFakeableLeptons.at(0).id();
        multilepton_Lepton2_P4 = vFakeableLeptons.at(1).p4();
        multilepton_Lepton2_Id = vFakeableLeptons.at(1).id();
    }

    if (vFakeableLeptons.size()>=3)
    {
        multilepton_Lepton3_P4 = vFakeableLeptons.at(2).p4();
        multilepton_Lepton3_Id = vFakeableLeptons.at(2).id();
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
                lep1_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vFakeableLeptons.at(0).eta(), vFakeableLeptons.at(0).phi() );
                if( lep1_dr_gen < lep1_dr_gen_min)
                {
                    lep1_dr_gen_min = lep1_dr_gen;  lep1_matched = itruth;
                }

                lep2_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vFakeableLeptons.at(1).eta(), vFakeableLeptons.at(1).phi() );
                if( lep2_dr_gen < lep2_dr_gen_min)
                {   lep2_dr_gen_min = lep2_dr_gen;  lep2_matched = itruth;  }

                if(vFakeableLeptons.size()>=3)
                {
                    lep3_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vFakeableLeptons.at(2).eta(), vFakeableLeptons.at(2).phi() );
                    if( lep3_dr_gen < lep3_dr_gen_min)
                    {   lep3_dr_gen_min = lep3_dr_gen;  lep3_matched = itruth;  }
                }

                if(vFakeableLeptons.size()>=4)
                {
                    lep4_dr_gen = GetDeltaR(vTruth->at(0).mc_truth_eta().at(itruth),  vTruth->at(0).mc_truth_phi().at(itruth), vFakeableLeptons.at(3).eta(), vFakeableLeptons.at(3).phi() );
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
            if(vFakeableLeptons.size()>=3)
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
            if(vFakeableLeptons.size()>=4)
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

//Decides if given sample is to be used for training or not (else, don't call training selection function)
bool TTbarHiggsMultileptonAnalysis::Sample_isUsed_forTraining()
{
    if(_sampleName.Contains("ttZJets_13TeV_madgraphMLM") || _sampleName.Contains("ttWJets_13TeV_madgraphMLM") || _sampleName.Contains("THQ") || _sampleName.Contains("TTJets") )
    {
        return true;
    }

    return false;
}

//Decides if given sample is to be used for build gamma-conversion sample or not
bool TTbarHiggsMultileptonAnalysis::Sample_isUsed_forGammaConv()
{
    if(_sampleName.Contains("TTGJets") || _sampleName.Contains("WGToLNuG") || _sampleName.Contains("ZGTo2LG") || _sampleName.Contains("TGJets") )
    {
        return true;
    }

    return false;
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
