#ifndef MULTILEPTONTREE_H
#define MULTILEPTONTREE_H

const string EtaRangeLabel[] = {"00eta08", "08eta16", "16eta24"};
const string EnergyRangeLabel[] = {"25E50", "50E80", "80E120", "120E200", "200E300", "EGT300"};
const float EtaRange[] = {0.0, 0.8, 1.6, 2.4};
const float EnergyRange[] = {25, 50, 80, 120, 200, 300, 7000};

const string MetRangeLabel[] = {"mET0E100Sum0E1200", "mET0E100Sum1200E1600", "mET0E100SumEGT1600", "mETGT100Sum0E1200", "mETGT100Sum1200E1600", "mETGT100SumEGT1600"};
const float MetRange[] = {0, 100, 7000};
const float MetSumRange[] = {0, 1200, 1600, 7000};

class MultiLeptonTree {

  public:

  MultiLeptonTree();
  ~MultiLeptonTree();

  void initializeOutputTree();
  void fillOutputTree();
  void FillJetInfoOutputTree(int* tree_Id, int id, TLorentzVector* tree_P4, Jet* jet);
  //void FillJetInfoOutputTree(int*, int, TLorentzVector*, TLorentzVector, float*, float, float*, float*, float, float*, float*, float, float, float);
 
  TTree* tOutput;

	Float_t weight;
        Float_t mc_weight;                                                                           // weight MC
        Float_t weight_PV;                                                                           // PU reweighting from PV distribution

	Int_t mc_ttZhypAllowed;
	Int_t catJets;

	Int_t is_2lss_TTH_SR, is_3l_TTH_SR, is_3l_TTZ_CR;

        Int_t           multilepton_Lepton1_Id,             multilepton_Lepton2_Id,                 multilepton_Lepton3_Id,             multilepton_Lepton4_Id;
        TLorentzVector  multilepton_Lepton1_P4,             multilepton_Lepton2_P4,                 multilepton_Lepton3_P4,             multilepton_Lepton4_P4;

        Float_t         multilepton_Lepton1_DeltaR_Matched,  multilepton_Lepton2_DeltaR_Matched,    multilepton_Lepton3_DeltaR_Matched,  multilepton_Lepton4_DeltaR_Matched;
        Int_t           multilepton_Lepton1_Label_Matched,   multilepton_Lepton2_Label_Matched,     multilepton_Lepton3_Label_Matched,   multilepton_Lepton4_Label_Matched;
        Int_t           multilepton_Lepton1_Id_Matched,      multilepton_Lepton2_Id_Matched,        multilepton_Lepton3_Id_Matched,      multilepton_Lepton4_Id_Matched;
        TLorentzVector  multilepton_Lepton1_P4_Matched,      multilepton_Lepton2_P4_Matched,        multilepton_Lepton3_P4_Matched,      multilepton_Lepton4_P4_Matched;

        Int_t multilepton_Bjet1_Id, multilepton_Bjet2_Id;
        TLorentzVector multilepton_Bjet1_P4, multilepton_Bjet2_P4;
        Float_t multilepton_Bjet1_CSV, multilepton_Bjet2_CSV;
        Float_t multilepton_Bjet1_JEC_Up, multilepton_Bjet1_JEC_Down, multilepton_Bjet2_JEC_Up, multilepton_Bjet2_JEC_Down;
        Float_t multilepton_Bjet1_JER_Up, multilepton_Bjet1_JER_Down, multilepton_Bjet2_JER_Up, multilepton_Bjet2_JER_Down;

        Float_t         multilepton_Bjet1_DeltaR_Matched,   multilepton_Bjet2_DeltaR_Matched;
        Int_t           multilepton_Bjet1_Label_Matched,    multilepton_Bjet2_Label_Matched;
        Int_t           multilepton_Bjet1_Id_Matched,       multilepton_Bjet2_Id_Matched;
        TLorentzVector  multilepton_Bjet1_P4_Matched,       multilepton_Bjet2_P4_Matched;

        Int_t multilepton_JetHighestPt1_Id, multilepton_JetHighestPt2_Id, multilepton_JetClosestMw1_Id, multilepton_JetClosestMw2_Id, multilepton_JetLowestMjj1_Id, multilepton_JetLowestMjj2_Id;
        TLorentzVector multilepton_JetHighestPt1_P4, multilepton_JetHighestPt2_P4, multilepton_JetClosestMw1_P4, multilepton_JetClosestMw2_P4, multilepton_JetLowestMjj1_P4, multilepton_JetLowestMjj2_P4;
        Float_t multilepton_JetHighestPt1_CSV, multilepton_JetHighestPt2_CSV, multilepton_JetClosestMw1_CSV, multilepton_JetClosestMw2_CSV, multilepton_JetLowestMjj1_CSV, multilepton_JetLowestMjj2_CSV;
        Float_t multilepton_JetHighestPt1_JEC_Up, multilepton_JetHighestPt2_JEC_Up, multilepton_JetClosestMw1_JEC_Up, multilepton_JetClosestMw2_JEC_Up, multilepton_JetLowestMjj1_JEC_Up, multilepton_JetLowestMjj2_JEC_Up;
        Float_t multilepton_JetHighestPt1_JEC_Down, multilepton_JetHighestPt2_JEC_Down, multilepton_JetClosestMw1_JEC_Down, multilepton_JetClosestMw2_JEC_Down, multilepton_JetLowestMjj1_JEC_Down, multilepton_JetLowestMjj2_JEC_Down;
        Float_t multilepton_JetHighestPt1_JER_Up, multilepton_JetHighestPt2_JER_Up, multilepton_JetClosestMw1_JER_Up, multilepton_JetClosestMw2_JER_Up, multilepton_JetLowestMjj1_JER_Up, multilepton_JetLowestMjj2_JER_Up;
        Float_t multilepton_JetHighestPt1_JER_Down, multilepton_JetHighestPt2_JER_Down, multilepton_JetClosestMw1_JER_Down, multilepton_JetClosestMw2_JER_Down, multilepton_JetLowestMjj1_JER_Down, multilepton_JetLowestMjj2_JER_Down;

        Float_t         multilepton_Jet1_DeltaR_Matched,   multilepton_Jet2_DeltaR_Matched;
        Int_t           multilepton_Jet1_Label_Matched,    multilepton_Jet2_Label_Matched;
        Int_t           multilepton_Jet1_Id_Matched,       multilepton_Jet2_Id_Matched;
        TLorentzVector  multilepton_Jet1_P4_Matched,       multilepton_Jet2_P4_Matched;

	TLorentzVector multilepton_W1_P4_Matched, multilepton_W2_P4_Matched;


        Int_t multilepton_JetHighestPt1_2ndPair_Id, multilepton_JetHighestPt2_2ndPair_Id, multilepton_JetClosestMw1_2ndPair_Id, multilepton_JetClosestMw2_2ndPair_Id, multilepton_JetLowestMjj1_2ndPair_Id, multilepton_JetLowestMjj2_2ndPair_Id;
        TLorentzVector  multilepton_JetHighestPt1_2ndPair_P4, multilepton_JetHighestPt2_2ndPair_P4, multilepton_JetClosestMw1_2ndPair_P4, multilepton_JetClosestMw2_2ndPair_P4, multilepton_JetLowestMjj1_2ndPair_P4, multilepton_JetLowestMjj2_2ndPair_P4;
        Float_t multilepton_JetHighestPt1_2ndPair_CSV, multilepton_JetHighestPt2_2ndPair_CSV, multilepton_JetClosestMw1_2ndPair_CSV, multilepton_JetClosestMw2_2ndPair_CSV, multilepton_JetLowestMjj1_2ndPair_CSV, multilepton_JetLowestMjj2_2ndPair_CSV;
        Float_t multilepton_JetHighestPt1_2ndPair_JEC_Up, multilepton_JetHighestPt2_2ndPair_JEC_Up, multilepton_JetClosestMw1_2ndPair_JEC_Up, multilepton_JetClosestMw2_2ndPair_JEC_Up, multilepton_JetLowestMjj1_2ndPair_JEC_Up, multilepton_JetLowestMjj2_2ndPair_JEC_Up;
        Float_t multilepton_JetHighestPt1_2ndPair_JEC_Down, multilepton_JetHighestPt2_2ndPair_JEC_Down, multilepton_JetClosestMw1_2ndPair_JEC_Down, multilepton_JetClosestMw2_2ndPair_JEC_Down, multilepton_JetLowestMjj1_2ndPair_JEC_Down, multilepton_JetLowestMjj2_2ndPair_JEC_Down;
        Float_t multilepton_JetHighestPt1_2ndPair_JER_Up, multilepton_JetHighestPt2_2ndPair_JER_Up, multilepton_JetClosestMw1_2ndPair_JER_Up, multilepton_JetClosestMw2_2ndPair_JER_Up, multilepton_JetLowestMjj1_2ndPair_JER_Up, multilepton_JetLowestMjj2_2ndPair_JER_Up;
        Float_t multilepton_JetHighestPt1_2ndPair_JER_Down, multilepton_JetHighestPt2_2ndPair_JER_Down, multilepton_JetClosestMw1_2ndPair_JER_Down, multilepton_JetClosestMw2_2ndPair_JER_Down, multilepton_JetLowestMjj1_2ndPair_JER_Down, multilepton_JetLowestMjj2_2ndPair_JER_Down;

        Int_t           multilepton_h0_Label,   multilepton_t1_Label,   multilepton_t2_Label;
        Int_t           multilepton_h0_Id,      multilepton_t1_Id,      multilepton_t2_Id;
        TLorentzVector  multilepton_h0_P4,      multilepton_t1_P4,      multilepton_t2_P4;

        TLorentzVector multilepton_mET, multilepton_Ptot;
        Float_t multilepton_mHT;
        Double_t multilepton_mETcov00;
        Double_t multilepton_mETcov01;
        Double_t multilepton_mETcov10;
        Double_t multilepton_mETcov11;

       TLorentzVector multilepton_mET_Matched;

        TH1F* TF_B[3][6];
        TH1F* TF_nonB[3][6];
	TH1F* TFratio_B[3][6];
	TH1F* TFratio_nonB[3][6];
	TH1F* TF_mET_Px[6];
        TH1F* TF_mET_Py[6];
        TH1F* TF_mET_Pt[6];
        TH1F* TF_mET_Phi[6];
	TH1F* Pdf_Ptot;
	TH1F* Pdf_PtotEta;
        TH1F* Pdf_PtotM;
	TH3F* Pdf_PtotP4;
};

MultiLeptonTree::MultiLeptonTree(){
}

MultiLeptonTree::~MultiLeptonTree(){
}

void MultiLeptonTree::initializeOutputTree(){

    tOutput = new TTree("Tree", "Tree");
    tOutput->Branch("weight",&weight,"weight/F");
    tOutput->Branch("mc_weight",&mc_weight,"mc_weight/F"); 
    tOutput->Branch("PV_weight",&weight_PV,"PV_weight/F");

    tOutput->Branch("mc_ttZhypAllowed",&mc_ttZhypAllowed,"mc_ttZhypAllowed/I");
    tOutput->Branch("catJets",&catJets,"catJets/I");

    tOutput->Branch("is_2lss_TTH_SR",&is_2lss_TTH_SR,"is_2lss_TTH_SR/B");
    tOutput->Branch("is_3l_TTH_SR",&is_3l_TTH_SR,"is_3l_TTH_SR/B");
    tOutput->Branch("is_3l_TTZ_CR",&is_3l_TTZ_CR,"is_3l_TTZ_CR/B");

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

    tOutput->Branch("multilepton_W1_P4_Matched", "TLorentzVector", &multilepton_W1_P4_Matched);
    tOutput->Branch("multilepton_W2_P4_Matched", "TLorentzVector", &multilepton_W2_P4_Matched);


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

    tOutput->Branch("multilepton_h0_Id",                    &multilepton_h0_Id,                 "multilepton_h0_Id/I");
    tOutput->Branch("multilepton_h0_P4",                    "TLorentzVector",                   &multilepton_h0_P4);
    tOutput->Branch("multilepton_h0_Label",                    &multilepton_h0_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t1_Id",                    &multilepton_t1_Id,                 "multilepton_t1_Id/I");
    tOutput->Branch("multilepton_t1_P4",                    "TLorentzVector",                   &multilepton_t1_P4);
    tOutput->Branch("multilepton_t1_Label",                    &multilepton_t1_Label,                 "multilepton_h0_Label/I");
    tOutput->Branch("multilepton_t2_Id",                    &multilepton_t2_Id,                 "multilepton_t2_Id/I");
    tOutput->Branch("multilepton_t2_P4",                    "TLorentzVector",                   &multilepton_t2_P4);
    tOutput->Branch("multilepton_t2_Label",                    &multilepton_t2_Label,                 "multilepton_h0_Label/I");

    tOutput->Branch("multilepton_mET","TLorentzVector",&multilepton_mET);
    tOutput->Branch("multilepton_mETcov00",&multilepton_mETcov00,"multilepton_mETcov00/D");
    tOutput->Branch("multilepton_mETcov01",&multilepton_mETcov01,"multilepton_mETcov01/D");
    tOutput->Branch("multilepton_mETcov10",&multilepton_mETcov10,"multilepton_mETcov10/D");
    tOutput->Branch("multilepton_mETcov11",&multilepton_mETcov11,"multilepton_mETcov11/D");

    tOutput->Branch("multilepton_mET_Matched", "TLorentzVector", &multilepton_mET_Matched);

    tOutput->Branch("multilepton_mHT",&multilepton_mHT,"multilepton_mHT/F");
    tOutput->Branch("multilepton_Ptot","TLorentzVector",&multilepton_Ptot);

    //TFratio_B = new TH1F("TFratio_B","TFratio_B", 100, 0, 5);
    //TFratio_nonB = new TH1F("TFratio_nonB", "TFratio_nonB", 100, 0, 5);

    Pdf_Ptot = new TH1F("Pdf_Ptot","Pdf_Ptot",1000,0,2000);
    Pdf_PtotEta = new TH1F("Pdf_PtotEta","Pdf_PtotEta",200, -7, 7);
    Pdf_PtotM = new TH1F("Pdf_PtotM","Pdf_PtotM",1000,0,2000);
    Pdf_PtotP4 = new TH3F("Pdf_PtotP4","Pdf_PtotP4",1000,0,2000, 200, -7, 7, 1000, 0, 2000);

    string histname;
    for (int ieta=0; ieta<3; ieta++)  {
      for (int ienergy=0; ienergy<6; ienergy++) {
        histname = "TF_B_" + EtaRangeLabel[ieta] + "_" + EnergyRangeLabel[ienergy];
        TF_B[ieta][ienergy] = new TH1F(histname.c_str(), histname.c_str(), 100, 0, 4);
        histname = "TF_nonB_" + EtaRangeLabel[ieta] + "_" + EnergyRangeLabel[ienergy];
        TF_nonB[ieta][ienergy] = new TH1F(histname.c_str(), histname.c_str(), 100, 0, 4);

	histname = "TFratio_B_" + EtaRangeLabel[ieta] + "_" + EnergyRangeLabel[ienergy];
	TFratio_B[ieta][ienergy] = new TH1F(histname.c_str(), histname.c_str(), 100, 0, 4);
	histname = "TFratio_nonB_" + EtaRangeLabel[ieta] + "_" + EnergyRangeLabel[ienergy];
	TFratio_nonB[ieta][ienergy] = new TH1F(histname.c_str(), histname.c_str(), 100, 0, 4);
      }
    }

    for (int imet=0; imet<6; imet++){
      histname = "TF_mET_Px_" + MetRangeLabel[imet];
      TF_mET_Px[imet] = new TH1F(histname.c_str(), histname.c_str(), 100, -150, 150);
      histname = "TF_mET_Py_" + MetRangeLabel[imet];
      TF_mET_Py[imet] = new TH1F(histname.c_str(), histname.c_str(), 100, -150, 150);
      histname = "TF_mET_Pt_" + MetRangeLabel[imet];
      TF_mET_Pt[imet] = new TH1F(histname.c_str(), histname.c_str(), 100, -150, 150);
      histname = "TF_mET_Phi_" + MetRangeLabel[imet];
      TF_mET_Phi[imet] = new TH1F(histname.c_str(), histname.c_str(), 100, -3.2, 3.2);
    }

}

void MultiLeptonTree::FillJetInfoOutputTree(int* tree_Id, int id, TLorentzVector* tree_P4, Jet* jet){

  *tree_Id = id;
  (*tree_P4).SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, 0);

  return;
}


#endif



