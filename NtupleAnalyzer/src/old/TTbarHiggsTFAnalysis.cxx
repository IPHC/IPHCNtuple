#include "../include/TTbarHiggsTFAnalysis.h"
#include "TSystem.h"

TTbarHiggsTFAnalysis::TTbarHiggsTFAnalysis() 
{

}


TTbarHiggsTFAnalysis::TTbarHiggsTFAnalysis(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName, bool isdata, float xsec, float lumi, int nowe, int nmax)
{    

    _isdata = isdata;
    _xsec = xsec;
    _lumi = lumi;
    _nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    _sampleName = sampleName;

    //
    //_file_PVreweighting = TFile::Open("/home-pbs/lebihan/someone/ttH_070116/ttH/NtupleAnalyzer/test/PUweight.root");
    //_h_PV = (TH1F*)_file_PVreweighting->Get("PU_reweighting");

    //
    tree = new TChain(treeName.Data());

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);
        tree->Add(fnameStr.c_str());
    }   
    infile.close();
    Init(tree);

    theHistoManager = new HistoManager();

    TString outputfileNameRoot = _outputFileName+".root";
    outputfile = new TFile(outputfileNameRoot.Data(), "recreate");  

}


void TTbarHiggsTFAnalysis::createHistograms()
{    
    outputfile->cd();
    initializeOutputTree();


//**********************************
// Data
//**********************************
 theHistoManager->addHisto("hData","nJets","All","TEST",61,-0.5,60.5);
 theHistoManager->addHisto("hData","nJets","All","",61,-0.5,60.5);
 theHistoManager->addHisto("hData","Run","All","",280,190000,204000);
 theHistoManager->addHisto("hData","JetPt","All","", 28, 20., 300.);
 theHistoManager->addHisto("hData","JetEta","All","",40, -4., 4.);

 theHistoManager->addHisto("hData","nJets","","",61,-0.5,60.5);
 theHistoManager->addHisto("hData","JetPt","","", 28, 20., 300.);
 theHistoManager->addHisto("hData","JetEta","","",40, -4., 4.);

 theHistoManager->addHisto("hData","metpx","","",100,-150.,150.);
 theHistoManager->addHisto("hData","metpy","","",100,-150.,150.);
 theHistoManager->addHisto("hData","metpt","","",100,0.,300.);
 theHistoManager->addHisto("hData","metphi","","",64,-3.2,3.2);
 theHistoManager->addHisto("hData","metsum","","",100,0.,3000.);
   
//**********************************
// All flavours in Monte Carlo
//**********************************
 theHistoManager->addHisto("hGenlep","nlep","","",9,0.5,9.5);
 theHistoManager->addHisto("hGenlep","emu q","","",9,-4.5,4.5);
 theHistoManager->addHisto("hGenlep","emu ptmax","","",30,0.,300.);
 theHistoManager->addHisto("hGenlep","nemu","","",9,0.5,9.5);
 theHistoManager->addHisto("hGenlep","dR","ele","",50,0.,0.5);
 theHistoManager->addHisto("hGenlep","dR","muo","",50,0.,0.5);
 theHistoManager->addHisto("hGenlep","dR","tau","",50,0.,0.5);

 theHistoManager->addHisto("hAllFlav","pthat","All","",50, 0., 500.);
 theHistoManager->addHisto("hAllFlav","Flavour","All","",22,-0.5,21.5);

 theHistoManager->addHisto("hAllFlav","Flavour","","",22,-0.5,21.5);

 theHistoManager->addHisto("hAllFlav","jes","etaLT18","20genpt30",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","etaLT18","30genpt50",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","etaLT18","50genpt80",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","etaLT18","80genpt120",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","etaLT18","120genpt200",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","etaLT18","genptGT200",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","20genpt30",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","30genpt50",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","50genpt80",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","80genpt120",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","120genpt200",75,0.,3.);
 theHistoManager->addHisto("hAllFlav","jes","18eta24","genptGT200",75,0.,3.);

 theHistoManager->addHisto("hAllFlav","TF","metpx","",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin1",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin1",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin1",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin1",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin2",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin2",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin2",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin2",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin3",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin3",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin3",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin3",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin4",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin4",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin4",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin4",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin5",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin5",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin5",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin5",320,-3.2,3.2);
 theHistoManager->addHisto("hAllFlav","TF","metpx","bin6",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpy","bin6",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metpt","bin6",150,-150.,150.);
 theHistoManager->addHisto("hAllFlav","TF","metphi","bin6",320,-3.2,3.2);

//  theHistoManager->addHisto("hAllFlav","TF","metpx","set1",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set1",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set1",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set1",320,-3.2,3.2);
//  theHistoManager->addHisto("hAllFlav","TF","metpx","set2",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set2",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set2",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set2",320,-3.2,3.2);
//  theHistoManager->addHisto("hAllFlav","TF","metpx","set3",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set3",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set3",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set3",320,-3.2,3.2);
//  theHistoManager->addHisto("hAllFlav","TF","metpx","set4",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set4",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set4",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set4",320,-3.2,3.2);
//  theHistoManager->addHisto("hAllFlav","TF","metpx","set5",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set5",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set5",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set5",320,-3.2,3.2);
//  theHistoManager->addHisto("hAllFlav","TF","metpx","set6",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpy","set6",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metpt","set6",150,-150.,150.);
//  theHistoManager->addHisto("hAllFlav","TF","metphi","set6",320,-3.2,3.2);

 theHistoManager->addHisto("hUFlav","TF","","",100,-200.,200.);
 theHistoManager->addHisto("hUFlav","TF","bin1","",100,-40.,100.);
 theHistoManager->addHisto("hUFlav","TF","bin2","",100,-60.,120.);
 theHistoManager->addHisto("hUFlav","TF","bin3","",100,-80.,140.);
 theHistoManager->addHisto("hUFlav","TF","bin4","",100,-150.,150.);
 theHistoManager->addHisto("hUFlav","TF","bin5","",100,-200.,200.);
 theHistoManager->addHisto("hUFlav","TF","bin6","",100,-600.,600.);
 theHistoManager->addHisto("hUFlav","TX","",""    ,200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin1","",200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin2","",200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin3","",200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin4","",200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin5","",200,0.,4.);
 theHistoManager->addHisto("hUFlav","TX","bin6","",200,0.,4.);

 theHistoManager->addHisto("hBFlav","TF","","",100,-200.,200.);
 theHistoManager->addHisto("hBFlav","TF","bin1","",100,-40.,100.);
 theHistoManager->addHisto("hBFlav","TF","bin2","",100,-60.,120.);
 theHistoManager->addHisto("hBFlav","TF","bin3","",100,-80.,140.);
 theHistoManager->addHisto("hBFlav","TF","bin4","",100,-150.,150.);
 theHistoManager->addHisto("hBFlav","TF","bin5","",100,-200.,200.);
 theHistoManager->addHisto("hBFlav","TF","bin6","",100,-600.,600.);
 theHistoManager->addHisto("hBFlav","TX","",""    ,200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin1","",200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin2","",200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin3","",200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin4","",200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin5","",200,0.,4.);
 theHistoManager->addHisto("hBFlav","TX","bin6","",200,0.,4.);
}


void TTbarHiggsTFAnalysis::writeHistograms()
{  
    outputfile->cd();

    std::vector<TH1F*> the1DHisto =  theHistoManager->getHisto1D_list();
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();

    for(unsigned int i=0; i<the1DHisto.size(); i++)  the1DHisto[i]->Write();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();

    tOutput->Write();
    outputfile->Close();
}


void TTbarHiggsTFAnalysis::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

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


void TTbarHiggsTFAnalysis::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries << std::endl;
    std::cout << "Number of processed events = " << nentries_max << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max;jentry++) 
    {

        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
        //std::cout << "number of processed events " << jentry << std::endl;

        //if(jentry > 100000) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
	//
	int pvn = vEvent->at(0).pv_n();
        theHistoManager->fillHisto("NumberOfPrimaryVertex", "noSel", "", "",  pvn, 1);
	
	if ( !_isdata )
	{
          weight = _lumi*_xsec/_nowe;
          mc_weight = vEvent->at(0).mc_weight();
          //weight_PV = _h_PV->GetBinContent(pvn);
          weight = weight * mc_weight; //*weight_PV;
	}
	else 
	{
	    weight    = 1.;
        mc_weight = 1.;
        weight_PV = 1.; 
	}

        //---------------------------
        // initialisation
        //---------------------------
        vLeptons.clear();
        vSelectedMuons.clear();
        vSelectedElectrons.clear();
        vSelectedLeptons.clear();
        vFakeMuons.clear();
        vFakeElectrons.clear();
        vFakeLeptons.clear();	
        vSelectedJets.clear();

  theHistoManager->fillHisto("hAllFlav","pthat","All","", vEvent->at(0).mc_ptHat(), weight);
  theHistoManager->fillHisto("hData","Run","All","", vEvent->at(0).run(), weight);

//*********************************
// primary leptons

  int nlep = 0, nemu = 0, qemu = 0;
  float lep_pt[20], lep_eta[20], lep_phi[20], lep_id[20];
  float ptmax = -1.;

// std::cout << std::endl;
  for (unsigned int i=0; i < vTruth->at(0).mc_truth_n() ; i++)
  {
    int ID = abs( vTruth->at(0).mc_truth_id().at(i) );
    int label = vTruth->at(0).mc_truth_label().at(i);
    
    if ( nlep < 20 && (ID == 11 || ID == 13 || ID == 15)  
        && (label == 120 || label == 124 || label == 1240 // from h->WW
         || label == 130 || label == 134 || label == 1340 // from h->WW
         || label == 140 || label == 141 || label == 144 || label == 145 || label == 1440 || label == 1450 // from h->ZZ 
         || label == 150 || label == 151 || label == 154 || label == 155 || label == 1540 || label == 1550 // from h->ZZ 
         || label ==  16 || label ==  17 || label == 160 || label == 170 // from tau
         || label ==  22 || label == 210 || label == 2220 // from top->W 
         || label ==  32 || label == 310 || label == 3220 // from top->W 
         || label ==  40 || label ==  43 || label == 430 // from W 
         || label ==  50 || label ==  51 || label ==  52 || label ==  53 || label == 520 || label == 530) // from Z 
	&& vTruth->at(0).mc_truth_pt().at(i) > 10. 
	&& abs(vTruth->at(0).mc_truth_eta().at(i)) < 2.4 ) {
      lep_pt[nlep]  = vTruth->at(0).mc_truth_pt().at(i);
      lep_eta[nlep] = vTruth->at(0).mc_truth_eta().at(i);
      lep_phi[nlep] = vTruth->at(0).mc_truth_phi().at(i);
      lep_id[nlep]  = vTruth->at(0).mc_truth_id().at(i);
//$$
// std::cout << nlep << " " << lep_pt[nlep] << " " << lep_eta[nlep]
// 	       << " " << lep_phi[nlep] << " " << lep_id[nlep] << " " << label << std::endl; 
//$$
      
      if ( abs(lep_id[nlep]) == 11 || abs(lep_id[nlep]) == 13 ) nemu++;
      if ( lep_id[nlep] == 11 || lep_id[nlep] == 13 ) qemu+=-1;
      if ( lep_id[nlep] ==-11 || lep_id[nlep] ==-13 ) qemu+= 1;
      if ( (abs(lep_id[nlep]) == 11 || abs(lep_id[nlep]) == 13) &&  
           lep_pt[nlep] > ptmax ) ptmax = lep_pt[nlep];

      nlep++;
    }
  }

  theHistoManager->fillHisto("hGenlep","nlep","","", nlep, weight);
  theHistoManager->fillHisto("hGenlep","emu q","","", qemu, weight);
  theHistoManager->fillHisto("hGenlep","emu ptmax","","", ptmax, weight);
  if ( nemu >= 2 && ptmax > 20. && (nemu >= 3 || qemu != 0) )
  theHistoManager->fillHisto("hGenlep","nemu","","", nemu, weight);

//$$
//$$  if ( nemu < 3 || ptmax < 20. ) continue;
//$$  if ( nemu != 2 || qemu == 0 || ptmax < 20. ) continue;
//$$


//*********************************
// loop on jets

  int njetall = 0, njet = 0;

  for (unsigned int ijet=0; ijet < vJet->size() ; ijet++)
  {
   njetall++;
   float ptjet = vJet->at(ijet).pt();
   float Ejet = vJet->at(ijet).E();
   float etajet = vJet->at(ijet).eta();

   float genpt = vJet->at(ijet).jet_genJet_pt();
   float jes = 0.;
   if ( genpt > 0. ) jes = ptjet / genpt;
//$$
   int flavour = abs(vJet->at(ijet).jet_partonFlavour());
   if ( flavour >= 1 && flavour <= 3 ) flavour = 1; 
   if ( flavour >= 4 && flavour <= 5 ) flavour = 0; 
   if ( abs(vJet->at(ijet).jet_hadronFlavour()) >= 4 && 
        abs(vJet->at(ijet).jet_hadronFlavour()) <= 5 ) flavour = abs(vJet->at(ijet).jet_hadronFlavour()); 
//$$

   theHistoManager->fillHisto("hAllFlav","Flavour","All","", flavour, weight);

   theHistoManager->fillHisto("hData","JetPt","All","", ptjet, weight);
   theHistoManager->fillHisto("hData","JetEta","All","", etajet, weight);

//$$
   if ( !(ptjet >= 25. && ptjet < 1000.) ) continue;
   if ( !(fabs(etajet) < 2.4) ) continue;
//$$   if ( !(fabs(etajet) < 0.8) ) continue;
//$$   if ( !(fabs(etajet) >= 0.8 && fabs(etajet) < 1.6) ) continue;
//$$   if ( !(fabs(etajet) >= 1.6 && fabs(etajet) < 2.4) ) continue;
   if ( flavour == 0 ) continue;
   if ( genpt < 10. ) continue;
//$$ 

   if ( ptjet >= 300.) ptjet = 299.;


//*********************************
// jets due to leptons

   float dRlj;
   float dRele = 1000., dRmuo = 1000., dRtau = 1000.;
   for (int k = 0; k < nlep; k++) {
     int lepid = abs(lep_id[k]);
     dRlj = GetDeltaR( etajet, vJet->at(ijet).phi(), lep_eta[k], lep_phi[k] );
     if ( lepid == 11 && dRlj < dRele ) dRele = dRlj;
     if ( lepid == 13 && dRlj < dRmuo ) dRmuo = dRlj;
     if ( lepid == 15 && dRlj < dRtau ) dRtau = dRlj;
   }
   if ( dRele < 0.5 ) theHistoManager->fillHisto("hGenlep","dR","ele","", dRele, weight);
   if ( dRmuo < 0.5 && dRele > 0.4 ) theHistoManager->fillHisto("hGenlep","dR","muo","", dRmuo, weight);
   if ( dRtau < 0.5 && dRmuo > 0.4 && dRele > 0.4 ) theHistoManager->fillHisto("hGenlep","dR","tau","", dRtau, weight);
//$$ 
   if ( dRele < 0.4 || dRmuo < 0.4 || dRtau < 0.4 ) flavour = -1;
//$$ 

//$$ 
 if ( flavour <= 0 ) continue;
//$$ 
   njet++;

//*********************************

   theHistoManager->fillHisto("hAllFlav","Flavour","","", flavour, weight);

   theHistoManager->fillHisto("hData","JetPt","","", ptjet, weight);
   theHistoManager->fillHisto("hData","JetEta","","", etajet, weight);

   if ( abs(etajet) < 1.8 ) {
     if ( genpt >= 20. && genpt < 30. )   theHistoManager->fillHisto("hAllFlav","jes","etaLT18","20genpt30", jes, weight);
     if ( genpt >= 30. && genpt < 50. )   theHistoManager->fillHisto("hAllFlav","jes","etaLT18","30genpt50", jes, weight);
     if ( genpt >= 50. && genpt < 80. )   theHistoManager->fillHisto("hAllFlav","jes","etaLT18","50genpt80", jes, weight);
     if ( genpt >= 80. && genpt < 120. )  theHistoManager->fillHisto("hAllFlav","jes","etaLT18","80genpt120", jes, weight);
     if ( genpt >= 120. && genpt < 200. ) theHistoManager->fillHisto("hAllFlav","jes","etaLT18","120genpt200", jes, weight);
     if ( genpt >= 200. )		  theHistoManager->fillHisto("hAllFlav","jes","etaLT18","genptGT200", jes, weight);
   }
   else if ( abs(etajet) < 2.4 ) {
     if ( genpt >= 20. && genpt < 30. )   theHistoManager->fillHisto("hAllFlav","jes","18eta24","20genpt30", jes, weight);
     if ( genpt >= 30. && genpt < 50. )   theHistoManager->fillHisto("hAllFlav","jes","18eta24","30genpt50", jes, weight);
     if ( genpt >= 50. && genpt < 80. )   theHistoManager->fillHisto("hAllFlav","jes","18eta24","50genpt80", jes, weight);
     if ( genpt >= 80. && genpt < 120. )  theHistoManager->fillHisto("hAllFlav","jes","18eta24","80genpt120", jes, weight);
     if ( genpt >= 120. && genpt < 200. ) theHistoManager->fillHisto("hAllFlav","jes","18eta24","120genpt200", jes, weight);
     if ( genpt >= 200. )		  theHistoManager->fillHisto("hAllFlav","jes","18eta24","genptGT200", jes, weight);
   }

//*********************************
// TF study
   
   bool hasParton = false;
   float TF_pt = -1000., TF_E = -1000., TX_E = -1000.; 

   if ( vJet->at(ijet).jet_genParton_pt() > 0. ) {
     int quaid = TMath::Abs( vJet->at(ijet).jet_genParton_id() );
     if ( (flavour != 5 && quaid != 5) || (flavour == 5 && quaid == 5) ) {
       if ( quaid > 0 && quaid <= 5 ) hasParton = true;
     }
   }
   if ( hasParton ) {
     TF_E  = vJet->at(ijet).jet_genParton_E() - vJet->at(ijet).E();
     TX_E  = vJet->at(ijet).E() / vJet->at(ijet).jet_genParton_E();
   }

   int binpt = -1;
   if      ( Ejet < 50. ) binpt = 1;
   else if ( Ejet < 80. ) binpt = 2;
   else if ( Ejet <120. ) binpt = 3;
   else if ( Ejet <200. ) binpt = 4;
   else if ( Ejet <300. ) binpt = 5;
   else                   binpt = 6;

   if ( TX_E > 0 ) {
     if ( flavour == 5 ) {
       theHistoManager->fillHisto("hBFlav","TF","","", TF_E, weight);
       theHistoManager->fillHisto("hBFlav","TX","","", TX_E, weight);
       if ( binpt == 1 ) theHistoManager->fillHisto("hBFlav","TF","bin1","", TF_E, weight);
       if ( binpt == 2 ) theHistoManager->fillHisto("hBFlav","TF","bin2","", TF_E, weight);
       if ( binpt == 3 ) theHistoManager->fillHisto("hBFlav","TF","bin3","", TF_E, weight);
       if ( binpt == 4 ) theHistoManager->fillHisto("hBFlav","TF","bin4","", TF_E, weight);
       if ( binpt == 5 ) theHistoManager->fillHisto("hBFlav","TF","bin5","", TF_E, weight);
       if ( binpt == 6 ) theHistoManager->fillHisto("hBFlav","TF","bin6","", TF_E, weight);
       if ( binpt == 1 ) theHistoManager->fillHisto("hBFlav","TX","bin1","", TX_E, weight);
       if ( binpt == 2 ) theHistoManager->fillHisto("hBFlav","TX","bin2","", TX_E, weight);
       if ( binpt == 3 ) theHistoManager->fillHisto("hBFlav","TX","bin3","", TX_E, weight);
       if ( binpt == 4 ) theHistoManager->fillHisto("hBFlav","TX","bin4","", TX_E, weight);
       if ( binpt == 5 ) theHistoManager->fillHisto("hBFlav","TX","bin5","", TX_E, weight);
       if ( binpt == 6 ) theHistoManager->fillHisto("hBFlav","TX","bin6","", TX_E, weight);
     }
     else {
       theHistoManager->fillHisto("hUFlav","TF","","", TF_E, weight);
       theHistoManager->fillHisto("hUFlav","TX","","", TX_E, weight);
       if ( binpt == 1 ) theHistoManager->fillHisto("hUFlav","TF","bin1","", TF_E, weight);
       if ( binpt == 2 ) theHistoManager->fillHisto("hUFlav","TF","bin2","", TF_E, weight);
       if ( binpt == 3 ) theHistoManager->fillHisto("hUFlav","TF","bin3","", TF_E, weight);
       if ( binpt == 4 ) theHistoManager->fillHisto("hUFlav","TF","bin4","", TF_E, weight);
       if ( binpt == 5 ) theHistoManager->fillHisto("hUFlav","TF","bin5","", TF_E, weight);
       if ( binpt == 6 ) theHistoManager->fillHisto("hUFlav","TF","bin6","", TF_E, weight);
       if ( binpt == 1 ) theHistoManager->fillHisto("hUFlav","TX","bin1","", TX_E, weight);
       if ( binpt == 2 ) theHistoManager->fillHisto("hUFlav","TX","bin2","", TX_E, weight);
       if ( binpt == 3 ) theHistoManager->fillHisto("hUFlav","TX","bin3","", TX_E, weight);
       if ( binpt == 4 ) theHistoManager->fillHisto("hUFlav","TX","bin4","", TX_E, weight);
       if ( binpt == 5 ) theHistoManager->fillHisto("hUFlav","TX","bin5","", TX_E, weight);
       if ( binpt == 6 ) theHistoManager->fillHisto("hUFlav","TX","bin6","", TX_E, weight);
     }
   }

  } // end loop on jets

  theHistoManager->fillHisto("hData","nJets","All","TEST", njetall, weight);
  theHistoManager->fillHisto("hData","nJets","All","", njetall, weight);
  theHistoManager->fillHisto("hData","nJets","","", njet, weight);

//*********************************
// TF MET study
   
  float MET = vEvent->at(0).metpt();
  float METphi = vEvent->at(0).metphi();
  float METx = MET * TMath::Cos(METphi);
  float METy = MET * TMath::Sin(METphi);
  float METsum = vEvent->at(0).metsumet();

  float METGEN = vTruth->at(0).metGen_pt();
  float METGENphi = vTruth->at(0).metGen_phi();
  float METGENx = vTruth->at(0).metGen_px();
  float METGENy = vTruth->at(0).metGen_py();

  float dphi =  GetSignedDPhi( METGENphi , METphi );

//   int binpt = -1;
//   if ( MET < 20. )	 binpt = 1;
//   else if ( MET < 40. )	 binpt = 2;
//   else if ( MET < 60. )	 binpt = 3;
//   else if ( MET < 100. ) binpt = 4;
//   else if ( MET < 150. ) binpt = 5;
//   else  		 binpt = 6;
// 
//   theHistoManager->fillHisto("hAllFlav","TF","metpx","",METGENx - METx, weight);
//   theHistoManager->fillHisto("hAllFlav","TF","metpy","",METGENy - METy, weight);
//   theHistoManager->fillHisto("hAllFlav","TF","metpt","",METGEN - MET, weight);
//   theHistoManager->fillHisto("hAllFlav","TF","metphi","",dphi, weight);
//   
//   if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin1",METGENx - METx, weight);
//   if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin1",METGENy - METy, weight);
//   if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin1",METGEN - MET, weight);
//   if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin1",dphi, weight);
//   if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin2",METGENx - METx, weight);
//   if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin2",METGENy - METy, weight);
//   if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin2",METGEN - MET, weight);
//   if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin2",dphi, weight);
//   if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin3",METGENx - METx, weight);
//   if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin3",METGENy - METy, weight);
//   if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin3",METGEN - MET, weight);
//   if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin3",dphi, weight);
//   if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin4",METGENx - METx, weight);
//   if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin4",METGENy - METy, weight);
//   if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin4",METGEN - MET, weight);
//   if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin4",dphi, weight);
//   if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin5",METGENx - METx, weight);
//   if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin5",METGENy - METy, weight);
//   if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin5",METGEN - MET, weight);
//   if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin5",dphi, weight);
//   if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin6",METGENx - METx, weight);
//   if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin6",METGENy - METy, weight);
//   if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin6",METGEN - MET, weight);
//   if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin6",dphi, weight);
// 
//   if ( MET >= 300. ) MET = 299.;
//   if ( METx >= 150. ) METx = 149.;
//   if ( METy >= 150. ) METy = 149.;
//   if ( METx <=-150. ) METx =-149.;
//   if ( METy <=-150. ) METy =-149.;
//   if ( MET >= 3000. ) METsum = 2999.;
//   
//   theHistoManager->fillHisto("hData","metpx","","", METx, weight);
//   theHistoManager->fillHisto("hData","metpy","","", METy, weight);
//   theHistoManager->fillHisto("hData","metpt","","", MET, weight);
//   theHistoManager->fillHisto("hData","metphi","","", METphi, weight);
//   theHistoManager->fillHisto("hData","metsum","","", METsum, weight);
//    
//   int sumet = -1;
//   if ( METsum < 1000. )	     sumet = 1;
//   else if ( METsum < 1200. ) sumet = 2;
//   else if ( METsum < 1600. ) sumet = 3;
//   else if ( METsum < 2000. ) sumet = 4;
//   else if ( METsum < 2500. ) sumet = 5;
//   else  		     sumet = 6;
// 
//   if ( sumet == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set1",METGENx - METx, weight);
//   if ( sumet == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set1",METGENy - METy, weight);
//   if ( sumet == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set1",METGEN - MET, weight);
//   if ( sumet == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set1",dphi, weight);
//   if ( sumet == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set2",METGENx - METx, weight);
//   if ( sumet == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set2",METGENy - METy, weight);
//   if ( sumet == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set2",METGEN - MET, weight);
//   if ( sumet == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set2",dphi, weight);
//   if ( sumet == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set3",METGENx - METx, weight);
//   if ( sumet == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set3",METGENy - METy, weight);
//   if ( sumet == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set3",METGEN - MET, weight);
//   if ( sumet == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set3",dphi, weight);
//   if ( sumet == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set4",METGENx - METx, weight);
//   if ( sumet == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set4",METGENy - METy, weight);
//   if ( sumet == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set4",METGEN - MET, weight);
//   if ( sumet == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set4",dphi, weight);
//   if ( sumet == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set5",METGENx - METx, weight);
//   if ( sumet == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set5",METGENy - METy, weight);
//   if ( sumet == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set5",METGEN - MET, weight);
//   if ( sumet == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set5",dphi, weight);
//   if ( sumet == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","set6",METGENx - METx, weight);
//   if ( sumet == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","set6",METGENy - METy, weight);
//   if ( sumet == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","set6",METGEN - MET, weight);
//   if ( sumet == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","set6",dphi, weight);
// 

  int binpt = -1;
  if      ( MET < 100. && METsum < 1200 ) binpt = 1;
  else if ( MET < 100. && METsum < 1600 ) binpt = 2;
  else if ( MET < 100. )	          binpt = 3;
  else if (               METsum < 1200 ) binpt = 4;
  else if (               METsum < 1600 ) binpt = 5;
  else  		                  binpt = 6;

  theHistoManager->fillHisto("hAllFlav","TF","metpx","",METGENx - METx, weight);
  theHistoManager->fillHisto("hAllFlav","TF","metpy","",METGENy - METy, weight);
  theHistoManager->fillHisto("hAllFlav","TF","metpt","",METGEN - MET, weight);
  theHistoManager->fillHisto("hAllFlav","TF","metphi","",dphi, weight);
  
  if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin1",METGENx - METx, weight);
  if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin1",METGENy - METy, weight);
  if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin1",METGEN - MET, weight);
  if ( binpt == 1 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin1",dphi, weight);
  if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin2",METGENx - METx, weight);
  if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin2",METGENy - METy, weight);
  if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin2",METGEN - MET, weight);
  if ( binpt == 2 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin2",dphi, weight);
  if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin3",METGENx - METx, weight);
  if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin3",METGENy - METy, weight);
  if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin3",METGEN - MET, weight);
  if ( binpt == 3 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin3",dphi, weight);
  if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin4",METGENx - METx, weight);
  if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin4",METGENy - METy, weight);
  if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin4",METGEN - MET, weight);
  if ( binpt == 4 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin4",dphi, weight);
  if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin5",METGENx - METx, weight);
  if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin5",METGENy - METy, weight);
  if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin5",METGEN - MET, weight);
  if ( binpt == 5 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin5",dphi, weight);
  if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpx","bin6",METGENx - METx, weight);
  if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpy","bin6",METGENy - METy, weight);
  if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metpt","bin6",METGEN - MET, weight);
  if ( binpt == 6 ) theHistoManager->fillHisto("hAllFlav","TF","metphi","bin6",dphi, weight);

  if ( MET >= 300. ) MET = 299.;
  if ( METx >= 150. ) METx = 149.;
  if ( METy >= 150. ) METy = 149.;
  if ( METx <=-150. ) METx =-149.;
  if ( METy <=-150. ) METy =-149.;
  if ( MET >= 3000. ) METsum = 2999.;
  
  theHistoManager->fillHisto("hData","metpx","","", METx, weight);
  theHistoManager->fillHisto("hData","metpy","","", METy, weight);
  theHistoManager->fillHisto("hData","metpt","","", MET, weight);
  theHistoManager->fillHisto("hData","metphi","","", METphi, weight);
  theHistoManager->fillHisto("hData","metsum","","", METsum, weight);
   
 } // end loop on events

}

void TTbarHiggsTFAnalysis::initializeOutputTree()
{

    outputfile->cd();
    tOutput = new TTree("Tree", "Tree");

    return;
}


float TTbarHiggsTFAnalysis::GetSignedDPhi(float phi1,float phi2)
{
    float deltaPhi = phi1 - phi2;
    if (deltaPhi > 3.141593 ) deltaPhi -= 2.*3.141593;
    if (deltaPhi <-3.141593 ) deltaPhi -=-2.*3.141593;
    return deltaPhi;
}


float TTbarHiggsTFAnalysis::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
    float DeltaPhi = TMath::Abs(phi2 - phi1);
    if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
    return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
