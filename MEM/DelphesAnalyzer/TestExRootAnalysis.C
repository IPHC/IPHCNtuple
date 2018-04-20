#include <iostream>

#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#endif

#include "TestExRootAnalysis.h"
#include "MultiLeptonTree.h"

#define kCat_3l_2b_2j   0
#define kCat_3l_1b_2j   1
#define kCat_3l_2b_1j   2
#define kCat_3l_1b_1j   3
#define kCat_3l_2b_0j   4
#define kCat_4l_2b      5
#define kCat_4l_1b      6
#define kCat_2lss_2b_4j 7
#define kCat_2lss_1b_4j 8
#define kCat_2lss_2b_3j 9
#define kCat_2lss_1b_3j 10
#define kCat_2lss_2b_2j 11

//#include "ExRootAnalysis/ExRootTreeReader.h"
//#include "ExRootAnalysis/ExRootClasses.h"
//#include "/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/ExRootAnalysis/ExRootAnalysis/ExRootTreeReader.h"
//#include "/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"
//#include "/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
//#include "/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/Delphes/external/ExRootAnalysis/ExRootClasses.h"

using namespace std;

void TestExRootAnalysis()
{

  // Create chain of root trees
  TChain chain("Delphes");

  //First 100k
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_0.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_1.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_2.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_3.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_4.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_5.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_6.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v1/tag_1_delphes_events_7.root");

  //Next 375k
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_0.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_1.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_3.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_4.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_5.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_6.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_7.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_8.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_10.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_11.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_12.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_15.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_16.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_18.root");
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_1_delphes_events_19.root");
  //Next 400k
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v2/tag_2_delphes_events.root");
  //Next 400k
  chain.Add("root://eoscms//eos/cms/store/user/chanon/TTH/Delphes/TTZ_v3/tag_3_delphes_events.root");


  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  //Long64_t numberOfEntries = 10000;
  cout << "nEntries="<<numberOfEntries<<endl;

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  //Long64_t numberOfEntries = branchEvent->Event_size();
  //cout << "nEntries="<<numberOfEntries<<endl;

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");

  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");

  std::vector<Jet*>	  vSelectedJets;
  std::vector<Lepton*> vSelectedLeptons;
  std::vector<GenParticle*> vSelectedPartons;
  std::vector<GenParticle*> vSelectedNeutrinos;
  std::vector<GenParticle*> vSelectedGenLeptons;
  std::vector<GenParticle*> vSelectedW;
  std::vector<GenParticle*> vSelectedNeutrinosFromTopLep;
  std::vector<GenParticle*> vSelectedBquark;

  std::vector<std::pair<Jet*, GenParticle*>> MatchedBJets;
  std::vector<std::pair<Jet*, GenParticle*>> MatchedJets;

  MissingET * mET = NULL;
  ScalarHT * msumET = NULL;
 
  bool doSelectOnlyBjets = true;
  int ib1=-1, ib2=-1;
  TLorentzVector Bjet1, Bjet2;
  double Ptot = 0;
  double Ptot_Px=0, Ptot_Py=0, Ptot_Pz=0, Ptot_Eta=0, Ptot_E=0, Ptot_M=0;
  TLorentzVector Ptot_P4(0,0,0,0);

  unsigned int nSelectedEventsTTZwithGen = 0;
  unsigned int nSelectedEventsTTZ = 0;
  unsigned int nSelectedEvent = 0; 

  TFile* fOutput = new TFile("output.root","RECREATE");
  MultiLeptonTree tree;
  tree.initializeOutputTree();

  bool doTTZselection = true;
  bool doTFselection = false;
  bool storeGenLevel = true;

  //int NEND_bis = NEND;
  //if (NEND > numberOfEntries) NEND_bis = numberOfEntries; 

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if (entry % 1000 == 0) cout << "Event "<<entry<<endl;
    //cout << "Event "<<entry<<endl; 

    vSelectedJets.clear();
    vSelectedGenLeptons.clear();
    vSelectedLeptons.clear();
    vSelectedPartons.clear();
    vSelectedW.clear();
    vSelectedNeutrinosFromTopLep.clear();
    vSelectedBquark.clear();
    MatchedBJets.clear();
    MatchedJets.clear();
    vSelectedNeutrinos.clear();
    ib1=-1, ib2=-1;
    bool foundHadronicTop = false;
    int iLeptonicTop = -1;
    int iHadronicTop = -1;
    int iTop = 0;

    GenParticle *genparticle;
    Ptot_Px=0, Ptot_Py=0, Ptot_Pz=0, Ptot_Eta=0, Ptot_E=0, Ptot_M=0;
    //Checkin gen particles
    if(branchGenParticle->GetEntries() > 0) {
      for (int ig=0; ig<branchGenParticle->GetEntries(); ig++){
        genparticle = (GenParticle*) branchGenParticle->At(ig);

	if (abs(genparticle->PID)==12 || abs(genparticle->PID)==14 || abs(genparticle->PID)==16) {
	  //cout << "TEST GenParticle "<<ig<<" PID="<<genparticle->PID<<" status="<<genparticle->Status << endl;
	  int m = genparticle->M1;
          GenParticle *genmother = (GenParticle*) branchGenParticle->At(m);
	  if (abs(genmother->PID)==24){
	    vSelectedW.push_back(genmother);
            vSelectedNeutrinosFromTopLep.push_back(genparticle);
            //cout << "nNeutFromTopLep size="<<vSelectedNeutrinosFromTopLep.size()<<endl;
	    GenParticle * gentmp = genmother;
	    GenParticle *gengranma = genmother;
	    bool doGoBackHistory = true;
	    while (doGoBackHistory){
	      int gm = gentmp->M1;
              gengranma = (GenParticle*) branchGenParticle->At(gm);
	      if (abs(gengranma->PID)==6){
		doGoBackHistory = false;
	      }
	      else gentmp = gengranma;
	    }
	    int md1 = gengranma->D1;
            int md2 = gengranma->D2;
            GenParticle *genmomdaughter1 = (GenParticle*) branchGenParticle->At(md1);
            GenParticle *genmomdaughter2 = (GenParticle*) branchGenParticle->At(md2);
            if (abs(genmomdaughter1->PID)==5) vSelectedBquark.push_back(genmomdaughter1);
            else if (abs(genmomdaughter2->PID)==5) vSelectedBquark.push_back(genmomdaughter2);
	    //cout << "vSelectedBquark size="<<vSelectedBquark.size()<<endl;
	    iLeptonicTop = iTop;
	    iTop++;
	  }
	}

	if (abs(genparticle->PID)<=4 && !foundHadronicTop){
	  int m = genparticle->M1;
          GenParticle *genmother = (GenParticle*) branchGenParticle->At(m);
	  if (abs(genmother->PID)==24){
	    vSelectedW.push_back(genmother);
            GenParticle * gentmp = genmother;
            GenParticle *gengranma = genmother;
            bool doGoBackHistory = true;
            while (doGoBackHistory){
              int gm = gentmp->M1;
              gengranma = (GenParticle*) branchGenParticle->At(gm);
              if (abs(gengranma->PID)==6){
                doGoBackHistory = false;
              }
              else gentmp = gengranma;
            }
            int md1 = gengranma->D1;
            int md2 = gengranma->D2;
            GenParticle *genmomdaughter1 = (GenParticle*) branchGenParticle->At(md1);
            GenParticle *genmomdaughter2 = (GenParticle*) branchGenParticle->At(md2);
            if (abs(genmomdaughter1->PID)==5) vSelectedBquark.push_back(genmomdaughter1);
            else if (abs(genmomdaughter2->PID)==5) vSelectedBquark.push_back(genmomdaughter2);
            //cout << "vSelectedBquark size="<<vSelectedBquark.size()<<endl;
	    foundHadronicTop = true;
            iHadronicTop = iTop;
            iTop++;
	  }
	}

	if (genparticle->Status==23 || genparticle->Status==24) {
	  if (abs(genparticle->PID)<=5) {
	    vSelectedPartons.push_back(genparticle);
	    Ptot_Px += genparticle->Px;
	    Ptot_Py += genparticle->Py;
	    Ptot_Pz += genparticle->Pz;
	    Ptot_E += genparticle->E;
	    //cout << "GenParticle "<<ig<<" PID="<<genparticle->PID<<" status="<<genparticle->Status << endl;
	  }
	  if (abs(genparticle->PID)==12 || abs(genparticle->PID)==14 || abs(genparticle->PID)==16){
	    //cout << "genparticle status="<<genparticle->Status<<endl;
	    vSelectedNeutrinos.push_back(genparticle);
	    //cout << "nSelectedNeutrinos="<<vSelectedNeutrinos.size()<<endl;
	    //int m = genparticle->M1;
	    //cout << "Mother num="<<m<<endl;
            //GenParticle *genmother = (GenParticle*) branchGenParticle->At(m);
	    //cout << "Mother pid="<<genmother->PID<<" status="<< genmother->Status<<endl;
	    //vSelectedNeutrinosFromTopLep.push_back(genparticle);
	    //cout << "nNeutFromTopLep size="<<vSelectedNeutrinosFromTopLep.size()<<endl;
            Ptot_Px += genparticle->Px;
	    Ptot_Py += genparticle->Py;
	    Ptot_Pz += genparticle->Pz;
	    Ptot_E += genparticle->E;
	    //cout << "GenParticle "<<ig<<" PID="<<genparticle->PID<<" status="<<genparticle->Status << endl;
	  }
	  if (abs(genparticle->PID)==11 || abs(genparticle->PID)==13 || abs(genparticle->PID)==15){
	    vSelectedGenLeptons.push_back(genparticle);
	    Ptot_Px += genparticle->Px;
	    Ptot_Py += genparticle->Py;
	    Ptot_Pz += genparticle->Pz;
	    Ptot_E += genparticle->E;
	  }
	}
	if (genparticle->Status==22) {
	  if (abs(genparticle->PID)==24){
	    //vSelectedW.push_back(genparticle);
	    int m = genparticle->M1;
	    GenParticle *genmother = (GenParticle*) branchGenParticle->At(m);
	    int md1 = genmother->D1;
	    int md2 = genmother->D2;
	    GenParticle *genmomdaughter1 = (GenParticle*) branchGenParticle->At(md1);
	    GenParticle *genmomdaughter2 = (GenParticle*) branchGenParticle->At(md2);
	    //if (abs(genmomdaughter1->PID)==5) vSelectedBquark.push_back(genmomdaughter1); 
	    //else if (abs(genmomdaughter2->PID)==5) vSelectedBquark.push_back(genmomdaughter2);

	   /* 
	   GenParticle *gentmp = genparticle;
	   int doSearchNeutrino = 1;
	    while (doSearchNeutrino>=1){
	      //cout << "doSearchNeutrino="<<doSearchNeutrino<<endl;
	      int d1 = gentmp->D1;
	      int d2 = gentmp->D2;
	      GenParticle *gendaughter1 = (GenParticle*) branchGenParticle->At(d1);
	      GenParticle *gendaughter2 = (GenParticle*) branchGenParticle->At(d2);
	      if (gendaughter1->Status==23 || gendaughter2->Status==23){
	        //cout << "Wdaughter1 pid="<<gendaughter1->PID<<" status="<<gendaughter1->Status<<endl;
                //cout << "Wdaughter2 pid="<<gendaughter2->PID<<" status="<<gendaughter2->Status<<endl;
	        if (abs(gendaughter1->PID)==12 || abs(gendaughter1->PID)==14 || abs(gendaughter1->PID)==16) { vSelectedNeutrinosFromTopLep.push_back(gendaughter1); doSearchNeutrino=0; }
                else if (abs(gendaughter2->PID)==12 || abs(gendaughter2->PID)==14 || abs(gendaughter2->PID)==16) { vSelectedNeutrinosFromTopLep.push_back(gendaughter2); doSearchNeutrino=0; }
	        else doSearchNeutrino++;
	        gentmp = gendaughter1;
		//cout << "vNeutFromTopLep="<<vSelectedNeutrinosFromTopLep.size()<<endl;
	      }
	      else doSearchNeutrino = 0;
	    }
	    */ 
	  }
	}
      }
    }
    Ptot_P4.SetPxPyPzE(Ptot_Px, Ptot_Py, Ptot_Pz, Ptot_E);
    Ptot = sqrt(Ptot_Px*Ptot_Px+Ptot_Py*Ptot_Py);
    //cout << "Ptot="<<Ptot<<endl;
 
    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0) {

      for (int ij=0; ij<branchJet->GetEntries(); ij++){
	Jet *jet = (Jet*) branchJet->At(ij);
	if (jet->PT > 25 && fabs(jet->Eta)<2.4) {
	  vSelectedJets.push_back(jet);
	  //cout << "Selected jet PT="<<jet->PT<<endl;
	}
      }

      ib1=-1, ib2=-1;
      selectBjets(vSelectedJets, "BtagHighestPt", &ib1, &ib2, doSelectOnlyBjets);
      //cout << "ib1="<<ib1<< " ib2="<<ib2<<endl; 

    }

    //Jet matching
    if (vSelectedJets.size()>0 && vSelectedPartons.size()>0){
      for (unsigned int i=0; i<vSelectedJets.size(); i++){
	TLorentzVector Pjet; Pjet.SetPtEtaPhiM(vSelectedJets.at(i)->PT, vSelectedJets.at(i)->Eta, vSelectedJets.at(i)->Phi, vSelectedJets.at(i)->Mass);

	for (unsigned int j=0; j<vSelectedPartons.size(); j++){
	  TLorentzVector Ppart; Ppart.SetPtEtaPhiM(vSelectedPartons.at(j)->PT, vSelectedPartons.at(j)->Eta, vSelectedPartons.at(j)->Phi, vSelectedPartons.at(j)->Mass);
	  double deltaR = Pjet.DeltaR(Ppart);
	  if (deltaR<0.2) {
	    //if (abs(vSelectedPartons.at(j)->PID)==5 && vSelectedJets.at(i)->BTag==1) MatchedBJets.push_back(std::make_pair(vSelectedJets.at(i), vSelectedPartons.at(j)));
	    if (abs(vSelectedPartons.at(j)->PID)<5) MatchedJets.push_back(std::make_pair(vSelectedJets.at(i), vSelectedPartons.at(j)));
	    //cout << "Jet"<<i<<" Parton"<<j<<" deltaR="<<deltaR<<endl;
	  }
	}
	for (unsigned int j=0; j<vSelectedBquark.size(); j++){
	  TLorentzVector Ppart; Ppart.SetPtEtaPhiM(vSelectedBquark.at(j)->PT, vSelectedBquark.at(j)->Eta, vSelectedBquark.at(j)->Phi, vSelectedBquark.at(j)->Mass);
	  double deltaR = Pjet.DeltaR(Ppart);
	  if (deltaR<0.2) {
	    if (abs(vSelectedBquark.at(j)->PID)==5 && vSelectedJets.at(i)->BTag==1) MatchedBJets.push_back(std::make_pair(vSelectedJets.at(i), vSelectedBquark.at(j)));
	  }
	}
      }
    }


    //Selecting leptons

    if (branchElectron->GetEntries() > 0) {
      for (int ie=0; ie<branchElectron->GetEntries(); ie++){
        Lepton * elec = new Lepton((Electron*)branchElectron->At(ie));//static_cast<Candidate*>( branchElectron->At(ie) );
	if (elec->PT > 15 && fabs(elec->Eta)<2.4){
	  vSelectedLeptons.push_back(elec);
	  //cout << "Selected Electron pT="<<elec->PT<<endl;
	}
      }
    }

    if (branchMuon->GetEntries() > 0) {
      for (int ie=0; ie<branchMuon->GetEntries(); ie++){
        Lepton * muon = new Lepton((Muon*)branchMuon->At(ie));//static_cast<Candidate*>( branchElectron->At(ie) );
        if (muon->PT > 15 && fabs(muon->Eta)<2.4){
          vSelectedLeptons.push_back(muon);
          //cout << "Selected Muon pT="<<muon->PT<<endl;
        }
      }
    }

    //cout << "Lepton size="<<vSelectedLeptons.size()<<endl;
    std::sort(vSelectedLeptons.begin(), vSelectedLeptons.end(), SortingLeptonPt);

    if (branchMissingET->GetEntries() == 1){
        mET = (MissingET*) branchMissingET->At(0);
    }

    if (branchScalarHT->GetEntries() == 1){
	msumET = (ScalarHT*) branchScalarHT->At(0);
    }

    
    int tot_charge = 0;
    int tot_id = 0;
    if (vSelectedLeptons.size()>=4){
      for (unsigned int i=0; i<vSelectedLeptons.size(); i++) {
        tot_charge += vSelectedLeptons.at(i)->Charge;
        tot_id += vSelectedLeptons.at(i)->Id;
      }
    }

    bool is3l = false;
    if (vSelectedLeptons.size()==3) is3l = true;

    //Mll > 12
    bool pass_invmasscut = true;
    if (vSelectedLeptons.size()>=2){
      for(unsigned int i=0; i<vSelectedLeptons.size()-1; i++)
      {
        TLorentzVector P1; P1.SetPtEtaPhiM(vSelectedLeptons.at(i)->PT, vSelectedLeptons.at(i)->Eta, vSelectedLeptons.at(i)->Phi, 0);
        for(unsigned int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            TLorentzVector P2; P2.SetPtEtaPhiM(vSelectedLeptons.at(j)->PT, vSelectedLeptons.at(j)->Eta, vSelectedLeptons.at(j)->Phi, 0);
            if ( fabs( ( P1 + P2).M() ) < 12 )
            { pass_invmasscut = false ;}
        }
      } 
    }

    //SFOS pair
    bool isSFOS = false;
    for(unsigned int i=0; i<vSelectedLeptons.size(); i++) {
        for(unsigned int j=0; j<vSelectedLeptons.size(); j++) {
            if (  ( i                           != j                            )
                    && ( vSelectedLeptons.at(i)->Id == -vSelectedLeptons.at(j)->Id ) )
            { isSFOS = true ;}
        }
    }

    //Z veto 
    bool pass_Zveto = true;
    if (vSelectedLeptons.size()>=2){
      for(unsigned int i=0; i<vSelectedLeptons.size()-1; i++)
      {
        TLorentzVector P1; P1.SetPtEtaPhiM(vSelectedLeptons.at(i)->PT, vSelectedLeptons.at(i)->Eta, vSelectedLeptons.at(i)->Phi, 0);
        for(unsigned int j=i+1; j<vSelectedLeptons.size(); j++)
        {
            TLorentzVector P2; P2.SetPtEtaPhiM(vSelectedLeptons.at(j)->PT, vSelectedLeptons.at(j)->Eta, vSelectedLeptons.at(j)->Phi, 0);
            if (  ( vSelectedLeptons.at(i)->Id     == -vSelectedLeptons.at(j)->Id                             )
               && ( vSelectedLeptons.at(i)->Charge == -vSelectedLeptons.at(j)->Charge                         )
                    && ( fabs( ( P1+P2).M() - 91.188) < 10 ) )
            { pass_Zveto = false ;}
        }
      }
    }

    //MET LD
    float jet_px = 0, jet_py = 0, lepton_px = 0, lepton_py = 0, tau_px = 0, tau_py = 0, MHT = 0, met_ld = 0;
    TLorentzVector jetp4;
    for(unsigned int i=0; i<vSelectedJets.size(); i++)
    {
        jetp4.SetPtEtaPhiM(vSelectedJets.at(i)->PT, vSelectedJets.at(i)->Eta, vSelectedJets.at(i)->Phi, vSelectedJets.at(i)->Mass);
        jet_px = jet_px + jetp4.Px();
        jet_py = jet_py + jetp4.Py();
    }
    for(unsigned int i=0; i<vSelectedLeptons.size(); i++)
    {
        TLorentzVector P1; P1.SetPtEtaPhiM(vSelectedLeptons.at(i)->PT, vSelectedLeptons.at(i)->Eta, vSelectedLeptons.at(i)->Phi, 0);
        lepton_px = lepton_px + P1.Px();
        lepton_py = lepton_py + P1.Py();
    }
    MHT = sqrt( (jet_px + lepton_px + tau_px) * (jet_px + lepton_px + tau_px) + (jet_py + lepton_py + tau_py) * (jet_py + lepton_py + tau_py) );
    met_ld = 0.00397 * mET->MET + 0.00265 * MHT;



    if (doTTZselection){

      if (vSelectedJets.size()<2) continue;
      if (ib1==-1) continue;
      //if (ib2==-1) continue;

      if (vSelectedLeptons.size()<3) continue;
      //if (vSelectedLeptons.at(0)->PT<25) continue;

      //if(!pass_invmasscut) continue;

      if (!isSFOS) continue;
      tree.mc_ttZhypAllowed = 1;

      //if(pass_Zveto)       continue;
      tree.is_3l_TTZ_CR = true;

      if (vSelectedLeptons.at(0)->Charge*vSelectedLeptons.at(1)->Charge>0 && vSelectedLeptons.at(1)->Charge*vSelectedLeptons.at(2)->Charge>0) continue;

      nSelectedEventsTTZ++;

      //if(vSelectedJets.size() < 4 && (met_ld < (0.2 + 0.1 * isSFOS)) ) continue;

      if(storeGenLevel){
	
	//if (MatchedJets.size()!=2) continue;
	if (MatchedBJets.size()!=2) continue;
        //if (vSelectedJets.size()!=2) continue;
        //if (vSelectedBJets.size()!=2) continue;
        //if (vSelectedNeutrinos.size()>=1) continue;
	if (vSelectedNeutrinosFromTopLep.size()!=1) continue;
	if (vSelectedW.size()!=2) continue; 
        /*
        for (int ij=0; ij< MatchedJets.size(); ij++){
          TLorentzVector Pjet; Pjet.SetPtEtaPhiM(MatchedJets.at(ij).first->PT, MatchedJets.at(ij).first->Eta, MatchedJets.at(ij).first->Phi, MatchedJets.at(ij).first->Mass);
          TLorentzVector Ppart; Ppart.SetPtEtaPhiM(MatchedJets.at(ij).second->PT, MatchedJets.at(ij).second->Eta, MatchedJets.at(ij).second->Phi, MatchedJets.at(ij).second->Mass);
	  if (ij==0) {
            tree.multilepton_Jet1_DeltaR_Matched = Pjet.DeltaR(Ppart);
            tree.multilepton_Jet1_Id_Matched = 0;   
            tree.multilepton_Jet1_P4_Matched = Ppart;
	  }
          if (ij==1) {
            tree.multilepton_Jet2_DeltaR_Matched = Pjet.DeltaR(Ppart);
            tree.multilepton_Jet2_Id_Matched = 0;   
            tree.multilepton_Jet2_P4_Matched = Ppart;
          }
        }
        */
        for (int ij=0; ij< MatchedBJets.size(); ij++){
          TLorentzVector Pjet; Pjet.SetPtEtaPhiM(MatchedBJets.at(ij).first->PT, MatchedBJets.at(ij).first->Eta, MatchedBJets.at(ij).first->Phi, MatchedBJets.at(ij).first->Mass);
          TLorentzVector Ppart; Ppart.SetPtEtaPhiM(MatchedBJets.at(ij).second->PT, MatchedBJets.at(ij).second->Eta, MatchedBJets.at(ij).second->Phi, MatchedBJets.at(ij).second->Mass);
	  if (ij==iHadronicTop){
            tree.multilepton_Bjet1_DeltaR_Matched = Pjet.DeltaR(Ppart); 
            tree.multilepton_Bjet1_Id_Matched = 5;       
            tree.multilepton_Bjet1_P4_Matched = Ppart;       
	  }
          if (ij==iLeptonicTop){
            tree.multilepton_Bjet2_DeltaR_Matched = Pjet.DeltaR(Ppart);
            tree.multilepton_Bjet2_Id_Matched = 5;
            tree.multilepton_Bjet2_P4_Matched = Ppart;
          }
        }

	//cout << "nNeutrinosFromTopLep: "<<vSelectedNeutrinosFromTopLep.size()<<endl;
	//cout << "nW: "<<vSelectedW.size()<<endl;
	//cout << "nBquarks: " << vSelectedBquark.size() << " nBjetMatched="<< MatchedBJets.size()<<endl;

	if (vSelectedNeutrinosFromTopLep.size()>0){
          TLorentzVector Pneutrino; Pneutrino.SetPtEtaPhiM(vSelectedNeutrinosFromTopLep.at(0)->PT, vSelectedNeutrinosFromTopLep.at(0)->Eta, vSelectedNeutrinosFromTopLep.at(0)->Phi, 0);
	  tree.multilepton_mET_Matched = Pneutrino;
	}

	if (vSelectedW.size()>=2){
	  TLorentzVector PW1; PW1.SetPtEtaPhiM(vSelectedW.at(iHadronicTop)->PT, vSelectedW.at(iHadronicTop)->Eta, vSelectedW.at(iHadronicTop)->Phi, vSelectedW.at(iHadronicTop)->Mass);
          TLorentzVector PW2; PW2.SetPtEtaPhiM(vSelectedW.at(iLeptonicTop)->PT, vSelectedW.at(iLeptonicTop)->Eta, vSelectedW.at(iLeptonicTop)->Phi, vSelectedW.at(iLeptonicTop)->Mass); 
	  tree.multilepton_W1_P4_Matched = PW1;
          tree.multilepton_W2_P4_Matched = PW2;
	}
	
	nSelectedEventsTTZwithGen++;
      }
    }

    if (doTFselection){

      if (vSelectedJets.size()<1) continue;
      //TRefArray Constituents = vSelectedJets.at(0)->Constituents;
      //for (int i=0; i<Constituents.GetEntries(); i++) {
	//cout << "Constituents "<<i<< " " <<Constituents.At(i)<<endl;
      //}
      //TRefArray Particles = vSelectedJets.at(0)->Particles;
      //cout << "First jet Particles size="<<Particles.GetEntriesFast()<<endl;
      //for (int i=0; i<Particles.GetEntriesFast(); i++) {
	//GenParticle* particle = (GenParticle*) Particles.At(i);
	//cout << "Particle "<< i<<" Id=" <<particle->PID<<endl;
      //}

      //Analyze matched jets TF
      for (int ij=0; ij< MatchedJets.size(); ij++){
        TLorentzVector Pjet; Pjet.SetPtEtaPhiM(MatchedJets.at(ij).first->PT, MatchedJets.at(ij).first->Eta, MatchedJets.at(ij).first->Phi, MatchedJets.at(ij).first->Mass);
        TLorentzVector Ppart; Ppart.SetPtEtaPhiM(MatchedJets.at(ij).second->PT, MatchedJets.at(ij).second->Eta, MatchedJets.at(ij).second->Phi, MatchedJets.at(ij).second->Mass);
	double Erec = Pjet.E();
	double Epart = Ppart.E();
	//tree.TFratio_nonB->Fill(Erec/Epart);
        //cout << "Jets Erec/Epart="<<Erec/Epart<<endl;
	int iEta = GetTransferFunctionEta(Pjet.Eta());
	int iEnergy = GetTransferFunctionEnergy(Pjet.E());
	tree.TFratio_nonB[iEta][iEnergy]->Fill(Erec/Epart);
      }
      for (int ij=0; ij< MatchedBJets.size(); ij++){
        TLorentzVector Pjet; Pjet.SetPtEtaPhiM(MatchedBJets.at(ij).first->PT, MatchedBJets.at(ij).first->Eta, MatchedBJets.at(ij).first->Phi, MatchedBJets.at(ij).first->Mass);
        TLorentzVector Ppart; Ppart.SetPtEtaPhiM(MatchedBJets.at(ij).second->PT, MatchedBJets.at(ij).second->Eta, MatchedBJets.at(ij).second->Phi, MatchedBJets.at(ij).second->Mass);
        double Erec = Pjet.E();
        double Epart = Ppart.E();
        //tree.TFratio_B->Fill(Erec/Epart);
	//cout << "Bjets Erec/Epart="<<Erec/Epart<<endl;
        int iEta = GetTransferFunctionEta(Pjet.Eta());
        int iEnergy = GetTransferFunctionEnergy(Pjet.E());
        tree.TFratio_B[iEta][iEnergy]->Fill(Erec/Epart);
      }

      //Missing ET TF
      TLorentzVector GenNeutrinoSum; GenNeutrinoSum.SetPtEtaPhiM(0,0,0,0);
      for (int in=0; in<vSelectedNeutrinos.size(); in++){
	TLorentzVector Pneutrino; Pneutrino.SetPtEtaPhiM(vSelectedNeutrinos.at(in)->PT, vSelectedNeutrinos.at(in)->Eta, vSelectedNeutrinos.at(in)->Phi, 0);
	GenNeutrinoSum += Pneutrino;
      }
      double GenMet_Px = GenNeutrinoSum.Px();
      double GenMet_Py = GenNeutrinoSum.Py();
      double GenMet_Pt = GenNeutrinoSum.Pt();
      double GenMet_Phi = GenNeutrinoSum.Phi();
      double RecMet_Px = mET->MET*cos(mET->Phi);
      double RecMet_Py = mET->MET*sin(mET->Phi);
      double RecMet_Pt = mET->MET;
      double RecMet_Phi = mET->Phi;
      //cout << "Gen MET Px="<<GenMet_Px<<" Py="<<GenMet_Py<<" Pt="<<GenMet_Pt<<" GenMet_Phi="<<GenMet_Phi<<endl;
      //cout << "Rec MET Px="<<RecMet_Px<<" Py="<<RecMet_Py<<" Pt="<<RecMet_Pt<<" RecMet_Phi="<<RecMet_Phi<<endl;
      if (GenMet_Pt>0){
        int iMet = GetTransferFunctionMet(mET->MET, msumET->HT);
        tree.TF_mET_Px[iMet]->Fill(GenMet_Px-RecMet_Px);
        tree.TF_mET_Py[iMet]->Fill(GenMet_Py-RecMet_Py);
        tree.TF_mET_Pt[iMet]->Fill(GenMet_Pt-RecMet_Pt);
        tree.TF_mET_Phi[iMet]->Fill(GenMet_Phi-RecMet_Phi);
      }
    }


    //Fill Output Tree

    tree.Pdf_Ptot->Fill(Ptot);
    tree.Pdf_PtotEta->Fill(Ptot_P4.Eta());
    tree.Pdf_PtotM->Fill(Ptot_P4.M());
    tree.Pdf_PtotP4->Fill(Ptot_P4.Pt(), Ptot_P4.Eta(), Ptot_P4.M());

  if (doTTZselection){ 
    if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2>=2) tree.catJets = kCat_3l_2b_2j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1>=2) tree.catJets = kCat_3l_1b_2j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==1) tree.catJets = kCat_3l_2b_1j;
    else if (is3l && ib1!=-1 && ib2==-1 && vSelectedJets.size()-1==1) tree.catJets = kCat_3l_1b_1j;
    else if (is3l && ib1!=-1 && ib2!=-1 && vSelectedJets.size()-2==0) tree.catJets = kCat_3l_2b_0j;
    else tree.catJets = -1;

    if (tree.catJets != kCat_3l_2b_2j) continue;

    tree.multilepton_Lepton1_Id = -999;
    tree.multilepton_Lepton2_Id = -999;
    tree.multilepton_Lepton3_Id = -999;
    tree.multilepton_Lepton4_Id = -999;

    if (vSelectedLeptons.size()>=2){
        tree.multilepton_Lepton1_P4.SetPtEtaPhiM(vSelectedLeptons.at(0)->PT, vSelectedLeptons.at(0)->Eta, vSelectedLeptons.at(0)->Phi, 0);//vSelectedLeptons.at(0).p4();
        tree.multilepton_Lepton1_Id = vSelectedLeptons.at(0)->Id;
        tree.multilepton_Lepton2_P4.SetPtEtaPhiM(vSelectedLeptons.at(1)->PT, vSelectedLeptons.at(1)->Eta, vSelectedLeptons.at(1)->Phi, 0);
        tree.multilepton_Lepton2_Id = vSelectedLeptons.at(1)->Id;
    }

    if (vSelectedLeptons.size()>=3)
    {
        tree.multilepton_Lepton3_P4.SetPtEtaPhiM(vSelectedLeptons.at(2)->PT, vSelectedLeptons.at(2)->Eta, vSelectedLeptons.at(2)->Phi, 0);
        tree.multilepton_Lepton3_Id = vSelectedLeptons.at(2)->Id;
    }

    //if (vSelectedLeptons.size()>=4 && is4l)
    //{
    //    tree.multilepton_Lepton4_P4.SetPtEtaPhiM(vSelectedLeptons.at(3)->PT, vSelectedLeptons.at(3)->Eta, vSelectedLeptons.at(3)->Phi, 0);
    //    tree.multilepton_Lepton4_Id = vSelectedLeptons.at(3)->Id;
    //}

   tree.multilepton_Bjet1_Id = -999;
   tree.multilepton_Bjet2_Id = -999;

   if (ib1!=-1) tree.FillJetInfoOutputTree(&tree.multilepton_Bjet1_Id, 5, &tree.multilepton_Bjet1_P4, vSelectedJets.at(ib1)); 
   if (ib2!=-1) tree.FillJetInfoOutputTree(&tree.multilepton_Bjet2_Id, 5, &tree.multilepton_Bjet2_P4, vSelectedJets.at(ib2));

    tree.multilepton_JetHighestPt1_Id = -999;
    tree.multilepton_JetHighestPt2_Id = -999;
    tree.multilepton_JetClosestMw1_Id = -999;
    tree.multilepton_JetClosestMw2_Id = -999;
    tree.multilepton_JetLowestMjj1_Id = -999;
    tree.multilepton_JetLowestMjj2_Id = -999;

    tree.multilepton_JetHighestPt1_2ndPair_Id = -999;
    tree.multilepton_JetHighestPt2_2ndPair_Id = -999;
    tree.multilepton_JetClosestMw1_2ndPair_Id = -999;
    tree.multilepton_JetClosestMw2_2ndPair_Id = -999;
    tree.multilepton_JetLowestMjj1_2ndPair_Id = -999;
    tree.multilepton_JetLowestMjj2_2ndPair_Id = -999;

    //Choose 2 jets
    TLorentzVector Pjet1, Pjet2;
    float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
    float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
    for (unsigned int ij=0; ij<vSelectedJets.size(); ij++){
        if (ij==ib1 || ij==ib2) continue;
        if (vSelectedJets.at(ij)->PT > pt_max ) {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = vSelectedJets.at(ij)->PT;
            ij1 = ij;
        } 
        if (vSelectedJets.at(ij)->PT < pt_max && vSelectedJets.at(ij)->PT > pt_max2){
            pt_max2 = vSelectedJets.at(ij)->PT; 
            ij2 = ij; 
        } 
        for (unsigned int ik=0; ik<vSelectedJets.size(); ik++){
            if (ik==ij) continue;
            if (ik==ib1 || ik==ib2) continue;
            Pjet1.SetPtEtaPhiM(vSelectedJets.at(ij)->PT, vSelectedJets.at(ij)->Eta, vSelectedJets.at(ij)->Phi, vSelectedJets.at(ij)->Mass);
            Pjet2.SetPtEtaPhiM(vSelectedJets.at(ik)->PT, vSelectedJets.at(ik)->Eta, vSelectedJets.at(ik)->Phi, vSelectedJets.at(ik)->Mass); 
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

    if (ij1!=-1) tree.FillJetInfoOutputTree(&tree.multilepton_JetHighestPt1_Id, 1, &tree.multilepton_JetHighestPt1_P4, vSelectedJets.at(ij1)); 
    if (ij2!=-1) tree.FillJetInfoOutputTree(&tree.multilepton_JetHighestPt2_Id, 1, &tree.multilepton_JetHighestPt2_P4, vSelectedJets.at(ij2));
    if (ik1!=-1 && ik2!=-1){
      tree.FillJetInfoOutputTree(&tree.multilepton_JetClosestMw1_Id, 2, &tree.multilepton_JetClosestMw1_P4, vSelectedJets.at(ik1));
      tree.FillJetInfoOutputTree(&tree.multilepton_JetClosestMw2_Id, 2, &tree.multilepton_JetClosestMw2_P4, vSelectedJets.at(ik2));
    }
    if (il1!=-1 && il2!=-1){
      tree.FillJetInfoOutputTree(&tree.multilepton_JetLowestMjj1_Id, 3, &tree.multilepton_JetLowestMjj1_P4, vSelectedJets.at(il1));
      tree.FillJetInfoOutputTree(&tree.multilepton_JetLowestMjj2_Id, 3, &tree.multilepton_JetLowestMjj2_P4, vSelectedJets.at(il2));
    }
    /*
    //Choose 2 more jets
    int io1=-1, io2=-1, ip1=-1, ip2=-1, im1=-1, im2=-1;
    diffmass_min = 10000, mass_min = 10000, pt_max2 = 0, pt_max = 0;
    for (unsigned int im=0; im<vSelectedJets.size(); im++){
        if (im==ib1 || im==ib2 || im==ik1 || im==ik2) continue;
        if (vSelectedJets.at(im)->PT > pt_max ) {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = vSelectedJets.at(im)->PT;
            im1 = im;
        }
        if (vSelectedJets.at(im)->PT < pt_max && vSelectedJets.at(im)->PT > pt_max2){
            pt_max2 = vSelectedJets.at(im)->PT;
            im2 = im;
        }
        for (unsigned int in=0; in<vSelectedJets.size(); in++){
            if (in==ib1 || in==ib2 || in==ik1 || in==ik2 || in==im) continue;
            Pjet1.SetPtEtaPhiM(vSelectedJets.at(im)->PT, vSelectedJets.at(im)->Eta, vSelectedJets.at(im)->Phi, vSelectedJets.at(im)->Mass);
            Pjet2.SetPtEtaPhiM(vSelectedJets.at(in)->PT, vSelectedJets.at(in)->Eta, vSelectedJets.at(in)->Phi, vSelectedJets.at(in)->Mass);
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
    */

    tree.multilepton_mET.SetPtEtaPhiE(mET->MET, 0, mET->Phi, mET->MET);
    tree.multilepton_mHT = msumET->HT;

    tree.tOutput->Fill();

  }

    nSelectedEvent++;

  }
  
  cout << "---RESULTS---"<<endl;
  cout << "nSelectedEvent="<<nSelectedEvent<<endl;
  cout << "nSelectedEventsTTZ="<<nSelectedEventsTTZ<<endl;
  cout << "nSelectedEventsTTZwithGen="<<nSelectedEventsTTZwithGen<<endl;

   double a = 0;
   a  = tree.Pdf_Ptot->Integral();
   tree.Pdf_Ptot->Scale(1./a);
   tree.Pdf_Ptot->Write();
 
   a  = tree.Pdf_PtotEta->Integral();
   tree.Pdf_PtotEta->Scale(1./a);
   tree.Pdf_PtotEta->Write();

   a  = tree.Pdf_PtotM->Integral();
   tree.Pdf_PtotM->Scale(1./a);
   tree.Pdf_PtotM->Write();
 
   a  = tree.Pdf_PtotP4->Integral();
   tree.Pdf_PtotP4->Scale(1./a);
   tree.Pdf_PtotP4->Write();

    for (int ieta=0; ieta<3; ieta++)  {
      for (int ienergy=0; ienergy<6; ienergy++) {
	a = tree.TF_B[ieta][ienergy]->Integral();
	if (a>0) tree.TF_B[ieta][ienergy]->Scale(1./a);
        tree.TF_B[ieta][ienergy]->Write();

	a = tree.TF_nonB[ieta][ienergy]->Integral();
	if (a>0) tree.TF_nonB[ieta][ienergy]->Scale(1./a);
        tree.TF_nonB[ieta][ienergy]->Write();
	
	a = tree.TFratio_B[ieta][ienergy]->Integral();
	if (a>0) tree.TFratio_B[ieta][ienergy]->Scale(1./a);
        tree.TFratio_B[ieta][ienergy]->Write();

	a = tree.TFratio_nonB[ieta][ienergy]->Integral();
	if (a>0) tree.TFratio_nonB[ieta][ienergy]->Scale(1./a);
        tree.TFratio_nonB[ieta][ienergy]->Write();
      }
    }
    for (int imet=0; imet<6; imet++){
      a = tree.TF_mET_Px[imet]->Integral();
      if (a>0) tree.TF_mET_Px[imet]->Scale(1./a); 
      tree.TF_mET_Px[imet]->Write();

      a = tree.TF_mET_Py[imet]->Integral();
      if (a>0) tree.TF_mET_Py[imet]->Scale(1./a);
      tree.TF_mET_Py[imet]->Write();

      a = tree.TF_mET_Pt[imet]->Integral();
      if (a>0) tree.TF_mET_Pt[imet]->Scale(1./a); 
      tree.TF_mET_Pt[imet]->Write();

      a = tree.TF_mET_Phi[imet]->Integral();
      if (a>0) tree.TF_mET_Phi[imet]->Scale(1./a);
      tree.TF_mET_Phi[imet]->Write();
    }


   tree.tOutput->Write();
   fOutput->Close();
}




