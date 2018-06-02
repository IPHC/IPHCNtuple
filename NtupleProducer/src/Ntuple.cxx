#include "include/Ntuple.h"

Ntuple::Ntuple(std::string name)
{
   fname_out = name;
}

Ntuple::~Ntuple()
{   
   m_file->Write();
   m_file->Close();
}

void Ntuple::Init()
{
   m_file = new TFile(fname_out.c_str(),"RECREATE");
   m_tree = new TTree("Nt","Ntuple");
}

void Ntuple::setBranchAddress()
{
   m_tree->Branch("Event",             "std::vector<Event>",       (NtEvent),32000,1);
   m_tree->Branch("ElectronLoose",     "std::vector<Electron>",    (NtElectronLoose),32000,1);
   m_tree->Branch("ElectronFakeable",  "std::vector<Electron>",    (NtElectronFakeable),32000,1);
   m_tree->Branch("ElectronTight",     "std::vector<Electron>",    (NtElectronTight),32000,1);
   m_tree->Branch("MuonLoose",         "std::vector<Muon>",        (NtMuonLoose),32000,1);
   m_tree->Branch("MuonFakeable",      "std::vector<Muon>",        (NtMuonFakeable),32000,1);
   m_tree->Branch("MuonTight",         "std::vector<Muon>",        (NtMuonTight),32000,1);
   m_tree->Branch("TauFakeable",       "std::vector<Tau>",         (NtTauFakeable),32000,1);
   m_tree->Branch("TauLoose",          "std::vector<Tau>",         (NtTauLoose),32000,1);
   m_tree->Branch("TauMedium",         "std::vector<Tau>",         (NtTauMedium),32000,1);
   m_tree->Branch("JetLoose",          "std::vector<Jet>",         (NtJetLoose),32000,1);
   m_tree->Branch("Truth",             "std::vector<Truth>",       (NtTruth),32000,1);
   m_tree->Branch("GenJet",            "std::vector<GenJet>",      (NtGenJet),32000,1);
   m_tree->Branch("TriggerObj",        "std::vector<TriggerObj>",  (NtTriggerObj),32000,1);
}

void Ntuple::createVar()
{
   NtEvent              = new std::vector<Event>;
   NtElectronLoose      = new std::vector<Electron>;
   NtElectronFakeable   = new std::vector<Electron>;
   NtElectronTight      = new std::vector<Electron>;
   NtMuonLoose          = new std::vector<Muon>;
   NtMuonFakeable       = new std::vector<Muon>;
   NtMuonTight          = new std::vector<Muon>;
   NtTauFakeable        = new std::vector<Tau>;
   NtTauLoose           = new std::vector<Tau>;
   NtTauMedium          = new std::vector<Tau>;
   NtJetLoose           = new std::vector<Jet>;
   NtTruth              = new std::vector<Truth>;
   NtGenJet             = new std::vector<GenJet>; 
   NtTriggerObj         = new std::vector<TriggerObj>;

   NtEventExt              = new std::vector<EventExt>;
   NtElectronLooseExt      = new std::vector<ElectronExt>;
   NtElectronFakeableExt   = new std::vector<ElectronExt>;
   NtElectronTightExt      = new std::vector<ElectronExt>;
   NtMuonLooseExt          = new std::vector<MuonExt>;
   NtMuonFakeableExt       = new std::vector<MuonExt>;
   NtMuonTightExt          = new std::vector<MuonExt>;
   NtTauFakeableExt        = new std::vector<TauExt>;
   NtTauLooseExt           = new std::vector<TauExt>;
   NtTauMediumExt          = new std::vector<TauExt>;
   NtJetLooseExt           = new std::vector<JetExt>;
   NtTruthExt              = new std::vector<TruthExt>;
   NtGenJetExt             = new std::vector<GenJetExt>; 
   NtTriggerObjExt         = new std::vector<TriggerObjExt>;
}

void Ntuple::clearVar()
{
   NtEvent->clear();
   NtElectronLoose->clear();
   NtElectronFakeable->clear();
   NtElectronTight->clear();
   NtMuonLoose->clear();
   NtMuonFakeable->clear();
   NtMuonTight->clear();
   NtTauFakeable->clear();
   NtTauLoose->clear();
   NtTauMedium->clear();
   NtJetLoose->clear();
   NtTruth->clear();
   NtGenJet->clear(); 
   NtTriggerObj->clear();

   NtEventExt->clear();
   NtElectronLooseExt->clear();
   NtElectronFakeableExt->clear();
   NtElectronTightExt->clear();
   NtMuonLooseExt->clear();
   NtMuonFakeableExt->clear();
   NtMuonTightExt->clear();
   NtTauFakeableExt->clear();
   NtTauLooseExt->clear();
   NtTauMediumExt->clear();
   NtJetLooseExt->clear();
   NtTruthExt->clear();
   NtGenJetExt->clear(); 
   NtTriggerObjExt->clear();
}

void Ntuple::fill()
{
   convert();
}

void Ntuple::write()
{
   m_tree->Fill();
}

void Ntuple::convert()
{
   for(int i=0;i<NtEventExt->size();i++)
     {
	Event* event = dynamic_cast<Event*>(&(NtEventExt->at(i)));
	NtEvent->push_back(*event);
     }

   for(int i=0;i<NtElectronLooseExt->size();i++)
     {
	Electron* electron = dynamic_cast<Electron*>(&(NtElectronLooseExt->at(i)));
	NtElectronLoose->push_back(*electron);
     }
   for(int i=0;i<NtElectronFakeableExt->size();i++)
     {
	Electron* electron = dynamic_cast<Electron*>(&(NtElectronFakeableExt->at(i)));
	NtElectronFakeable->push_back(*electron);
     }
   for(int i=0;i<NtElectronTightExt->size();i++)
     {
	Electron* electron = dynamic_cast<Electron*>(&(NtElectronTightExt->at(i)));
	NtElectronTight->push_back(*electron);
     }
   
   for(int i=0;i<NtMuonLooseExt->size();i++)
     {
	Muon* muon = dynamic_cast<Muon*>(&(NtMuonLooseExt->at(i)));
	NtMuonLoose->push_back(*muon);
     }
   for(int i=0;i<NtMuonFakeableExt->size();i++)
     {
	Muon* muon = dynamic_cast<Muon*>(&(NtMuonFakeableExt->at(i)));
	NtMuonFakeable->push_back(*muon);
     }
   for(int i=0;i<NtMuonTightExt->size();i++)
     {
	Muon* muon = dynamic_cast<Muon*>(&(NtMuonTightExt->at(i)));
	NtMuonTight->push_back(*muon);
     }

   for(int i=0;i<NtTauFakeableExt->size();i++)
     {
	Tau* tau = dynamic_cast<Tau*>(&(NtTauFakeableExt->at(i)));
	NtTauFakeable->push_back(*tau);
     }
   for(int i=0;i<NtTauLooseExt->size();i++)
     {
	Tau* tau = dynamic_cast<Tau*>(&(NtTauLooseExt->at(i)));
	NtTauLoose->push_back(*tau);
     }
   for(int i=0;i<NtTauMediumExt->size();i++)
     {
	Tau* tau = dynamic_cast<Tau*>(&(NtTauMediumExt->at(i)));
	NtTauMedium->push_back(*tau);
     }

   for(int i=0;i<NtJetLooseExt->size();i++)
     {
	Jet* jet = dynamic_cast<Jet*>(&(NtJetLooseExt->at(i)));
	NtJetLoose->push_back(*jet);
     }

   for(int i=0;i<NtTruthExt->size();i++)
     {
	Truth* truth = dynamic_cast<Truth*>(&(NtTruthExt->at(i)));
	NtTruth->push_back(*truth);
     }

   for(int i=0;i<NtGenJetExt->size();i++)
     {
	GenJet* genjet = dynamic_cast<GenJet*>(&(NtGenJetExt->at(i)));
	NtGenJet->push_back(*genjet);
     }

   for(int i=0;i<NtTriggerObjExt->size();i++)
     {
	TriggerObj* triggerobj = dynamic_cast<TriggerObj*>(&(NtTriggerObjExt->at(i)));
	NtTriggerObj->push_back(*triggerobj);
     }
}
