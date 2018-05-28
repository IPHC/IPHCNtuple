#include "include/Ntuple.h"

Ntuple::Ntuple(std::string fname_out)
{
   _fname_out = fname_out;
}

Ntuple::~Ntuple()
{   
   m_file->Write();
   m_file->Close();
}

void Ntuple::Init()
{
   m_file = new TFile(_fname_out.c_str(),"RECREATE");
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
   m_tree->Branch("TauTight",          "std::vector<Tau>",         (NtTauTight),32000,1);
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
   NtTauTight           = new std::vector<Tau>;
   NtJetLoose           = new std::vector<Jet>;
   NtTruth              = new std::vector<Truth>;
   NtGenJet             = new std::vector<GenJet>; 
   NtTriggerObj         = new std::vector<TriggerObj>;
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
   NtTauTight->clear();
   NtJetLoose->clear();
   NtTruth->clear();
   NtGenJet->clear(); 
   NtTriggerObj->clear();
}

void Ntuple::fill()
{
   m_tree->Fill();
}
