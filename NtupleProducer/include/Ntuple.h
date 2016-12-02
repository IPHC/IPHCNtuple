#ifndef NTUPLE_H
#define NTUPLE_H

#include "Event.h"
#include "Electron.h"
#include "Muon.h"
#include "Tau.h"
#include "Jet.h"
#include "Truth.h"
#include "GenJet.h"
#include "TriggerObj.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

class Ntuple
{
    public:

        Ntuple(std::string fname_out);
        virtual ~Ntuple();

        std::vector<Event>*      NtEvent;
        std::vector<Electron>*   NtElectron;
        std::vector<Muon>*       NtMuon;
        std::vector<Tau>*        NtTau;
        std::vector<Jet>*        NtJet;
        std::vector<Truth>*      NtTruth;
        std::vector<GenJet>*     NtGenJet;
        std::vector<TriggerObj>* NtTriggerObj;

        void Init();
 
        void setBranchAddress();
        void createVar();
        void clearVar();
        void fill();

        TFile*  m_file;

    private:

        TTree*  m_tree;
        TChain* m_chain;
		std::string _fname_out;
};

#endif
