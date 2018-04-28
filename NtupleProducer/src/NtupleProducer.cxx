#include "../include/NtupleProducer.h"

#include "Riostream.h"
#include "TSystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>

Tree *ntP;
TChain *ch;
Ntuple *nt;
Sync *sc;
std::vector<int> *evdebug;

unsigned int idx;

int main(int argc, char *argv[])
{
   if( argc < 6 )
     {
	std::cout << "NtupleProducer usage:" << std::endl;
	std::cout << "--file: input filename" << std::endl;
	std::cout << "--outfile: output filename" << std::endl;
	std::cout << "--tree: TTree name" << std::endl;
	std::cout << "--nmax: max number of events" << std::endl;
	std::cout << "--isdata: data or MC ?" << std::endl;
	std::cout << "--dosync: produce sync ntuple" << std::endl;
	exit(1);
     }
   
   const char *fname_str = "input.txt";
   const char *fname_out_str = "output.root";
   const char *stream_str = "FlatTree/tree";
   int nmax = -1;
   bool isdata = 0;
   bool dosync = 0;

   for(int i=0;i<argc;i++)
     {
        if( ! strcmp(argv[i],"--file") ) fname_str = argv[i+1];
        if( ! strcmp(argv[i],"--outfile") ) fname_out_str = argv[i+1];
        if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
        if( ! strcmp(argv[i],"--nmax") ) nmax = atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--isdata") ) isdata = (bool) atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--dosync") ) dosync = (bool) atoi(argv[i+1]);
     }
   
   std::string cmssw = std::string(getenv("CMSSW_BASE"));
   
   const char *fname = fname_str;
   const char *stream = stream_str;
   const char *fname_out = fname_out_str;
   
   std::cout << "--file="   << fname      << std::endl;
   std::cout << "--outfile="<< fname_out  << std::endl;
   std::cout << "--tree="   << stream     << std::endl;
   std::cout << "--nmax="   << nmax       << std::endl;
   std::cout << "--isdata=" << isdata     << std::endl;
   std::cout << "--dosync=" << dosync     << std::endl;
   
   Tree tree(0,const_cast<char*>(fname),stream);
   ntP = &tree;
   
   ch = tree.fChain;
   Long64_t nentries = ch->GetEntries();
   ntP->registerInputBranches(ch);

   TString fname_out_root = fname_out;
   fname_out_root += ".root";
   nt = new Ntuple(fname_out_root.Data());

   nt->Init();
   nt->createVar();
   nt->setBranchAddress();

   if( dosync )
     {	
	std::cout << "Running in sync mode" << std::endl;
	TString fname_sync_root = fname_out;
	fname_sync_root += "_sync.root";
	sc = new Sync(fname_sync_root.Data());
	sc->Init();
	sc->setBranchAddress();
     }
   
   Event     ev;
   Electron  el;
   Muon      mu;
   Tau      tau;
   Jet      jet;
   
   Truth  truth;
   GenJet genjet;
   TriggerObj trigObj;
   
   evdebug = new std::vector<int>();
   //   evdebug->push_back(120);

   int nlep = 0;
   int njet = 0;
   
   int n_presel_mu  = 0;
   int n_presel_el  = 0;
   int n_presel_tau = 0;
   int n_presel_jet = 0;
   
   JetCorrectionUncertainty *jesTotal;

   std::string jecFilesPath = cmssw+"/src/IPHCNtuple/NtupleProducer/data/jecFiles/";
   std::string jecMC = jecFilesPath+"Fall17_17Nov2017_V6_MC_UncertaintySources_AK4PFchs.txt";
   std::string jecData = jecFilesPath+"Fall17_17Nov2017F_V6_DATA_UncertaintySources_AK4PFchs.txt";
   
   if( isdata == false ) jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecMC.c_str(), "Total")));
   else	jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecData.c_str(), "Total")));   
   
   for(Long64_t i=0;i<nentries;i++)
    {
       if( i > nmax && nmax >= 0 ) break;

       ch->GetEntry(i);
       nt->clearVar();
       if( dosync ) sc->initVar();

        //	if( !(isHtoWW || isHtoZZ || isHtoTT) ) continue;

        // event
        ev.init();
        ev.read(isdata);

        //	if( isHtoWW ) ev._tth_channel = 0;
        //	else if( isHtoZZ ) ev._tth_channel = 1;
        //	else if( isHtoTT ) ev._tth_channel = 2;

        nt->NtEvent->push_back(ev);

        int n_mu_evt = 0;

        // muons
        for(int j=0;j<ntP->mu_n;j++)
        {
            idx = j;

            mu.init();
            mu.read();

            if( mu.sel())
            {
                nt->NtMuon->push_back(mu);
                n_mu_evt++;
            }
        }

        int n_el_evt = 0;

        // electrons
        for(int j=0;j<ntP->el_n;j++)
        {
            idx = j;

            el.init();
            el.read();

            if( el.sel()  )
            {
                nt->NtElectron->push_back(el);
                n_el_evt++;
            }
        }

       int n_tau_evt = 0;

        // taus
        for(int j=0;j<ntP->tau_n;j++)
        {
            idx = j;

            tau.init();
            tau.read();

            if (tau.sel())
            {
                nt->NtTau->push_back(tau);
                n_tau_evt++;
            }
        }

        int n_jet_evt = 0;

        // jets
        for(int j=0;j<ntP->jet_n;j++)
	 {
            idx = j;
	    
            jet.init();
            jet.read(isdata);
	    
            jesTotal->setJetPt(ntP->jet_pt->at(idx));
            jesTotal->setJetEta(ntP->jet_eta->at(idx));
	    
	    jet.setJESUncertainty(jesTotal->getUncertainty(true));
	    
	    if( jet.sel() )
	      {
		 nt->NtJet->push_back(jet);
		 n_jet_evt++;
	      }
        }

        //trigger objects
        /*for(int j=0;j<ntP->triggerobject_n;j++)
          {
          idx = j;

          trigObj.init();
          trigObj.read();

          if (trigObj.sel()) nt->NtTriggerObj->push_back(trigObj);
          }*/

        if (!isdata)
        {
            // genjets
            for(int j=0;j<ntP->genJet_n;j++)
            {
                idx = j;

                genjet.init();
                genjet.read();

                if (genjet.sel()) nt->NtGenJet->push_back(genjet);
            }

            // truth
            truth.init();
            truth.read();
            truth.readMultiLepton();

            nt->NtTruth->push_back(truth);
        }
       
       if( n_el_evt > 0 ) n_presel_el++;
       if( n_mu_evt > 0 ) n_presel_mu++;
       if( n_tau_evt > 0 ) n_presel_tau++;
       if( n_jet_evt > 0 ) n_presel_jet++;

       // Sync ntuple variables
       sc->nEvent = i;
       sc->ls = ntP->ev_lumi;
       sc->run = ntP->ev_run;
       sc->n_presel_mu = n_presel_mu;
       sc->n_fakeablesel_mu = -9999;
       sc->n_mvasel_mu = -9999;
       sc->n_presel_ele = n_presel_el;
       sc->n_fakeablesel_ele = -9999;
       sc->n_mvasel_ele = -9999;
       sc->n_presel_tau = n_presel_tau;
       sc->n_presel_jet = n_presel_jet;
       
       nt->fill();
       if( dosync ) sc->fill();
    }

   if( dosync )
     {	
	std::cout << "n_presel_mu = " << n_presel_mu  << std::endl;
	std::cout << "n_presel_el = " << n_presel_el  << std::endl;
	std::cout << "n_presel_tau = " << n_presel_tau << std::endl;
	std::cout << "n_presel_jet = " << n_presel_jet << std::endl;
     }   
   
   delete evdebug;
   delete nt;
   delete sc;
}
