#include "../include/NtupleProducer.h"

#include "TSystem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <sstream>

Tree *ntP;
TChain *ch;
Ntuple *nt;
Sync *sc;

unsigned int idx;

int main(int argc, char *argv[])
{
   if( argc < 7 )
     {
	std::cout << "NtupleProducer usage:" << std::endl;
	std::cout << "--file: input filename" << std::endl;
	std::cout << "--outfile: output filename" << std::endl;
	std::cout << "--tree: TTree name" << std::endl;
	std::cout << "--nmax: max number of events" << std::endl;
	std::cout << "--isdata: data or MC ?" << std::endl;
	std::cout << "--year: 2016, 2017 or 2018" << std::endl;
	std::cout << "--sync: produce sync ntuple (1-object,2-event)" << std::endl;
	exit(1);
     }

   const char *fname_str = "input.txt";
   const char *fname_out_str = "output.root";
   const char *stream_str = "FlatTree/tree";
   int nmax = -1;
   bool isdata = 0;
   int sync = 0;
   int year = 0;

   for(int i=0;i<argc;i++)
     {
        if( ! strcmp(argv[i],"--file") ) fname_str = argv[i+1];
        if( ! strcmp(argv[i],"--outfile") ) fname_out_str = argv[i+1];
        if( ! strcmp(argv[i],"--tree") ) stream_str = argv[i+1];
        if( ! strcmp(argv[i],"--nmax") ) nmax = atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--isdata") ) isdata = (bool) atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--year") ) year = atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--sync") ) sync = atoi(argv[i+1]);
     }

   std::string cmssw = "/storage_mnt/storage/user/kskovpen/analysis/Multilepton/CMSSW_10_2_13";

   const char *fname = fname_str;
   const char *stream = stream_str;
   const char *fname_out = fname_out_str;

   std::cout << "--file="   << fname      << std::endl;
   std::cout << "--outfile="<< fname_out  << std::endl;
   std::cout << "--tree="   << stream     << std::endl;
   std::cout << "--nmax="   << nmax       << std::endl;
   std::cout << "--isdata=" << isdata     << std::endl;
   std::cout << "--year="   << year     << std::endl;
   std::cout << "--sync="   << sync       << std::endl;

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

   //TString fname_sync_root = fname_out;
   //fname_sync_root += "_sync.root";
   TString fname_sync_root = "tHq_syncTree_Object.root";
   if(sync == 2) {fname_sync_root = "tHq_syncTree_Selection.root";}

   sc = new Sync(fname_sync_root.Data(),sync);
   sc->Init();
   sc->setBranchAddress();

   EventExt     ev;
   ElectronExt  el;
   MuonExt      mu;
   TauExt      tau;
   JetExt      jet;

   TruthExt  truth;
   GenJetExt genjet;
   TriggerObjExt trigObj;

//-- DEBUG : here, can read a list of event ids directly from a file, and then look specifically for them
//===============================
   std::vector<long int> evdebugid;
   std::vector<long int> evdebuglumi;
   /*
   ifstream file_in("/home-pbs/ntonon/tHq/IPHCNtuple_2017/CMSSW_9_4_3/src/IPHCNtuple/NtupleProducer/test/sync/input_3l_SR_METinTH_noTriggerMatch_tH_only_not_in_tHq.txt");  
   //if (file_in.is_open()) {std::cout <<"!"<<file_in.rdbuf();}
   
   string line;
   while(getline(file_in, line))
   {
   	//std::cout<<"line = "<<line<<std::endl;
   	std::stringstream ss(line);
	Int_t run, lumi, ev_id;
	char c1, c2;
	ss>>run>>c1>>lumi>>c2>>ev_id;
	evdebugid.push_back(ev_id);
	evdebugid.push_back(lumi);
	std::cout<<"DEBUG Ev ID = "<<ev_id<<std::endl;
   }*/
   //==========================
   
   evdebugid.push_back(2511395);
   
   int nlep = 0;
   int njet = 0;

   int n_presel_mu  = 0;
   int n_presel_el  = 0;
   int n_presel_tau = 0;
   int n_presel_jet = 0;
   int n_presel_jetFwd = 0;

   //-- JES //Ask for "Total" source of uncert.
   JetCorrectionUncertainty* jesTotal = 0;
   std::string jecFilesPath = cmssw+"/src/IPHCNtuple/NtupleProducer/data/jecFiles/";
   std::string jecMC = jecFilesPath+"Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt";
   std::string jecData = jecFilesPath+"Fall17_17Nov2017F_V32_DATA/Fall17_17Nov2017F_V32_DATA_UncertaintySources_AK4PFchs.txt";
   if(!isdata) {jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecMC.c_str(), "Total")));}
   else {jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters(jecData.c_str(), "Total")));}

   //-- JER
   //See : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#2017_data
   JME::JetResolution* JER_resolution = 0;
   JME::JetResolutionScaleFactor* JER_scalefactor = 0;
   std::string jerFilesPath = cmssw+"/src/IPHCNtuple/NtupleProducer/data/jerFiles/";
   std::string resolution_file = jerFilesPath+"Fall17_V3_MC_PtResolution_AK4PFchs.txt";
   std::string jer_sf_file = jerFilesPath+"Fall17_V3_MC_SF_AK4PFchs.txt";
   if(!isdata)
   {
      JER_resolution = new JME::JetResolution(resolution_file);
      JER_scalefactor = new JME::JetResolutionScaleFactor(jer_sf_file);
   }
   
   for(Long64_t i=0;i<nentries;i++)
     {
	if( i > nmax && nmax >= 0 ) break;

	if(i%10000 == 0) {cout<<"-- "<<i<<" / "<<nentries<<endl;}
	
	//cout<<"ientry "<<i<<endl;

	ch->GetEntry(i);
	
	nt->clearVar();
	sc->initVar();
	
	//Only keep events to debug -- disactivate if not debugging specific events //FIXME
	bool is_debug_event = false;
//-------------------
        ev.init();
	
	if( is_debug_event )
	  {
	     std::cout << "Event " << ev.id << std::endl;
	  }	
        ev.read(isdata,year,is_debug_event);

	if( is_debug_event )
	  {
	     std::cout << "Event " << ev.id << std::endl;
	  }	
	
	//cout<<"Event id = "<<ev.id<<endl;
	for(int k = 0; k < evdebugid.size(); k++)
	{
	   if(ev.id == evdebugid.at(k)) {is_debug_event = true;} //match ID
//		if(evdebuglumi.size() > 0 && ev.lumi != evdebuglumi.at(k)) {is_debug_event = false;} //also match lumi
	}
/*	if(!is_debug_event && evdebugid.size() > 0) {continue;} //only debug
	else if(evdebugid.size() > 0) {cout<<endl<<endl<<endl<<endl<<setprecision(12)<<"== Run/Lumi/Event "<<ev.run<<" / "<<ev.lumi<<" / "<<ev.id<<" =="<<endl;}
	*/
//-------------------
	
        int n_mu_evt = 0, n_mu_fakeable_evt = 0;
	
        // muons
        for(int j=0;j<ntP->mu_n;j++)
	  {
	     idx = j;

	     mu.init();
	     mu.read(isdata);
	     mu.sel(is_debug_event,year);

	     if( mu.isLooseTTH )
	       {
		  nt->NtMuonLooseExt->push_back(mu);
		  n_mu_evt++;
	       }
	     if( mu.isFakeableTTH )
	       {
		  nt->NtMuonFakeableExt->push_back(mu);
		  n_mu_fakeable_evt++;
	       }
	     if( mu.isTightTTH )
	       {
		  nt->NtMuonTightExt->push_back(mu);
	       }
	  }

        int n_el_evt = 0, n_el_fakeable_evt = 0;

	// electrons
        bool pass_eventVeto = true;
        for(int j=0;j<ntP->el_n;j++)
	  {
	     idx = j;

	     el.init();
	     el.read(isdata);
	     el.sel(is_debug_event,year);

	     if( el.isLooseTTH  )
	       {
		  nt->NtElectronLooseExt->push_back(el);
		  n_el_evt++;
	       }
	     if( el.isFakeableTTH  )
	       {
		  nt->NtElectronFakeableExt->push_back(el);
		  n_el_fakeable_evt++;
	       }
	     if( el.isTightTTH  )
	       {
		  nt->NtElectronTightExt->push_back(el);
	       }
	  }

	int n_tau_evt = 0;

        // taus
        for(int j=0;j<ntP->tau_n;j++)
	  {
	     idx = j;

	     tau.init();
	     tau.read(isdata);
	     tau.sel(is_debug_event,year);

	      if( tau.isLooseTTH )
	      {
		  nt->NtTauLooseExt->push_back(tau);

	      }
	     if( tau.isFakeableTTH )
	       {
		  nt->NtTauFakeableExt->push_back(tau);
  		  n_tau_evt++; //Sync : 'presel' taus = FO...?
	       }
	     if( tau.isMediumTTH )
	       {
		  nt->NtTauMediumExt->push_back(tau);
	       }
	  }

        int n_jet_evt = 0;
	int n_jetFwd_evt = 0;
	int nBL = 0;
	
        // jets
        for(int j=0;j<ntP->jet_n;j++)
	{
	     idx = j;

	     jet.init();
	     jet.read(isdata);

         //JER (4-momentum smearing + uncertainties)
         if(!isdata)
         {
            JME::JetParameters param_JER;
            param_JER.setJetPt(jet.pt);
            param_JER.setJetEta(jet.eta);
            param_JER.setRho(ev.rho);

            float JER_res = JER_resolution->getResolution(param_JER);
            float JER_sf = JER_scalefactor->getScaleFactor(param_JER);
            float JER_sf_up = JER_scalefactor->getScaleFactor(param_JER, Variation::UP);
            float JER_sf_down = JER_scalefactor->getScaleFactor(param_JER, Variation::DOWN);

            if(sync == 0) //Don't apply JER for sync
            {
               jet.apply_JER_smearing(isdata, JER_res, JER_sf, JER_sf_up, JER_sf_down); //Fill JER variables and apply JER SF
            }
         }

         //JES uncertainties
	 //NB : also stored in FlatTrees (ntP->jet_Unc) ; re-computed here, since should apply JER first
	 jesTotal->setJetPt(ntP->jet_pt->at(idx));
	 jesTotal->setJetEta(ntP->jet_eta->at(idx));

         if(sync == 0) //Don't apply JES for sync
         {
	    float jes_unc = jesTotal->getUncertainty(true);
            jet.setJESUncertainty(isdata, jes_unc); //Fill pt_jes_up & pt_jes_down
         }

	 jet.sel(sync, is_debug_event, year);

	 if( jet.isLooseTTH )
	 {
	 	nt->NtJetLooseExt->push_back(jet);
		n_jet_evt++;

		if(jet.isLooseBTag) nBL++;
	 }
//	   else if(jet.isSoftLooseTTH) nt->NtJetLooseSoftExt->push_back(jet);
	   else if ( jet.isLooseFwdTTH )
	 {
	        nt->NtJetLooseFwdExt->push_back(jet);
		n_jetFwd_evt++;
	 }	   

	}


        if( !isdata )
	  {
	     // genjets
	     for(int j=0;j<ntP->genJet_n;j++)
	       {
		  idx = j;

		  genjet.init();
		  genjet.read();

		  if (genjet.sel()) nt->NtGenJetExt->push_back(genjet);
	       }

	     // truth
	     truth.init();
	     truth.read();
	     truth.readMultiLepton();

	     nt->NtTruthExt->push_back(truth);
	  }

	if( n_el_evt > 0 ) n_presel_el++;
	if( n_mu_evt > 0 ) n_presel_mu++;
	if( n_tau_evt > 0 ) n_presel_tau++;
	if( n_jet_evt > 0 ) n_presel_jet++;
	if( n_jetFwd_evt > 0 ) n_presel_jetFwd++;
	
	sc->get(nt,n_el_evt,n_mu_evt,n_tau_evt,n_jet_evt,n_jetFwd_evt,n_el_fakeable_evt,n_mu_fakeable_evt, nBL);

	bool pass = sc->fill(nt,&ev, is_debug_event);

	nt->NtEventExt->push_back(ev);

	nt->fill();

	if( pass ) nt->write();
     }

   if( sync )
     {
	std::cout << "n_presel_mu = " << n_presel_mu  << std::endl;
	std::cout << "n_presel_el = " << n_presel_el  << std::endl;
	std::cout << "n_presel_tau = " << n_presel_tau << std::endl;
	std::cout << "n_presel_jet = " << n_presel_jet << std::endl;
	std::cout << "n_presel_jetFwd = " << n_presel_jetFwd << std::endl;
     }

   delete nt;
   delete sc;

   delete jesTotal;
   
   // if(!isdata)
   // {
   //    delete JER_resolution;
   //    delete JER_scalefactor;
   // }

  return 0;
}
