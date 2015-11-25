#include "include/NtupleProducer.h"

ClassImp(Truth)

Truth::Truth()
{
}

Truth::~Truth()
{
}

void Truth::read()
{
    int UNINT = -666;

    int gen_n = ntP->gen_n;
    
    int gen_n_sel = 0;
    for(int i=0;i<gen_n;i++)
    {
        int status = ntP->gen_status->at(i);
        if( status != 1 && status != 3 ) continue;
        _gen_pt.push_back(ntP->gen_pt->at(i));
        _gen_eta.push_back(ntP->gen_eta->at(i));
        _gen_phi.push_back(ntP->gen_phi->at(i));
        _gen_m.push_back(ntP->gen_m->at(i));
        _gen_id.push_back(ntP->gen_id->at(i));
        _gen_status.push_back(status);

        int mom = ntP->gen_mother_index->at(i);

        if( mom >= ntP->gen_index->size() )
        {

            std::cout << "Problem with gen length" << std::endl;
            exit(1);
        }

        int mother_id = ntP->gen_id->at(mom);

        _gen_mother_id.push_back(mother_id);

        gen_n_sel++;
    }
    _gen_PVz = ntP->gen_PVz;
   
    _mc_truth_n = 0;

    int mc_truth_h0W1_id = ntP->mc_truth_h0W1_id;
    if( mc_truth_h0W1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0W1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0W1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0W1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0W1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0W1_E);
        _mc_truth_label.push_back(12);
        _mc_truth_n++;
    }
    int mc_truth_h0Wl1_id = ntP->mc_truth_h0Wl1_id;
    if( mc_truth_h0Wl1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wl1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wl1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wl1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wl1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wl1_E);
        _mc_truth_label.push_back(120);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq11_id = ntP->mc_truth_h0Wq11_id;
    if( mc_truth_h0Wq11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq11_E);
        _mc_truth_label.push_back(122);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq21_id = ntP->mc_truth_h0Wq21_id;
    if( mc_truth_h0Wq21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq21_E);
        _mc_truth_label.push_back(123);
        _mc_truth_n++;
    }
    int mc_truth_h0W2_id = ntP->mc_truth_h0W2_id;
    if( mc_truth_h0W2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0W2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0W2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0W2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0W2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0W2_E);
        _mc_truth_label.push_back(13);
        _mc_truth_n++;
    }
    int mc_truth_h0Wl2_id = ntP->mc_truth_h0Wl2_id;
    if( mc_truth_h0Wl2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wl2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wl2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wl2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wl2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wl2_E);
        _mc_truth_label.push_back(130);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq12_id = ntP->mc_truth_h0Wq12_id;
    if( mc_truth_h0Wq12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq12_E);
        _mc_truth_label.push_back(132);
        _mc_truth_n++;
    }
    int mc_truth_h0Wq22_id = ntP->mc_truth_h0Wq22_id;
    if( mc_truth_h0Wq22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Wq22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Wq22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Wq22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Wq22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Wq22_E);
        _mc_truth_label.push_back(133);
        _mc_truth_n++;
    }
    int mc_truth_h0Z1_id = ntP->mc_truth_h0Z1_id;
    if( mc_truth_h0Z1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Z1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Z1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Z1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Z1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Z1_E);
        _mc_truth_label.push_back(14);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl11_id = ntP->mc_truth_h0Zl11_id;
    if( mc_truth_h0Zl11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl11_E);
        _mc_truth_label.push_back(140);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl21_id = ntP->mc_truth_h0Zl21_id;
    if( mc_truth_h0Zl21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl21_E);
        _mc_truth_label.push_back(141);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq11_id = ntP->mc_truth_h0Zq11_id;
    if( mc_truth_h0Zq11_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq11_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq11_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq11_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq11_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq11_E);
        _mc_truth_label.push_back(142);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq21_id = ntP->mc_truth_h0Zq21_id;
    if( mc_truth_h0Zq21_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq21_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq21_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq21_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq21_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq21_E);
        _mc_truth_label.push_back(143);
        _mc_truth_n++;
    }
    int mc_truth_h0Z2_id = ntP->mc_truth_h0Z2_id;
    if( mc_truth_h0Z2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Z2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Z2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Z2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Z2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Z2_E);
        _mc_truth_label.push_back(15);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl12_id = ntP->mc_truth_h0Zl12_id;
    if( mc_truth_h0Zl12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl12_E);
        _mc_truth_label.push_back(150);
        _mc_truth_n++;
    }
    int mc_truth_h0Zl22_id = ntP->mc_truth_h0Zl22_id;
    if( mc_truth_h0Zl22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zl22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zl22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zl22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zl22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zl22_E);
        _mc_truth_label.push_back(151);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq12_id = ntP->mc_truth_h0Zq12_id;
    if( mc_truth_h0Zq12_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq12_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq12_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq12_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq12_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq12_E);
        _mc_truth_label.push_back(152);
        _mc_truth_n++;
    }
    int mc_truth_h0Zq22_id = ntP->mc_truth_h0Zq22_id;
    if( mc_truth_h0Zq22_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0Zq22_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0Zq22_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0Zq22_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0Zq22_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0Zq22_E);
        _mc_truth_label.push_back(153);
        _mc_truth_n++;
    }
    int mc_truth_h0tau1_id = ntP->mc_truth_h0tau1_id;
    if( mc_truth_h0tau1_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0tau1_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0tau1_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0tau1_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0tau1_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0tau1_E);
        _mc_truth_label.push_back(16);
        _mc_truth_n++;
    }
    int mc_truth_h0tau2_id = ntP->mc_truth_h0tau2_id;
    if( mc_truth_h0tau2_id != UNINT )
    {
        _mc_truth_id.push_back(mc_truth_h0tau2_id);
        _mc_truth_pt.push_back(ntP->mc_truth_h0tau2_pt);
        _mc_truth_eta.push_back(ntP->mc_truth_h0tau2_eta);
        _mc_truth_phi.push_back(ntP->mc_truth_h0tau2_phi);
        _mc_truth_E.push_back(ntP->mc_truth_h0tau2_E);
        _mc_truth_label.push_back(17);
        _mc_truth_n++;	
    }
   
    _metGen_px = ntP->metGen_px;
    _metGen_py = ntP->metGen_py;
    _metGen_pt = ntP->metGen_pt;
    _metGen_phi = ntP->metGen_phi;
    _metGen_sumet = ntP->metGen_sumet;
    _metGen_MuonEt = ntP->metGen_MuonEt;
	 
}

void Truth::readMultiLepton()
{  
  int UNINT = -666;

  if (ntP->mc_truth_t1_id==-666 || ntP->mc_truth_t2_id==-666)  return;
  
  if (ntP->mc_truth_h0_id>-666 && ntP->mc_truth_h0Wl1_id>-666 && ntP->mc_truth_h0Wl2_id>-666)
  {
  
      _Leptons_id.push_back(ntP->mc_truth_h0Wl1_id);
      _Leptons_id.push_back(ntP->mc_truth_h0Wl2_id);
      _Leptons_pt.push_back(ntP->mc_truth_h0Wl1_pt);
      _Leptons_pt.push_back(ntP->mc_truth_h0Wl2_pt);
      _Leptons_eta.push_back(ntP->mc_truth_h0Wl1_eta);
      _Leptons_eta.push_back(ntP->mc_truth_h0Wl2_eta);
      _Leptons_phi.push_back(ntP->mc_truth_h0Wl1_phi);
      _Leptons_phi.push_back(ntP->mc_truth_h0Wl2_phi);
      _Leptons_E.push_back(ntP->mc_truth_h0Wl1_E);
      _Leptons_E.push_back(ntP->mc_truth_h0Wl2_E);
      
      _boson_decay = 0;
      
  }
  
  
  else if (ntP->mc_truth_h0_id>-666 && ntP->mc_truth_h0Wq11_id>-666 && ntP->mc_truth_h0Wl2_id>-666)
  {
        
    _Leptons_id.push_back(ntP->mc_truth_h0Wl2_id);
    _Leptons_pt.push_back(ntP->mc_truth_h0Wl2_pt);
    _Leptons_eta.push_back(ntP->mc_truth_h0Wl2_eta);
    _Leptons_phi.push_back(ntP->mc_truth_h0Wl2_phi);
    _Leptons_E.push_back(ntP->mc_truth_h0Wl2_E);
   
    _QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq11_id);
    _QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq11_pt);
    _QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq11_eta);
    _QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq11_phi);
    _QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq11_E);
    
    _QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
    _QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq12_pt);
    _QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq12_eta);
    _QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq12_phi);
    _QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq12_E);
   
    _boson_decay = 1;   
    
    float dr1 = 999; 
    float dr2 = 999;
    float drMax1 = 0.4;
    float drMax2 = 0.4;
    int iMax1 = -1;
    int iMax2 = -1;
    
    for (int i=0; i<nt->NtGenJet->size(); i++)
    {
      dr1 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_h0Wq11_eta, ntP->mc_truth_h0Wq11_phi);
      dr2 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_h0Wq21_eta, ntP->mc_truth_h0Wq21_phi);
      if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
      if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
    }
    
    if (iMax1 == -1)
    {
      _JetsFromWs_id.push_back(UNINT);
      _JetsFromWs_pt.push_back(UNINT);
      _JetsFromWs_eta.push_back(UNINT);
      _JetsFromWs_phi.push_back(UNINT);
      _JetsFromWs_E.push_back(UNINT);
    }
    else 
    {
      _JetsFromWs_id.push_back(ntP->mc_truth_h0Wq11_id);
      _JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt());
      _JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta());
      _JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi());
      _JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E());
    }  
    if (iMax2 == -1)
    {
      _JetsFromWs_id.push_back(UNINT);
      _JetsFromWs_pt.push_back(UNINT);
      _JetsFromWs_eta.push_back(UNINT);
      _JetsFromWs_phi.push_back(UNINT);
      _JetsFromWs_E.push_back(UNINT);
    }
    else 
    {
      _JetsFromWs_id.push_back(ntP->mc_truth_h0Wq21_id);
      _JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt());
      _JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta());
      _JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi());
      _JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E());
    }
      
  }
  
  
  else if (ntP->mc_truth_h0_id>-666 && ntP->mc_truth_h0Wq22_id>-666 && ntP->mc_truth_h0Wl1_id>-666)
  {
                 
      _Leptons_id.push_back(ntP->mc_truth_h0Wl1_id);
      _Leptons_pt.push_back(ntP->mc_truth_h0Wl1_pt);
      _Leptons_eta.push_back(ntP->mc_truth_h0Wl1_eta);
      _Leptons_phi.push_back(ntP->mc_truth_h0Wl1_phi);
      _Leptons_E.push_back(ntP->mc_truth_h0Wl1_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq12_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq12_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq12_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq12_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq22_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq22_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq22_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq22_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq22_E);
      
      _boson_decay = 1;
     
      float dr1 = 999; 
      float dr2 = 999;
      float drMax1 = 0.4;
      float drMax2 = 0.4;
      int iMax1 = -1;
      int iMax2 = -1;
    
      for (int i=0; i<nt->NtGenJet->size(); i++)
      {
     	dr1 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_h0Wq12_eta, ntP->mc_truth_h0Wq12_phi);
     	dr2 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_h0Wq22_eta, ntP->mc_truth_h0Wq22_phi);
     	if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
     	if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
      }
      
      if (iMax1 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
     	_JetsFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
     	_JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt());
     	_JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta());
     	_JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi());
     	_JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E());
      }
      if (iMax2 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
        _JetsFromWs_id.push_back(ntP->mc_truth_h0Wq22_id);
        _JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt());
        _JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta());
        _JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi());
        _JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E());
      }
            
    }
    
    
    else if (ntP->mc_truth_W_id>-666 && ntP->mc_truth_Wl_id>-666)
    {
       
      _Leptons_id.push_back(ntP->mc_truth_Wl_id);
      _Leptons_pt.push_back(ntP->mc_truth_Wl_pt);
      _Leptons_eta.push_back(ntP->mc_truth_Wl_eta);
      _Leptons_phi.push_back(ntP->mc_truth_Wl_phi);
      _Leptons_E.push_back(ntP->mc_truth_Wl_E);
      
      _boson_decay = 3;    
     
    }


    else if (ntP->mc_truth_Z_id>-666 && ntP->mc_truth_Zl1_id>-666 && ntP->mc_truth_Zl2_id>-666)
    {
      
      _Leptons_id.push_back(ntP->mc_truth_Zl1_id);      
      _Leptons_pt.push_back(ntP->mc_truth_Zl1_pt);
      _Leptons_eta.push_back(ntP->mc_truth_Zl1_eta);
      _Leptons_phi.push_back(ntP->mc_truth_Zl1_phi);
      _Leptons_E.push_back(ntP->mc_truth_Zl1_E);
      
      _Leptons_id.push_back(ntP->mc_truth_Zl2_id);      
      _Leptons_pt.push_back(ntP->mc_truth_Zl2_pt);
      _Leptons_eta.push_back(ntP->mc_truth_Zl2_eta);
      _Leptons_phi.push_back(ntP->mc_truth_Zl2_phi);
      _Leptons_E.push_back(ntP->mc_truth_Zl2_E);
      
      _boson_decay = 2;      
      
    }
    
    else if (ntP->mc_truth_gammal1_id>-666 && ntP->mc_truth_gammal2_id>-666)
    {
      
      _Leptons_id.push_back(ntP->mc_truth_gammal1_id);      
      _Leptons_pt.push_back(ntP->mc_truth_gammal1_pt);
      _Leptons_eta.push_back(ntP->mc_truth_gammal1_eta);
      _Leptons_phi.push_back(ntP->mc_truth_gammal1_phi);
      _Leptons_E.push_back(ntP->mc_truth_gammal1_E);
      
      _Leptons_id.push_back(ntP->mc_truth_gammal2_id);      
      _Leptons_pt.push_back(ntP->mc_truth_gammal2_pt);
      _Leptons_eta.push_back(ntP->mc_truth_gammal2_eta);
      _Leptons_phi.push_back(ntP->mc_truth_gammal2_phi);
      _Leptons_E.push_back(ntP->mc_truth_gammal2_E);
    
      _boson_decay = 2;    
      
    }
    
    else return;


    if (ntP->mc_truth_tWq11_id>-666 && ntP->mc_truth_tWl2_id>-666)
    {
     
      _Leptons_id.push_back(ntP->mc_truth_tWl2_id);      
      _Leptons_pt.push_back(ntP->mc_truth_tWl2_pt);
      _Leptons_eta.push_back(ntP->mc_truth_tWl2_eta);
      _Leptons_phi.push_back(ntP->mc_truth_tWl2_phi);
      _Leptons_E.push_back(ntP->mc_truth_tWl2_E);
      
      _Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb1_E);
      
      _Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb2_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_tWq11_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_tWq11_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_tWq11_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_tWq11_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_tWq11_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_tWq21_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_tWq21_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_tWq21_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_tWq21_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_tWq21_E);
        
      _ttbar_decay = 1;
            
      float dr1 = 999; 
      float dr2 = 999;
      float drMax1 = 0.4;
      float drMax2 = 0.4;
      int iMax1 = -1;
      int iMax2 = -1;
    
      for (int i=0; i<nt->NtGenJet->size(); i++)
      {
     	dr1 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_tWq11_eta, ntP->mc_truth_tWq11_phi);
     	dr2 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_tWq22_eta, ntP->mc_truth_tWq22_phi);
     	if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
     	if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
      }
      
      if (iMax1 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
     	_JetsFromWs_id.push_back(ntP->mc_truth_tWq11_id);
     	_JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt());
     	_JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta());
     	_JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi());
     	_JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E());
      }
      if (iMax2 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
        _JetsFromWs_id.push_back(ntP->mc_truth_tWq22_id);
        _JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt());
        _JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta());
        _JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi());
        _JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E());
      }
       
    }
    
    
    else if (ntP->mc_truth_tWq22_id>-666 && ntP->mc_truth_tWl1_id>-666)
    {      
      
      _Leptons_id.push_back(ntP->mc_truth_tWl1_id);      
      _Leptons_pt.push_back(ntP->mc_truth_tWl1_pt);
      _Leptons_eta.push_back(ntP->mc_truth_tWl1_eta);
      _Leptons_phi.push_back(ntP->mc_truth_tWl1_phi);
      _Leptons_E.push_back(ntP->mc_truth_tWl1_E);
      
      _Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb1_E);
      
      _Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb2_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_tWq12_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_tWq12_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_tWq12_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_tWq12_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_tWq12_E);
      
      _QuarksFromWs_id.push_back(ntP->mc_truth_tWq22_id);
      _QuarksFromWs_pt.push_back(ntP->mc_truth_tWq22_pt);
      _QuarksFromWs_eta.push_back(ntP->mc_truth_tWq22_eta);
      _QuarksFromWs_phi.push_back(ntP->mc_truth_tWq22_phi);
      _QuarksFromWs_E.push_back(ntP->mc_truth_tWq22_E);
     
      float dr1 = 999; 
      float dr2 = 999;
      float drMax1 = 0.4;
      float drMax2 = 0.4;
      int iMax1 = -1;
      int iMax2 = -1;
    
      for (int i=0; i<nt->NtGenJet->size(); i++)
      {
     	dr1 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_tWq12_eta, ntP->mc_truth_tWq12_phi);
     	dr2 = GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(), ntP->mc_truth_tWq22_eta, ntP->mc_truth_tWq22_phi);
     	if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
     	if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
      }
      
      if (iMax1 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
     	_JetsFromWs_id.push_back(ntP->mc_truth_tWq12_id);
     	_JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt());
     	_JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta());
     	_JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi());
     	_JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E());
      }
      if (iMax2 == -1)
      {
     	_JetsFromWs_id.push_back(UNINT);
     	_JetsFromWs_pt.push_back(UNINT);
     	_JetsFromWs_eta.push_back(UNINT);
     	_JetsFromWs_phi.push_back(UNINT);
     	_JetsFromWs_E.push_back(UNINT);
      }
      else 
      {
        _JetsFromWs_id.push_back(ntP->mc_truth_tWq22_id);
        _JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt());
        _JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta());
        _JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi());
        _JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E());
      }
      
      _ttbar_decay = 1;      
      
   }
   
   else if (ntP->mc_truth_tWl1_id>-666 && ntP->mc_truth_tWl2_id>-666 && (_boson_decay==1 || _boson_decay==3))
   {
   
      _Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb1_E);
      
      _Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
      _Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
      _Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
      _Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
      _Bjets_E.push_back(ntP->mc_truth_tb2_E);
      
      _Leptons_id.push_back(ntP->mc_truth_tWl1_id);      
      _Leptons_pt.push_back(ntP->mc_truth_tWl1_pt);
      _Leptons_eta.push_back(ntP->mc_truth_tWl1_eta);
      _Leptons_phi.push_back(ntP->mc_truth_tWl1_phi);
      _Leptons_E.push_back(ntP->mc_truth_tWl1_E);
      
      _Leptons_id.push_back(ntP->mc_truth_tWl2_id);      
      _Leptons_pt.push_back(ntP->mc_truth_tWl2_pt);
      _Leptons_eta.push_back(ntP->mc_truth_tWl2_eta);
      _Leptons_phi.push_back(ntP->mc_truth_tWl2_phi);
      _Leptons_E.push_back(ntP->mc_truth_tWl2_E);     
       
      _ttbar_decay = 2;
     
   }
   
   else return;

   
   for (int i=0; i<nt->NtGenJet->size(); i++)
   { 
     if ( fabs(nt->NtGenJet->at(i).eta())>2.5 ) continue;
     if ( nt->NtGenJet->at(i).pt()<25 ) continue;
     
     if ( GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(),_Bjets_eta[0],_Bjets_phi[0]) < 0.4 ||
          GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(),_Bjets_eta[1],_Bjets_phi[1]) < 0.4 || 
	  GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(),_Leptons_eta[0],_Leptons_phi[0]) < 0.4 ||
          GetDeltaR(nt->NtGenJet->at(i).eta(),nt->NtGenJet->at(i).phi(),_Leptons_eta[1],_Leptons_phi[1]) < 0.4 ) continue;
	   
     _AllJets_id.push_back(0);
     _AllJets_pt.push_back(nt->NtGenJet->at(i).pt());
     _AllJets_eta.push_back(nt->NtGenJet->at(i).eta());
     _AllJets_phi.push_back(nt->NtGenJet->at(i).phi());
     _AllJets_E.push_back(nt->NtGenJet->at(i).E());
    
   }
      
   if (_AllJets_E.size()==1) 
   {
   
    _Jets_id.push_back(0);
    _Jets_pt.push_back(_AllJets_pt[0]);
    _Jets_eta.push_back(_AllJets_eta[0]);
    _Jets_phi.push_back(_AllJets_phi[0]);
    _Jets_E.push_back(_AllJets_E[0]);
   
   }
   else if (_AllJets_E.size()>=2) 
   {
      //mc_3l_category = 2;

      float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
      float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
      
      for (unsigned int ij=0; ij<_AllJets_E.size(); ij++){
        if (_AllJets_pt[ij] > pt_max ) {
           pt_max2 = pt_max;
           ij2 = ij1;
           pt_max = _AllJets_pt[ij];
           ij1 = ij;
         }
         if ( _AllJets_pt[ij] < pt_max && _AllJets_pt[ij] > pt_max2){
           pt_max2 = _AllJets_pt[ij];
           ij2 = ij;
         }
         for (unsigned int ik=0; ik<_AllJets_E.size(); ik++){
           if (ik==ij) continue;
	   
	   TLorentzVector Jij;
	   Jij.SetPtEtaPhiE(_AllJets_pt[ij], _AllJets_eta[ij], _AllJets_phi[ij], _AllJets_E[ij]);
	   TLorentzVector Jik;
	   Jik.SetPtEtaPhiE(_AllJets_pt[ik], _AllJets_eta[ik], _AllJets_phi[ik], _AllJets_E[ik]);
	   
           if (fabs((Jij+Jik).M()-80.419)<diffmass_min){
             ik1=ij;
             ik2=ik;
             diffmass_min = fabs((Jij+Jik).M()-80.419);
           }
	   if ((Jij+Jik).M()<mass_min){
	     il1=ij;
             il2=ik;
	     mass_min = (Jij+Jik).M();
	   }
         }
      }
      if (ij1!=-1 && ij2!=-1) {
      
	_JetsHighestPt_id.push_back(1);
        _JetsHighestPt_pt.push_back(_AllJets_pt[ij1]);
        _JetsHighestPt_eta.push_back(_AllJets_eta[ij1]);
        _JetsHighestPt_phi.push_back(_AllJets_phi[ij1]);
        _JetsHighestPt_E.push_back(_AllJets_E[ij1]);
	
	_JetsHighestPt_id.push_back(1);
        _JetsHighestPt_pt.push_back(_AllJets_pt[ij2]);
        _JetsHighestPt_eta.push_back(_AllJets_eta[ij2]);
        _JetsHighestPt_phi.push_back(_AllJets_phi[ij2]);
        _JetsHighestPt_E.push_back(_AllJets_E[ij2]);
		
      }
      if (ik1!=-1 && ik2!=-1){
      
	_JetsClosestMw_id.push_back(2);
        _JetsClosestMw_pt.push_back(_AllJets_pt[ik1]);
        _JetsClosestMw_eta.push_back(_AllJets_eta[ik1]);
        _JetsClosestMw_phi.push_back(_AllJets_phi[ik1]);
        _JetsClosestMw_E.push_back(_AllJets_E[ik1]);
	
	_JetsClosestMw_id.push_back(2);
        _JetsClosestMw_pt.push_back(_AllJets_pt[ik2]);
        _JetsClosestMw_eta.push_back(_AllJets_eta[ik2]);
        _JetsClosestMw_phi.push_back(_AllJets_phi[ik2]);
        _JetsClosestMw_E.push_back(_AllJets_E[ik2]);
	
      }
      if (il1!=-1 && il2!=-1){
       
        _JetsLowestMjj_id.push_back(3);
        _JetsLowestMjj_pt.push_back(_AllJets_pt[il1]);
        _JetsLowestMjj_eta.push_back(_AllJets_eta[il1]);
        _JetsLowestMjj_phi.push_back(_AllJets_phi[il1]);
        _JetsLowestMjj_E.push_back(_AllJets_E[il1]);
	
	_JetsLowestMjj_id.push_back(3);
        _JetsLowestMjj_pt.push_back(_AllJets_pt[il2]);
        _JetsLowestMjj_eta.push_back(_AllJets_eta[il2]);
        _JetsLowestMjj_phi.push_back(_AllJets_phi[il2]);
        _JetsLowestMjj_E.push_back(_AllJets_E[il2]);
	
      }
   }
   
   /*
   mc_totp4_px = Ptot.Px();
   mc_totp4_py = Ptot.Py();
   mc_totp4_pt = Ptot.Pt();
   mc_met = PtotNeut.Pt();

   mc_njets25 = (*multiLepton).AllJets.size();

   (*multiLepton).Ptot = Ptot;
   (*multiLepton).mET = PtotNeut;
   */
}

void Truth::init()
{   
    _gen_PVz = -666.;
    
    _metGen_px = -666.;
    _metGen_py = -666.;
    _metGen_pt = -666.;
    _metGen_phi = -666.;
    _metGen_sumet = -666.;
    _metGen_MuonEt = -666.;
	
    _mc_truth_n = 0;

    _mc_truth_id.clear();
    _mc_truth_pt.clear();
    _mc_truth_eta.clear();
    _mc_truth_phi.clear();
    _mc_truth_E.clear();
    _mc_truth_label.clear();
    
    _gen_n = 0;
  
    _gen_pt.clear();
    _gen_eta.clear();
    _gen_phi.clear();
    _gen_m.clear();
    _gen_id.clear();
    _gen_status.clear();
    _gen_mother_id.clear();
    
    _Bjets_id.clear();
    _Leptons_id.clear();
    _Jets_id.clear();
    _AllJets_id.clear();
    _JetsHighestPt_id.clear();
    _JetsClosestMw_id.clear();
    _JetsLowestMjj_id.clear();
    _QuarksFromWs_id.clear();
    _JetsFromWs_id.clear();
        
    _Bjets_pt.clear();
    _Leptons_pt.clear();
    _Jets_pt.clear();
    _AllJets_pt.clear();
    _JetsHighestPt_pt.clear();
    _JetsClosestMw_pt.clear();
    _JetsLowestMjj_pt.clear();
    _QuarksFromWs_pt.clear();
    _JetsFromWs_pt.clear();
    
    _Bjets_eta.clear();
    _Leptons_eta.clear();
    _Jets_eta.clear();
    _AllJets_eta.clear();
    _JetsHighestPt_eta.clear();
    _JetsClosestMw_eta.clear();
    _JetsLowestMjj_eta.clear();
    _QuarksFromWs_eta.clear();
    _JetsFromWs_eta.clear();
    
    _Bjets_phi.clear();
    _Leptons_phi.clear();
    _Jets_phi.clear();
    _AllJets_phi.clear();
    _JetsHighestPt_phi.clear();
    _JetsClosestMw_phi.clear();
    _JetsLowestMjj_phi.clear();
    _QuarksFromWs_phi.clear();
    _JetsFromWs_phi.clear();
    
    _Bjets_E.clear();
    _Leptons_E.clear();
    _Jets_E.clear();
    _AllJets_E.clear();
    _JetsHighestPt_E.clear();
    _JetsClosestMw_E.clear();
    _JetsLowestMjj_E.clear();
    _QuarksFromWs_E.clear();
    _JetsFromWs_E.clear();
    
    _boson_decay = -1; // 0:Higgs->dilep, 1:Higgs->semilep, 2:Z->ll, 3:W->lnu 
    _ttbar_decay = -1; // 1:ttbar->semilep, 2:ttbar->dilep
        
}

