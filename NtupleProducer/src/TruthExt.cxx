#include "include/NtupleProducer.h"

ClassImp(TruthExt)

TruthExt::TruthExt()
{
}

TruthExt::~TruthExt()
{
}

void TruthExt::read()
{
   int UNINT = -1000;
   
   int gen_n = ntP->gen_n;
   
   bool  print_gen = false;
   
   int gen_n_sel = 0;
   for(int i=0;i<gen_n;i++)
     {
        int status = ntP->gen_status->at(i);
        if( status != 1 && status != 3 ) continue;
        gen_pt.push_back(ntP->gen_pt->at(i));
        gen_eta.push_back(ntP->gen_eta->at(i));
        gen_phi.push_back(ntP->gen_phi->at(i));
        gen_m.push_back(ntP->gen_m->at(i));
        gen_id.push_back(ntP->gen_id->at(i));
        gen_status.push_back(status);
	
        int mom = ntP->gen_mother_index->at(i);
	
        if( mom >= ntP->gen_index->size() )
	  {
	     
            std::cout << "Problem with gen length" << std::endl;
            exit(1);
	  }
	
        int mother_id = ntP->gen_id->at(mom);
	
        gen_mother_id.push_back(mother_id);
	
        gen_n_sel++;
     }
   gen_PVz = ntP->gen_PVz;
   
   mc_truth_n = 0;
   
   int mc_truth_h0_id = ntP->mc_truth_h0_id;
   if(print_gen) std::cout << "mc_truth_h0_id : " << mc_truth_h0_id << std::endl;
   if( mc_truth_h0_id != UNINT )
     {
        if(print_gen) std::cout << "is HO" << std::endl;
        mc_truth_id.push_back(mc_truth_h0_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0_E);
        mc_truth_label.push_back(1);
        mc_truth_n++;
     }
   int mc_truth_h0W1_id = ntP->mc_truth_h0W1_id;
   if(print_gen) std::cout << "mc_truth_h0W1_id : " << mc_truth_h0W1_id << std::endl;
   if( mc_truth_h0W1_id != UNINT )
     {
        if(print_gen) std::cout << "is HO to W" << std::endl;
        mc_truth_id.push_back(mc_truth_h0W1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0W1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0W1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0W1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0W1_E);
        mc_truth_label.push_back(12);
        mc_truth_n++;
     }
   int mc_truth_h0Wl1_id = ntP->mc_truth_h0Wl1_id;
   if( mc_truth_h0Wl1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wl1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wl1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wl1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wl1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wl1_E);
        mc_truth_label.push_back(120);
        mc_truth_n++;
     }
   int mc_truth_h0Wq11_id = ntP->mc_truth_h0Wq11_id;
   if( mc_truth_h0Wq11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wq11_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wq11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wq11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wq11_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wq11_E);
        mc_truth_label.push_back(122);
        mc_truth_n++;
     }
   int mc_truth_h0Wq21_id = ntP->mc_truth_h0Wq21_id;
   if( mc_truth_h0Wq21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wq21_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wq21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wq21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wq21_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wq21_E);
        mc_truth_label.push_back(123);
        mc_truth_n++;
     }
   int mc_truth_h0Wtau1_id = ntP->mc_truth_h0Wtau1_id;
   if( mc_truth_h0Wtau1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wtau1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wtau1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wtau1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wtau1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wtau1_E);
        mc_truth_label.push_back(124);
        mc_truth_n++;
     }
   int mc_truth_h0Wtaul1_id = ntP->mc_truth_h0Wtaul1_id;
   if( mc_truth_h0Wtaul1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wtaul1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wtaul1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wtaul1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wtaul1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wtaul1_E);
        mc_truth_label.push_back(1240);
        mc_truth_n++;
     }
   int mc_truth_h0W2_id = ntP->mc_truth_h0W2_id;
   if( mc_truth_h0W2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0W2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0W2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0W2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0W2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0W2_E);
        mc_truth_label.push_back(13);
        mc_truth_n++;
     }
   int mc_truth_h0Wl2_id = ntP->mc_truth_h0Wl2_id;
   if( mc_truth_h0Wl2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wl2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wl2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wl2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wl2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wl2_E);
        mc_truth_label.push_back(130);
        mc_truth_n++;
     }
   int mc_truth_h0Wq12_id = ntP->mc_truth_h0Wq12_id;
   if( mc_truth_h0Wq12_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wq12_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wq12_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wq12_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wq12_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wq12_E);
        mc_truth_label.push_back(132);
        mc_truth_n++;
     }
   int mc_truth_h0Wq22_id = ntP->mc_truth_h0Wq22_id;
   if( mc_truth_h0Wq22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wq22_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wq22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wq22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wq22_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wq22_E);
        mc_truth_label.push_back(133);
        mc_truth_n++;
     }
   int mc_truth_h0Wtau2_id = ntP->mc_truth_h0Wtau2_id;
   if( mc_truth_h0Wtau2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wtau2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wtau2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wtau2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wtau2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wtau2_E);
        mc_truth_label.push_back(134);
        mc_truth_n++;
     }
   int mc_truth_h0Wtaul2_id = ntP->mc_truth_h0Wtaul2_id;
   if( mc_truth_h0Wtaul2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Wtaul2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Wtaul2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Wtaul2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Wtaul2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Wtaul2_E);
        mc_truth_label.push_back(1340);
        mc_truth_n++;
     }
   int mc_truth_h0Z1_id = ntP->mc_truth_h0Z1_id;
   if(print_gen) std::cout << "mc_truth_h0Z1_id : " << mc_truth_h0Z1_id << std::endl;
   if( mc_truth_h0Z1_id != UNINT )
     {
        if(print_gen) std::cout << "is HO Z" << std::endl;
        mc_truth_id.push_back(mc_truth_h0Z1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Z1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Z1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Z1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Z1_E);
        mc_truth_label.push_back(14);
        mc_truth_n++;
     }
   int mc_truth_h0Zl11_id = ntP->mc_truth_h0Zl11_id;
   if( mc_truth_h0Zl11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zl11_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zl11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zl11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zl11_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zl11_E);
        mc_truth_label.push_back(140);
        mc_truth_n++;
     }
   int mc_truth_h0Zl21_id = ntP->mc_truth_h0Zl21_id;
   if( mc_truth_h0Zl21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zl21_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zl21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zl21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zl21_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zl21_E);
        mc_truth_label.push_back(141);
        mc_truth_n++;
     }
   int mc_truth_h0Zq11_id = ntP->mc_truth_h0Zq11_id;
   if( mc_truth_h0Zq11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zq11_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zq11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zq11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zq11_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zq11_E);
        mc_truth_label.push_back(142);
        mc_truth_n++;
     }
   int mc_truth_h0Zq21_id = ntP->mc_truth_h0Zq21_id;
   if( mc_truth_h0Zq21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zq21_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zq21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zq21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zq21_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zq21_E);
        mc_truth_label.push_back(143);
        mc_truth_n++;
     }
   //AC
   int mc_truth_h0Ztau11_id = ntP->mc_truth_h0Ztau11_id;
   if( mc_truth_h0Ztau11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztau11_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztau11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztau11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztau11_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztau11_E);
        mc_truth_label.push_back(144);
        mc_truth_n++;
     }
   int mc_truth_h0Ztaul11_id = ntP->mc_truth_h0Ztaul11_id;
   if( mc_truth_h0Ztaul11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztaul11_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztaul11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztaul11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztaul11_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztaul11_E);
        mc_truth_label.push_back(1440);
        mc_truth_n++;
     }
   int mc_truth_h0Ztau21_id = ntP->mc_truth_h0Ztau21_id;
   if( mc_truth_h0Ztau21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztau21_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztau21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztau21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztau21_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztau21_E);
        mc_truth_label.push_back(145);
        mc_truth_n++;
     }
   int mc_truth_h0Ztaul21_id = ntP->mc_truth_h0Ztaul21_id;
   if( mc_truth_h0Ztaul21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztaul21_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztaul21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztaul21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztaul21_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztaul21_E);
        mc_truth_label.push_back(1450);
        mc_truth_n++;
     }
   //AC
   int mc_truth_h0Z2_id = ntP->mc_truth_h0Z2_id;
   if( mc_truth_h0Z2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Z2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Z2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Z2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Z2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Z2_E);
        mc_truth_label.push_back(15);
        mc_truth_n++;
     }
   int mc_truth_h0Zl12_id = ntP->mc_truth_h0Zl12_id;
   if( mc_truth_h0Zl12_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zl12_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zl12_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zl12_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zl12_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zl12_E);
        mc_truth_label.push_back(150);
        mc_truth_n++;
     }
   int mc_truth_h0Zl22_id = ntP->mc_truth_h0Zl22_id;
   if( mc_truth_h0Zl22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zl22_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zl22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zl22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zl22_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zl22_E);
        mc_truth_label.push_back(151);
        mc_truth_n++;
     }
   int mc_truth_h0Zq12_id = ntP->mc_truth_h0Zq12_id;
   if( mc_truth_h0Zq12_id != UNINT )
     {
	mc_truth_id.push_back(mc_truth_h0Zq12_id);
	mc_truth_pt.push_back(ntP->mc_truth_h0Zq12_pt);
	mc_truth_eta.push_back(ntP->mc_truth_h0Zq12_eta);
	mc_truth_phi.push_back(ntP->mc_truth_h0Zq12_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zq12_E);
        mc_truth_label.push_back(152);
        mc_truth_n++;
     }
   int mc_truth_h0Zq22_id = ntP->mc_truth_h0Zq22_id;
   if( mc_truth_h0Zq22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Zq22_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Zq22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Zq22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Zq22_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Zq22_E);
        mc_truth_label.push_back(153);
        mc_truth_n++;
     }
   //AC
   int mc_truth_h0Ztau12_id = ntP->mc_truth_h0Ztau12_id;
   if( mc_truth_h0Ztau12_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztau12_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztau12_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztau12_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztau12_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztau12_E);
        mc_truth_label.push_back(154);
        mc_truth_n++;
     }
   int mc_truth_h0Ztaul12_id = ntP->mc_truth_h0Ztaul12_id;
   if( mc_truth_h0Ztaul12_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztaul12_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztaul12_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztaul12_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztaul12_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztaul12_E);
        mc_truth_label.push_back(1540);
        mc_truth_n++;
     }
   int mc_truth_h0Ztau22_id = ntP->mc_truth_h0Ztau22_id;
   if( mc_truth_h0Ztau22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztau22_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztau22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztau22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztau22_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztau22_E);
        mc_truth_label.push_back(155);
        mc_truth_n++;
     }
   int mc_truth_h0Ztaul22_id = ntP->mc_truth_h0Ztaul22_id;
   if( mc_truth_h0Ztaul22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0Ztaul22_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0Ztaul22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0Ztaul22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0Ztaul22_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0Ztaul22_E);
        mc_truth_label.push_back(1550);
        mc_truth_n++;
     }
   //AC
   int mc_truth_h0tau1_id = ntP->mc_truth_h0tau1_id;
   if(print_gen) std::cout << "mc_truth_h0tau1_id : " << mc_truth_h0tau1_id << std::endl;
   if( mc_truth_h0tau1_id != UNINT )
     {
        if(print_gen) std::cout << "is HO tau" << std::endl;
        mc_truth_id.push_back(mc_truth_h0tau1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0tau1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0tau1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0tau1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0tau1_E);
        mc_truth_label.push_back(16);
        mc_truth_n++;
     }
   int mc_truth_h0taul1_id = ntP->mc_truth_h0taul1_id;
   if( mc_truth_h0taul1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0taul1_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0taul1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0taul1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0taul1_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0taul1_E);
        mc_truth_label.push_back(160);
        mc_truth_n++;
     }
   int mc_truth_h0tau2_id = ntP->mc_truth_h0tau2_id;
   if( mc_truth_h0tau2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0tau2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0tau2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0tau2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0tau2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0tau2_E);
        mc_truth_label.push_back(17);
        mc_truth_n++;	
     }
   int mc_truth_h0taul2_id = ntP->mc_truth_h0taul2_id;
   if( mc_truth_h0taul2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_h0taul2_id);
        mc_truth_pt.push_back(ntP->mc_truth_h0taul2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_h0taul2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_h0taul2_phi);
        mc_truth_E.push_back(ntP->mc_truth_h0taul2_E);
        mc_truth_label.push_back(170);
        mc_truth_n++;	
     }   
   //
   int mc_truth_t1_id = ntP->mc_truth_t1_id;
   if( mc_truth_t1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_t1_id);
        mc_truth_pt.push_back(ntP->mc_truth_t1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_t1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_t1_phi);
        mc_truth_E.push_back(ntP->mc_truth_t1_E);
        mc_truth_label.push_back(2);
        mc_truth_n++;
     }
   int mc_truth_t2_id = ntP->mc_truth_t2_id;
   if( mc_truth_t2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_t2_id);
        mc_truth_pt.push_back(ntP->mc_truth_t2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_t2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_t2_phi);
        mc_truth_E.push_back(ntP->mc_truth_t2_E);
        mc_truth_label.push_back(3);
        mc_truth_n++;
     }
   int mc_truth_tb1_id = ntP->mc_truth_tb1_id;
   if(print_gen) std::cout << "mc_truth_tb1_id : " << mc_truth_tb1_id << std::endl;
   if( mc_truth_tb1_id != UNINT )
     {
        if(print_gen) std::cout << "is t b" << std::endl;
        mc_truth_id.push_back(mc_truth_tb1_id);
        if(print_gen) std::cout << "id: " << mc_truth_tb1_id << std::endl;
        mc_truth_pt.push_back(ntP->mc_truth_tb1_pt);
        if(print_gen) std::cout << "pt: " << ntP->mc_truth_tb1_pt << std::endl;
        mc_truth_eta.push_back(ntP->mc_truth_tb1_eta);
        if(print_gen) std::cout << "eta: " << ntP->mc_truth_tb1_eta << std::endl;
        mc_truth_phi.push_back(ntP->mc_truth_tb1_phi);
        if(print_gen) std::cout << "phi: " << ntP->mc_truth_tb1_phi << std::endl;
        mc_truth_E.push_back(ntP->mc_truth_tb1_E);
        if(print_gen) std::cout << "E: " << ntP->mc_truth_tb1_E << std::endl;
        mc_truth_label.push_back(20);
        mc_truth_n++;
     }
   int mc_truth_tb2_id = ntP->mc_truth_tb2_id;
   if( mc_truth_tb2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tb2_id);
        mc_truth_pt.push_back(ntP->mc_truth_tb2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tb2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tb2_phi);
        mc_truth_E.push_back(ntP->mc_truth_tb2_E);
        mc_truth_label.push_back(30);
        mc_truth_n++;
     }
   int mc_truth_tW1_id = ntP->mc_truth_tW1_id;
   if( mc_truth_tW1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tW1_id);
        mc_truth_pt.push_back(ntP->mc_truth_tW1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tW1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tW1_phi);
        mc_truth_E.push_back(ntP->mc_truth_tW1_E);
        mc_truth_label.push_back(21);
        mc_truth_n++;
     }
   int mc_truth_tW2_id = ntP->mc_truth_tW2_id;
   if( mc_truth_tW2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tW2_id);
        mc_truth_pt.push_back(ntP->mc_truth_tW2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tW2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tW2_phi);
        mc_truth_E.push_back(ntP->mc_truth_tW2_E);
        mc_truth_label.push_back(31);
        mc_truth_n++;
     }
   int mc_truth_tWl1_id = ntP->mc_truth_tWl1_id;
   if( mc_truth_tWl1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWl1_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWl1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWl1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWl1_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWl1_E);
        mc_truth_label.push_back(210);
        mc_truth_n++;
     }
   int mc_truth_tWl2_id = ntP->mc_truth_tWl2_id;
   if( mc_truth_tWl2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWl2_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWl2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWl2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWl2_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWl2_E);
        mc_truth_label.push_back(310);
        mc_truth_n++;
     }
   int mc_truth_tWtau1_id = ntP->mc_truth_tWtau1_id;
   if( mc_truth_tWtau1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWtau1_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWtau1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWtau1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWtau1_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWtau1_E);
        mc_truth_label.push_back(22);
        mc_truth_n++;
     }
   int mc_truth_tWtau2_id = ntP->mc_truth_tWtau2_id;
   if( mc_truth_tWtau2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWtau2_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWtau2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWtau2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWtau2_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWtau2_E);
        mc_truth_label.push_back(32);
        mc_truth_n++;
     }
   int mc_truth_tWtaul1_id = ntP->mc_truth_tWtaul1_id;
   if( mc_truth_tWtaul1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWtaul1_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWtaul1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWtaul1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWtaul1_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWtaul1_E);
        mc_truth_label.push_back(2220);
        mc_truth_n++;
     }
   int mc_truth_tWtaul2_id = ntP->mc_truth_tWtaul2_id;
   if( mc_truth_tWtaul2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWtaul2_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWtaul2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWtaul2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWtaul2_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWtaul2_E);
        mc_truth_label.push_back(3220);
        mc_truth_n++;
     }
   int mc_truth_tWq11_id = ntP->mc_truth_tWq11_id;
   if( mc_truth_tWq11_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWq11_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWq11_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWq11_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWq11_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWq11_E);
        mc_truth_label.push_back(223);
        mc_truth_n++;
     }
   int mc_truth_tWq21_id = ntP->mc_truth_tWq21_id;
   if( mc_truth_tWq21_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWq21_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWq21_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWq21_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWq21_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWq21_E);
        mc_truth_label.push_back(224);
        mc_truth_n++;
     }
   int mc_truth_tWq12_id = ntP->mc_truth_tWq12_id;
   if( mc_truth_tWq12_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWq12_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWq12_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWq12_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWq12_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWq12_E);
        mc_truth_label.push_back(323);
        mc_truth_n++;
     }
   int mc_truth_tWq22_id = ntP->mc_truth_tWq22_id;
   if( mc_truth_tWq22_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_tWq22_id);
        mc_truth_pt.push_back(ntP->mc_truth_tWq22_pt);
        mc_truth_eta.push_back(ntP->mc_truth_tWq22_eta);
        mc_truth_phi.push_back(ntP->mc_truth_tWq22_phi);
        mc_truth_E.push_back(ntP->mc_truth_tWq22_E);
        mc_truth_label.push_back(324);
        mc_truth_n++;
     }
   int mc_truth_W_id = ntP->mc_truth_W_id;
   if( mc_truth_W_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_W_id);
        mc_truth_pt.push_back(ntP->mc_truth_W_pt);
        mc_truth_eta.push_back(ntP->mc_truth_W_eta);
        mc_truth_phi.push_back(ntP->mc_truth_W_phi);
        mc_truth_E.push_back(ntP->mc_truth_W_E);
        mc_truth_label.push_back(4);
        mc_truth_n++;
     }
   int mc_truth_Wl_id = ntP->mc_truth_Wl_id;
   if( mc_truth_Wl_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Wl_id);
        mc_truth_pt.push_back(ntP->mc_truth_Wl_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Wl_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Wl_phi);
        mc_truth_E.push_back(ntP->mc_truth_Wl_E);
        mc_truth_label.push_back(40);
        mc_truth_n++;
     }
   int mc_truth_Wtau_id = ntP->mc_truth_Wtau_id;
   if( mc_truth_Wtau_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Wtau_id);
        mc_truth_pt.push_back(ntP->mc_truth_Wtau_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Wtau_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Wtau_phi);
        mc_truth_E.push_back(ntP->mc_truth_Wtau_E);
        mc_truth_label.push_back(43);
        mc_truth_n++;
     }
   int mc_truth_Wtaul_id = ntP->mc_truth_Wtaul_id;
   if( mc_truth_Wtaul_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Wtaul_id);
        mc_truth_pt.push_back(ntP->mc_truth_Wtaul_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Wtaul_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Wtaul_phi);
        mc_truth_E.push_back(ntP->mc_truth_Wtaul_E);
        mc_truth_label.push_back(430);
        mc_truth_n++;
     }
   int mc_truth_Wq1_id = ntP->mc_truth_Wq1_id;
   if( mc_truth_Wq1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Wq1_id);
        mc_truth_pt.push_back(ntP->mc_truth_Wq1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Wq1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Wq1_phi);
        mc_truth_E.push_back(ntP->mc_truth_Wq1_E);
        mc_truth_label.push_back(41);
        mc_truth_n++;
     }
   int mc_truth_Wq2_id = ntP->mc_truth_Wq2_id;
   if( mc_truth_Wq2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Wq2_id);
        mc_truth_pt.push_back(ntP->mc_truth_Wq2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Wq2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Wq2_phi);
        mc_truth_E.push_back(ntP->mc_truth_Wq2_E);
        mc_truth_label.push_back(42);
        mc_truth_n++;
     }
   int mc_truth_Z_id = ntP->mc_truth_Z_id;
   if( mc_truth_Z_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Z_id);
        mc_truth_pt.push_back(ntP->mc_truth_Z_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Z_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Z_phi);
        mc_truth_E.push_back(ntP->mc_truth_Z_E);
        mc_truth_label.push_back(5);
        mc_truth_n++;
     }
   int mc_truth_Zl1_id = ntP->mc_truth_Zl1_id;
   if( mc_truth_Zl1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Zl1_id);
        mc_truth_pt.push_back(ntP->mc_truth_Zl1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Zl1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Zl1_phi);
        mc_truth_E.push_back(ntP->mc_truth_Zl1_E);
        mc_truth_label.push_back(50);
        mc_truth_n++;
     }
   int mc_truth_Zl2_id = ntP->mc_truth_Zl2_id;
   if( mc_truth_Zl2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Zl2_id);
        mc_truth_pt.push_back(ntP->mc_truth_Zl2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Zl2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Zl2_phi);
        mc_truth_E.push_back(ntP->mc_truth_Zl2_E);
        mc_truth_label.push_back(51);
        mc_truth_n++;
     }
   int mc_truth_Ztau1_id = ntP->mc_truth_Ztau1_id;
   if( mc_truth_Ztau1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Ztau1_id);
        mc_truth_pt.push_back(ntP->mc_truth_Ztau1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Ztau1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Ztau1_phi);
        mc_truth_E.push_back(ntP->mc_truth_Ztau1_E);
        mc_truth_label.push_back(52);
        mc_truth_n++;
     }
   int mc_truth_Ztau2_id = ntP->mc_truth_Ztau2_id;
   if( mc_truth_Ztau2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Ztau2_id);
        mc_truth_pt.push_back(ntP->mc_truth_Ztau2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Ztau2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Ztau2_phi);
        mc_truth_E.push_back(ntP->mc_truth_Ztau2_E);
        mc_truth_label.push_back(53);
        mc_truth_n++;
     }
   int mc_truth_Ztaul1_id = ntP->mc_truth_Ztaul1_id;
   if( mc_truth_Ztaul1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Ztaul1_id);
        mc_truth_pt.push_back(ntP->mc_truth_Ztaul1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Ztaul1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Ztaul1_phi);
        mc_truth_E.push_back(ntP->mc_truth_Ztaul1_E);
        mc_truth_label.push_back(520);
        mc_truth_n++;
     }
   int mc_truth_Ztaul2_id = ntP->mc_truth_Ztaul2_id;
   if( mc_truth_Ztaul2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Ztaul2_id);
        mc_truth_pt.push_back(ntP->mc_truth_Ztaul2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Ztaul2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Ztaul2_phi);
        mc_truth_E.push_back(ntP->mc_truth_Ztaul2_E);
        mc_truth_label.push_back(530);
        mc_truth_n++;
     }
   int mc_truth_Zq1_id = ntP->mc_truth_Zq1_id;
   if( mc_truth_Zq1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Zq1_id);
        mc_truth_pt.push_back(ntP->mc_truth_Zq1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Zq1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Zq1_phi);
        mc_truth_E.push_back(ntP->mc_truth_Zq1_E);
        mc_truth_label.push_back(54);
        mc_truth_n++;
     }
   int mc_truth_Zq2_id = ntP->mc_truth_Zq2_id;
   if( mc_truth_Zq2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_Zq2_id);
        mc_truth_pt.push_back(ntP->mc_truth_Zq2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_Zq2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_Zq2_phi);
        mc_truth_E.push_back(ntP->mc_truth_Zq2_E);
        mc_truth_label.push_back(55);
        mc_truth_n++;
     }
   int mc_truth_gamma_id = ntP->mc_truth_gamma_id;
   if( mc_truth_gamma_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gamma_id);
        mc_truth_pt.push_back(ntP->mc_truth_gamma_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gamma_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gamma_phi);
        mc_truth_E.push_back(ntP->mc_truth_gamma_E);
        mc_truth_label.push_back(6);
        mc_truth_n++;
     }
   int mc_truth_gammal1_id = ntP->mc_truth_gammal1_id;
   if( mc_truth_gammal1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammal1_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammal1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammal1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammal1_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammal1_E);
        mc_truth_label.push_back(60);
        mc_truth_n++;
     }
   int mc_truth_gammal2_id = ntP->mc_truth_gammal2_id;
   if( mc_truth_gammal2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammal2_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammal2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammal2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammal2_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammal2_E);
        mc_truth_label.push_back(61);
        mc_truth_n++;
     }
   int mc_truth_gammatau1_id = ntP->mc_truth_gammatau1_id;
   if( mc_truth_gammatau1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammatau1_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammatau1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammatau1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammatau1_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammatau1_E);
        mc_truth_label.push_back(62);
        mc_truth_n++;
     }
   int mc_truth_gammatau2_id = ntP->mc_truth_gammatau2_id;
   if( mc_truth_gammatau2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammatau2_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammatau2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammatau2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammatau2_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammatau2_E);
        mc_truth_label.push_back(63);
        mc_truth_n++;
     }
   int mc_truth_gammataul1_id = ntP->mc_truth_gammataul1_id;
   if( mc_truth_gammataul1_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammataul1_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammataul1_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammataul1_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammataul1_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammataul1_E);
        mc_truth_label.push_back(620);
        mc_truth_n++;
     }
   int mc_truth_gammataul2_id = ntP->mc_truth_gammataul2_id;
   if( mc_truth_gammataul2_id != UNINT )
     {
        mc_truth_id.push_back(mc_truth_gammataul2_id);
        mc_truth_pt.push_back(ntP->mc_truth_gammataul2_pt);
        mc_truth_eta.push_back(ntP->mc_truth_gammataul2_eta);
        mc_truth_phi.push_back(ntP->mc_truth_gammataul2_phi);
        mc_truth_E.push_back(ntP->mc_truth_gammataul2_E);
        mc_truth_label.push_back(630);
        mc_truth_n++;
     }

   metGen_px = ntP->metGen_px;
   metGen_py = ntP->metGen_py;
   metGen_pt = ntP->metGen_pt;
   metGen_phi = ntP->metGen_phi;
   metGen_sumet = ntP->metGen_sumet;
   metGen_MuonEt = ntP->metGen_MuonEt;   
}

void TruthExt::readMultiLepton()
{  
   int UNINT = -100;
   
   if (ntP->mc_truth_t1_id==-100 || ntP->mc_truth_t2_id==-100)  return;
   
   if (ntP->mc_truth_h0_id>-100 && ntP->mc_truth_h0Wl1_id>-100 && ntP->mc_truth_h0Wl2_id>-100)
     {
        Leptons_id.push_back(ntP->mc_truth_h0Wl1_id);
        Leptons_id.push_back(ntP->mc_truth_h0Wl2_id);
        Leptons_pt.push_back(ntP->mc_truth_h0Wl1_pt);
        Leptons_pt.push_back(ntP->mc_truth_h0Wl2_pt);
        Leptons_eta.push_back(ntP->mc_truth_h0Wl1_eta);
        Leptons_eta.push_back(ntP->mc_truth_h0Wl2_eta);
        Leptons_phi.push_back(ntP->mc_truth_h0Wl1_phi);
        Leptons_phi.push_back(ntP->mc_truth_h0Wl2_phi);
        Leptons_E.push_back(ntP->mc_truth_h0Wl1_E);
        Leptons_E.push_back(ntP->mc_truth_h0Wl2_E);
	
        boson_decay = 0;
     }      
   else if (ntP->mc_truth_h0_id>-100 && ntP->mc_truth_h0Wq11_id>-100 && ntP->mc_truth_h0Wl2_id>-100)
     {
        Leptons_id.push_back(ntP->mc_truth_h0Wl2_id);
        Leptons_pt.push_back(ntP->mc_truth_h0Wl2_pt);
        Leptons_eta.push_back(ntP->mc_truth_h0Wl2_eta);
        Leptons_phi.push_back(ntP->mc_truth_h0Wl2_phi);
        Leptons_E.push_back(ntP->mc_truth_h0Wl2_E);
	
        QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq11_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq11_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq11_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq11_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq11_E);

        QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq12_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq12_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq12_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq12_E);

        boson_decay = 1;   

        float dr1 = 999; 
        float dr2 = 999;
        float drMax1 = 0.4;
        float drMax2 = 0.4;
        int iMax1 = -1;
        int iMax2 = -1;
	
        for (int i=0; i<nt->NtGenJet->size(); i++)
	  {
	     dr1 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_h0Wq11_eta, ntP->mc_truth_h0Wq11_phi);
	     dr2 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_h0Wq21_eta, ntP->mc_truth_h0Wq21_phi);
	     if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
	     if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
	  }
	
        if (iMax1 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_h0Wq11_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E);
	  }  
        if (iMax2 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_h0Wq21_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E);
	  }
     }   
   else if (ntP->mc_truth_h0_id>-100 && ntP->mc_truth_h0Wq22_id>-100 && ntP->mc_truth_h0Wl1_id>-100)
     {
	
        Leptons_id.push_back(ntP->mc_truth_h0Wl1_id);
        Leptons_pt.push_back(ntP->mc_truth_h0Wl1_pt);
        Leptons_eta.push_back(ntP->mc_truth_h0Wl1_eta);
        Leptons_phi.push_back(ntP->mc_truth_h0Wl1_phi);
        Leptons_E.push_back(ntP->mc_truth_h0Wl1_E);
	
        QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq12_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq12_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq12_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq12_E);

        QuarksFromWs_id.push_back(ntP->mc_truth_h0Wq22_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_h0Wq22_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_h0Wq22_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_h0Wq22_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_h0Wq22_E);

        boson_decay = 1;

        float dr1 = 999; 
        float dr2 = 999;
        float drMax1 = 0.4;
        float drMax2 = 0.4;
        int iMax1 = -1;
        int iMax2 = -1;

        for (int i=0; i<nt->NtGenJet->size(); i++)
	  {
	     dr1 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_h0Wq12_eta, ntP->mc_truth_h0Wq12_phi);
	     dr2 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_h0Wq22_eta, ntP->mc_truth_h0Wq22_phi);
	     if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
	     if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
	  }
	
        if (iMax1 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_h0Wq12_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E);
	  }
        if (iMax2 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_h0Wq22_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E);
	  }
	
     }      
   else if (ntP->mc_truth_W_id>-100 && ntP->mc_truth_Wl_id>-100)
     {
        Leptons_id.push_back(ntP->mc_truth_Wl_id);
        Leptons_pt.push_back(ntP->mc_truth_Wl_pt);
        Leptons_eta.push_back(ntP->mc_truth_Wl_eta);
        Leptons_phi.push_back(ntP->mc_truth_Wl_phi);
        Leptons_E.push_back(ntP->mc_truth_Wl_E);
	
        boson_decay = 3;    
     }   
   else if (ntP->mc_truth_Z_id>-100 && ntP->mc_truth_Zl1_id>-100 && ntP->mc_truth_Zl2_id>-100)
     {
	Leptons_id.push_back(ntP->mc_truth_Zl1_id);      
	Leptons_pt.push_back(ntP->mc_truth_Zl1_pt);
        Leptons_eta.push_back(ntP->mc_truth_Zl1_eta);
        Leptons_phi.push_back(ntP->mc_truth_Zl1_phi);
        Leptons_E.push_back(ntP->mc_truth_Zl1_E);

        Leptons_id.push_back(ntP->mc_truth_Zl2_id);      
        Leptons_pt.push_back(ntP->mc_truth_Zl2_pt);
        Leptons_eta.push_back(ntP->mc_truth_Zl2_eta);
        Leptons_phi.push_back(ntP->mc_truth_Zl2_phi);
        Leptons_E.push_back(ntP->mc_truth_Zl2_E);

        boson_decay = 2;      
     }
   else if (ntP->mc_truth_gammal1_id>-100 && ntP->mc_truth_gammal2_id>-100)
     {
        Leptons_id.push_back(ntP->mc_truth_gammal1_id);      
        Leptons_pt.push_back(ntP->mc_truth_gammal1_pt);
        Leptons_eta.push_back(ntP->mc_truth_gammal1_eta);
        Leptons_phi.push_back(ntP->mc_truth_gammal1_phi);
        Leptons_E.push_back(ntP->mc_truth_gammal1_E);

        Leptons_id.push_back(ntP->mc_truth_gammal2_id);      
        Leptons_pt.push_back(ntP->mc_truth_gammal2_pt);
        Leptons_eta.push_back(ntP->mc_truth_gammal2_eta);
        Leptons_phi.push_back(ntP->mc_truth_gammal2_phi);
        Leptons_E.push_back(ntP->mc_truth_gammal2_E);

        boson_decay = 2;    
     }   
   else return;

   if (ntP->mc_truth_tWq11_id>-100 && ntP->mc_truth_tWl2_id>-100)
     {
        Leptons_id.push_back(ntP->mc_truth_tWl2_id);      
        Leptons_pt.push_back(ntP->mc_truth_tWl2_pt);
        Leptons_eta.push_back(ntP->mc_truth_tWl2_eta);
        Leptons_phi.push_back(ntP->mc_truth_tWl2_phi);
        Leptons_E.push_back(ntP->mc_truth_tWl2_E);

        Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
        Bjets_E.push_back(ntP->mc_truth_tb1_E);

        Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
        Bjets_E.push_back(ntP->mc_truth_tb2_E);

        QuarksFromWs_id.push_back(ntP->mc_truth_tWq11_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_tWq11_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_tWq11_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_tWq11_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_tWq11_E);

        QuarksFromWs_id.push_back(ntP->mc_truth_tWq21_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_tWq21_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_tWq21_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_tWq21_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_tWq21_E);

        ttbar_decay = 1;

        float dr1 = 999; 
        float dr2 = 999;
        float drMax1 = 0.4;
        float drMax2 = 0.4;
        int iMax1 = -1;
        int iMax2 = -1;

        for (int i=0; i<nt->NtGenJet->size(); i++)
	  {
	     dr1 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_tWq11_eta, ntP->mc_truth_tWq11_phi);
	     dr2 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_tWq22_eta, ntP->mc_truth_tWq22_phi);
	     if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
	     if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
	  }
	
        if (iMax1 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_tWq11_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E);
	  }
        if (iMax2 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_tWq22_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E);
	  }
     }   
   else if (ntP->mc_truth_tWq22_id>-100 && ntP->mc_truth_tWl1_id>-100)
     {      
        Leptons_id.push_back(ntP->mc_truth_tWl1_id);      
        Leptons_pt.push_back(ntP->mc_truth_tWl1_pt);
        Leptons_eta.push_back(ntP->mc_truth_tWl1_eta);
        Leptons_phi.push_back(ntP->mc_truth_tWl1_phi);
        Leptons_E.push_back(ntP->mc_truth_tWl1_E);
	
        Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
        Bjets_E.push_back(ntP->mc_truth_tb1_E);

        Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
        Bjets_E.push_back(ntP->mc_truth_tb2_E);

        QuarksFromWs_id.push_back(ntP->mc_truth_tWq12_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_tWq12_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_tWq12_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_tWq12_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_tWq12_E);
	
        QuarksFromWs_id.push_back(ntP->mc_truth_tWq22_id);
        QuarksFromWs_pt.push_back(ntP->mc_truth_tWq22_pt);
        QuarksFromWs_eta.push_back(ntP->mc_truth_tWq22_eta);
        QuarksFromWs_phi.push_back(ntP->mc_truth_tWq22_phi);
        QuarksFromWs_E.push_back(ntP->mc_truth_tWq22_E);

        float dr1 = 999; 
        float dr2 = 999;
        float drMax1 = 0.4;
        float drMax2 = 0.4;
        int iMax1 = -1;
        int iMax2 = -1;

        for (int i=0; i<nt->NtGenJet->size(); i++)
	  {
	     dr1 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_tWq12_eta, ntP->mc_truth_tWq12_phi);
	     dr2 = GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi, ntP->mc_truth_tWq22_eta, ntP->mc_truth_tWq22_phi);
	     if (dr1 < drMax1) {iMax1 = i; drMax1 = dr1;}
	     if (dr2 < drMax2) {iMax2 = i; drMax2 = dr2;}
	  }
	
        if (iMax1 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_tWq12_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax1).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax1).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax1).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax1).E);
	  }
        if (iMax2 == -1)
	  {
	     JetsFromWs_id.push_back(UNINT);
	     JetsFromWs_pt.push_back(UNINT);
	     JetsFromWs_eta.push_back(UNINT);
	     JetsFromWs_phi.push_back(UNINT);
	     JetsFromWs_E.push_back(UNINT);
	  }
        else 
	  {
	     JetsFromWs_id.push_back(ntP->mc_truth_tWq22_id);
	     JetsFromWs_pt.push_back(nt->NtGenJet->at(iMax2).pt);
	     JetsFromWs_eta.push_back(nt->NtGenJet->at(iMax2).eta);
	     JetsFromWs_phi.push_back(nt->NtGenJet->at(iMax2).phi);
	     JetsFromWs_E.push_back(nt->NtGenJet->at(iMax2).E);
	  }
	
        ttbar_decay = 1;      
     }   
   else if (ntP->mc_truth_tWl1_id>-100 && ntP->mc_truth_tWl2_id>-100 && (boson_decay==1 || boson_decay==3))
     {
        Bjets_id.push_back(ntP->mc_truth_tb1_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb1_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb1_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb1_phi);
        Bjets_E.push_back(ntP->mc_truth_tb1_E);

        Bjets_id.push_back(ntP->mc_truth_tb2_id);	  
        Bjets_pt.push_back(ntP->mc_truth_tb2_pt);
        Bjets_eta.push_back(ntP->mc_truth_tb2_eta);
        Bjets_phi.push_back(ntP->mc_truth_tb2_phi);
        Bjets_E.push_back(ntP->mc_truth_tb2_E);

        Leptons_id.push_back(ntP->mc_truth_tWl1_id);      
        Leptons_pt.push_back(ntP->mc_truth_tWl1_pt);
        Leptons_eta.push_back(ntP->mc_truth_tWl1_eta);
        Leptons_phi.push_back(ntP->mc_truth_tWl1_phi);
        Leptons_E.push_back(ntP->mc_truth_tWl1_E);

        Leptons_id.push_back(ntP->mc_truth_tWl2_id);      
        Leptons_pt.push_back(ntP->mc_truth_tWl2_pt);
        Leptons_eta.push_back(ntP->mc_truth_tWl2_eta);
        Leptons_phi.push_back(ntP->mc_truth_tWl2_phi);
        Leptons_E.push_back(ntP->mc_truth_tWl2_E);     

        ttbar_decay = 2;
     }   
    else return;
   
   for (int i=0; i<nt->NtGenJet->size(); i++)
     { 
        if ( fabs(nt->NtGenJet->at(i).eta)>2.5 ) continue;
        if ( nt->NtGenJet->at(i).pt<25 ) continue;
	
        if ( GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi,Bjets_eta[0],Bjets_phi[0]) < 0.4 ||
	     GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi,Bjets_eta[1],Bjets_phi[1]) < 0.4 || 
	     GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi,Leptons_eta[0],Leptons_phi[0]) < 0.4 ||
	     GetDeltaR(nt->NtGenJet->at(i).eta,nt->NtGenJet->at(i).phi,Leptons_eta[1],Leptons_phi[1]) < 0.4 ) continue;
	
        AllJets_id.push_back(0);
        AllJets_pt.push_back(nt->NtGenJet->at(i).pt);
        AllJets_eta.push_back(nt->NtGenJet->at(i).eta);
        AllJets_phi.push_back(nt->NtGenJet->at(i).phi);
        AllJets_E.push_back(nt->NtGenJet->at(i).E);	
     }
   
   if (AllJets_E.size()==1) 
     {       
	Jets_id.push_back(0);
	Jets_pt.push_back(AllJets_pt[0]);
	Jets_eta.push_back(AllJets_eta[0]);
	Jets_phi.push_back(AllJets_phi[0]);
	Jets_E.push_back(AllJets_E[0]);
     }
   else if (AllJets_E.size()>=2) 
     {
        //mc_3l_category = 2;
	
        float pt_max=0, pt_max2=0; int ij1=-1, ij2=-1;
        float diffmass_min = 10000, mass_min = 10000; int ik1=-1, ik2=-1, il1=-1, il2=-1;
	
        for (unsigned int ij=0; ij<AllJets_E.size(); ij++){
	   if (AllJets_pt[ij] > pt_max ) {
	      pt_max2 = pt_max;
	      ij2 = ij1;
	      pt_max = AllJets_pt[ij];
	      ij1 = ij;
	   }
	   if ( AllJets_pt[ij] < pt_max && AllJets_pt[ij] > pt_max2){
	      pt_max2 = AllJets_pt[ij];
	      ij2 = ij;
	   }
	   for (unsigned int ik=0; ik<AllJets_E.size(); ik++){
	      if (ik==ij) continue;
	      
	      TLorentzVector Jij;
	      Jij.SetPtEtaPhiE(AllJets_pt[ij], AllJets_eta[ij], AllJets_phi[ij], AllJets_E[ij]);
	      TLorentzVector Jik;
	      Jik.SetPtEtaPhiE(AllJets_pt[ik], AllJets_eta[ik], AllJets_phi[ik], AllJets_E[ik]);

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
	   
	   JetsHighestPt_id.push_back(1);
	   JetsHighestPt_pt.push_back(AllJets_pt[ij1]);
	   JetsHighestPt_eta.push_back(AllJets_eta[ij1]);
	   JetsHighestPt_phi.push_back(AllJets_phi[ij1]);
	   JetsHighestPt_E.push_back(AllJets_E[ij1]);
	   
	   JetsHighestPt_id.push_back(1);
	   JetsHighestPt_pt.push_back(AllJets_pt[ij2]);
	   JetsHighestPt_eta.push_back(AllJets_eta[ij2]);
	   JetsHighestPt_phi.push_back(AllJets_phi[ij2]);
	   JetsHighestPt_E.push_back(AllJets_E[ij2]);
	   
        }
        if (ik1!=-1 && ik2!=-1){
	   
	   JetsClosestMw_id.push_back(2);
	   JetsClosestMw_pt.push_back(AllJets_pt[ik1]);
	   JetsClosestMw_eta.push_back(AllJets_eta[ik1]);
	   JetsClosestMw_phi.push_back(AllJets_phi[ik1]);
	   JetsClosestMw_E.push_back(AllJets_E[ik1]);
	   
	   JetsClosestMw_id.push_back(2);
	   JetsClosestMw_pt.push_back(AllJets_pt[ik2]);
	   JetsClosestMw_eta.push_back(AllJets_eta[ik2]);
	   JetsClosestMw_phi.push_back(AllJets_phi[ik2]);
	   JetsClosestMw_E.push_back(AllJets_E[ik2]);
	   
        }
        if (il1!=-1 && il2!=-1){

	   JetsLowestMjj_id.push_back(3);
	   JetsLowestMjj_pt.push_back(AllJets_pt[il1]);
	   JetsLowestMjj_eta.push_back(AllJets_eta[il1]);
	   JetsLowestMjj_phi.push_back(AllJets_phi[il1]);
	   JetsLowestMjj_E.push_back(AllJets_E[il1]);
	   
	   JetsLowestMjj_id.push_back(3);
	   JetsLowestMjj_pt.push_back(AllJets_pt[il2]);
	   JetsLowestMjj_eta.push_back(AllJets_eta[il2]);
	   JetsLowestMjj_phi.push_back(AllJets_phi[il2]);
	   JetsLowestMjj_E.push_back(AllJets_E[il2]);
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

void TruthExt::init()
{
   gen_PVz = -100.;
   
   metGen_px = -100.;
   metGen_py = -100.;
   metGen_pt = -100.;
   metGen_phi = -100.;
   metGen_sumet = -100.;
   metGen_MuonEt = -100.;
   
   mc_truth_n = 0;
   
   mc_truth_id.clear();
   mc_truth_pt.clear();
   mc_truth_eta.clear();
   mc_truth_phi.clear();
   mc_truth_E.clear();
   mc_truth_label.clear();
   
   gen_n = 0;
   
   gen_pt.clear();
   gen_eta.clear();
   gen_phi.clear();
   gen_m.clear();
   gen_id.clear();
   gen_status.clear();
   gen_mother_id.clear();
   
   Bjets_id.clear();
   Leptons_id.clear();
   Jets_id.clear();
   AllJets_id.clear();
   JetsHighestPt_id.clear();
   JetsClosestMw_id.clear();
   JetsLowestMjj_id.clear();
   QuarksFromWs_id.clear();
   JetsFromWs_id.clear();
   
   Bjets_pt.clear();
   Leptons_pt.clear();
   Jets_pt.clear();
   AllJets_pt.clear();
   JetsHighestPt_pt.clear();
   JetsClosestMw_pt.clear();
   JetsLowestMjj_pt.clear();
   QuarksFromWs_pt.clear();
   JetsFromWs_pt.clear();
   
   Bjets_eta.clear();
   Leptons_eta.clear();
   Jets_eta.clear();
   AllJets_eta.clear();
   JetsHighestPt_eta.clear();
   JetsClosestMw_eta.clear();
   JetsLowestMjj_eta.clear();
   QuarksFromWs_eta.clear();
   JetsFromWs_eta.clear();
   
   Bjets_phi.clear();
   Leptons_phi.clear();
   Jets_phi.clear();
   AllJets_phi.clear();
   JetsHighestPt_phi.clear();
   JetsClosestMw_phi.clear();
   JetsLowestMjj_phi.clear();
   QuarksFromWs_phi.clear();
   JetsFromWs_phi.clear();
   
   Bjets_E.clear();
   Leptons_E.clear();
   Jets_E.clear();
   AllJets_E.clear();
   JetsHighestPt_E.clear();
   JetsClosestMw_E.clear();
   JetsLowestMjj_E.clear();
   QuarksFromWs_E.clear();
   JetsFromWs_E.clear();
   
   boson_decay = -1; // 0:Higgs->dilep, 1:Higgs->semilep, 2:Z->ll, 3:W->lnu 
   ttbar_decay = -1; // 1:ttbar->semilep, 2:ttbar->dilep
}

