#ifndef HIST_H
#define HIST_H

#include "Analyzer.h"

class Hist
{
   
 public:
   
   Hist();
   virtual ~Hist();
   
   void setElectron(std::vector<Electron> *v)                  {_v_Electron = v;};
   void setMuon(std::vector<Muon> *v)                      {_v_Muon = v;};
   void setJet(std::vector<Jet> *v)                      {_v_Jet = v;};
   void setEvent(std::vector<Event> *v)                      {_v_Event = v;};
	
   void init();
   void fill();
   void close();
   bool printout(bool doPrint);

   std::vector<int> filterPt(double Pt1, double Pt2, double Pt3,
			     std::vector<Electron>* ntElectron,
			     std::vector<Muon>* ntMuon);

   unsigned int passSel;
   
   std::vector<TLorentzVector> lep_tlv(std::string sys,int llc,bool isTIGHT);
   std::vector<TLorentzVector> jet_tlv(std::string sys); 
	
   std::string dilchan[100];                                                                                           
   std::string dilchar[100];                                                                                           
   std::string jets[100];                                                                                              
   std::string sys[100];                                                                                               
   std::string sys_low[100];                                                                                           
   std::string sys_up[100];                                                                                            
   std::string histname_dilep[100];                                                                                    
   std::string histname_lep[100];                                                                                      
   std::string histname_jet[100];                                                                                      
   std::string type[100];                                                                                              
   std::string sel[100];                                                                                               
   
   std::string histname_lep_2d[100];                                                                                   
   
   std::vector<std::string> hname;
	
   int sel_n;
   int type_n;
   int sys_n;                                                                                                          
   int sys_low_n;
   int sys_up_n;
   int n_jets;
	
   bool fillThis;

   // Tree variables
   double m_weight;
   int m_channel;
   int m_sel;
   
   int m_l1_type;
   int m_l1_charge;
   double m_l1_pt;
   double m_l1_eta;
   double m_l1_phi;
   double m_l1_m;
   
   int m_l2_type;
   int m_l2_charge;
   double m_l2_pt;
   double m_l2_eta;
   double m_l2_phi;
   double m_l2_m;

   int m_l3_type;
   int m_l3_charge;
   double m_l3_pt;
   double m_l3_eta;
   double m_l3_phi;
   double m_l3_m;

   int m_l4_type;
   int m_l4_charge;
   double m_l4_pt;
   double m_l4_eta;
   double m_l4_phi;
   double m_l4_m;
   
   int m_j1_ntrk;
   double m_j1_btag;
   double m_j1_pt;
   double m_j1_eta;
   double m_j1_phi;
   double m_j1_m;

   int m_j2_ntrk;
   double m_j2_btag;
   double m_j2_pt;
   double m_j2_eta;
   double m_j2_phi;
   double m_j2_m;

   int m_j3_ntrk;
   double m_j3_btag;
   double m_j3_pt;
   double m_j3_eta;
   double m_j3_phi;
   double m_j3_m;

   int m_j4_ntrk;
   double m_j4_btag;
   double m_j4_pt;
   double m_j4_eta;
   double m_j4_phi;
   double m_j4_m;

   int m_j5_ntrk;
   double m_j5_btag;
   double m_j5_pt;
   double m_j5_eta;
   double m_j5_phi;
   double m_j5_m;

   int m_j6_ntrk;
   double m_j6_btag;
   double m_j6_pt;
   double m_j6_eta;
   double m_j6_phi;
   double m_j6_m;
   
   double m_met;
   double m_met_phi;
   
 protected:

   std::vector<Electron>             *_v_Electron;
   std::vector<Muon>             *_v_Muon;
   std::vector<Electron>             *_v_ElectronTight;
   std::vector<Muon>             *_v_MuonTight;
   
   std::vector<Event>             *_v_Event;
   std::vector<Jet>             *_v_Jet;
   std::vector<Jet>             *_v_JetTight;
  
   std::map<std::string, TH1D*> *_m1d_Event;
   std::map<std::string, TH1D*> *_m1d_Dilepton;
   std::map<std::string, TH1D*> *_m1d_Lepton;
   std::map<std::string, TH1D*> *_m1d_Jet;

   std::vector<std::pair<std::vector<std::string>,double*> > *_s_Dilepton;
	
 private:

   // [FLAVOUR][PT][ETA][BTAG][ADDSEL][VAR][2*(NSYS-1)+1]
   std::string histNAMES[5][19][1][10][1][10][25];
   std::string histNAMES_2d[5][19][1][10][1][10][25];
   std::string histNAMES_3d[5][19][1][10][1][10][25];
   
   double cJER[5];
   double cJER_down[5];
   double cJER_up[5];
	
   FILE *_fevc;
   std::ofstream _fevcVal;
	
   TRandom3 *rnd;
   
   TFile *_fout;
   TTree *_trout;
};

#endif
