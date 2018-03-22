#include "../include/Hist.h"
#include "../include/Ranges.h"

#include <boost/bind.hpp>                                                                                                   
#include <boost/function.hpp>                                                                                               
#include <boost/type_traits.hpp>

Hist::Hist()
{
    std::string foutlog = "log/list.txt";
    _fevc = fopen(foutlog.c_str(),"w");
    std::string foutlogVal = "log/list.val";
    _fevcVal.open(foutlogVal.c_str());

    _v_ElectronTight = new std::vector<Electron>;
    _v_MuonTight = new std::vector<Muon>;
    _v_JetTight = new std::vector<Jet>;
}

Hist::~Hist()
{
    delete _v_ElectronTight;
    delete _v_MuonTight;
    delete _v_JetTight;
}

void Hist::init()
{
    rnd = new TRandom3(666);

    _fout = new TFile("hist/output.root","RECREATE");

    _trout = new TTree("tr","tree");

    _trout->Branch("weight",&m_weight,"weight/D");
    _trout->Branch("channel",&m_channel,"channel/I");
    _trout->Branch("sel",&m_sel,"sel/I");

    _trout->Branch("l1_type",&m_l1_type,"l1_type/I");
    _trout->Branch("l1_charge",&m_l1_charge,"l1_charge/I");                                                  
    _trout->Branch("l1_pt",&m_l1_pt,"l1_pt/D");
    _trout->Branch("l1_eta",&m_l1_eta,"l1_eta/D");
    _trout->Branch("l1_phi",&m_l1_phi,"l1_phi/D");
    _trout->Branch("l1_m",&m_l1_m,"l1_m/D");

    _trout->Branch("l2_type",&m_l2_type,"l2_type/I");
    _trout->Branch("l2_charge",&m_l2_charge,"l2_charge/I");
    _trout->Branch("l2_pt",&m_l2_pt,"l2_pt/D");
    _trout->Branch("l2_eta",&m_l2_eta,"l2_eta/D");
    _trout->Branch("l2_phi",&m_l2_phi,"l2_phi/D");
    _trout->Branch("l2_m",&m_l2_m,"l2_m/D");

    _trout->Branch("l3_type",&m_l3_type,"l3_type/I");
    _trout->Branch("l3_charge",&m_l3_charge,"l3_charge/I");
    _trout->Branch("l3_pt",&m_l3_pt,"l3_pt/D");
    _trout->Branch("l3_eta",&m_l3_eta,"l3_eta/D");
    _trout->Branch("l3_phi",&m_l3_phi,"l3_phi/D");
    _trout->Branch("l3_m",&m_l3_m,"l3_m/D");

    _trout->Branch("l4_type",&m_l4_type,"l4_type/I");
    _trout->Branch("l4_charge",&m_l4_charge,"l4_charge/I");
    _trout->Branch("l4_pt",&m_l4_pt,"l4_pt/D");
    _trout->Branch("l4_eta",&m_l4_eta,"l4_eta/D");
    _trout->Branch("l4_phi",&m_l4_phi,"l4_phi/D");
    _trout->Branch("l4_m",&m_l4_m,"l4_m/D");

    _trout->Branch("j1_ntrk",&m_j1_ntrk,"j1_ntrk/I");
    _trout->Branch("j1_btag",&m_j1_btag,"j1_btag/D");
    _trout->Branch("j1_pt",&m_j1_pt,"j1_pt/D");
    _trout->Branch("j1_eta",&m_j1_eta,"j1_eta/D");
    _trout->Branch("j1_phi",&m_j1_phi,"j1_phi/D");
    _trout->Branch("j1_m",&m_j1_m,"j1_m/D");

    _trout->Branch("j2_ntrk",&m_j2_ntrk,"j2_ntrk/I");                                                        
    _trout->Branch("j2_btag",&m_j2_btag,"j2_btag/D");                                                        
    _trout->Branch("j2_pt",&m_j2_pt,"j2_pt/D");                                                              
    _trout->Branch("j2_eta",&m_j2_eta,"j2_eta/D");                                                           
    _trout->Branch("j2_phi",&m_j2_phi,"j2_phi/D");                                                           
    _trout->Branch("j2_m",&m_j2_m,"j2_m/D");

    _trout->Branch("j3_ntrk",&m_j3_ntrk,"j3_ntrk/I");                                                        
    _trout->Branch("j3_btag",&m_j3_btag,"j3_btag/D");
    _trout->Branch("j3_pt",&m_j3_pt,"j3_pt/D");                                                              
    _trout->Branch("j3_eta",&m_j3_eta,"j3_eta/D");                                                           
    _trout->Branch("j3_phi",&m_j3_phi,"j3_phi/D");                                                           
    _trout->Branch("j3_m",&m_j3_m,"j3_m/D");

    _trout->Branch("j4_ntrk",&m_j4_ntrk,"j4_ntrk/I");
    _trout->Branch("j4_btag",&m_j4_btag,"j4_btag/D");
    _trout->Branch("j4_pt",&m_j4_pt,"j4_pt/D");                                                              
    _trout->Branch("j4_eta",&m_j4_eta,"j4_eta/D");                                                           
    _trout->Branch("j4_phi",&m_j4_phi,"j4_phi/D");                                                           
    _trout->Branch("j4_m",&m_j4_m,"j4_m/D");

    _trout->Branch("j5_ntrk",&m_j5_ntrk,"j5_ntrk/I");
    _trout->Branch("j5_btag",&m_j5_btag,"j5_btag/D");
    _trout->Branch("j5_pt",&m_j5_pt,"j5_pt/D");                                                              
    _trout->Branch("j5_eta",&m_j5_eta,"j5_eta/D");                                                           
    _trout->Branch("j5_phi",&m_j5_phi,"j5_phi/D");                                                           
    _trout->Branch("j5_m",&m_j5_m,"j5_m/D");

    _trout->Branch("j6_ntrk",&m_j6_ntrk,"j6_ntrk/I");
    _trout->Branch("j6_btag",&m_j6_btag,"j6_btag/D");
    _trout->Branch("j6_pt",&m_j6_pt,"j6_pt/D");                                                              
    _trout->Branch("j6_eta",&m_j6_eta,"j6_eta/D");                                                           
    _trout->Branch("j6_phi",&m_j6_phi,"j6_phi/D");                                                           
    _trout->Branch("j6_m",&m_j6_m,"j6_m/D");

    _trout->Branch("met_phi",&m_met_phi,"met_phi/D");
    _trout->Branch("met",&m_met,"met/D");

    std::cout << "Initialisation done" << std::endl;
}

void Hist::fill()
{   
}

void Hist::close()
{
    _fout->Write();
    _fout->Close();
    _fevcVal.close();
}

bool Hist::printout(bool doPrint)
{
    m_weight = -666;

    m_channel = -666;
    m_sel = -666;

    m_l1_type = -666;
    m_l1_charge = -666;
    m_l1_pt = -666;
    m_l1_eta = -666;
    m_l1_phi = -666;
    m_l1_m = -666;

    m_l2_type = -666;
    m_l2_charge = -666;
    m_l2_pt = -666;
    m_l2_eta = -666;
    m_l2_phi = -666;
    m_l2_m = -666;

    m_l3_type = -666;
    m_l3_charge = -666;
    m_l3_pt = -666;
    m_l3_eta = -666;
    m_l3_phi = -666;
    m_l3_m = -666;

    m_l4_type = -666;
    m_l4_charge = -666;
    m_l4_pt = -666;
    m_l4_eta = -666;
    m_l4_phi = -666;
    m_l4_m = -666;

    m_j1_ntrk = -666;
    m_j1_btag = -666;
    m_j1_pt = -666;
    m_j1_eta = -666;
    m_j1_phi = -666;
    m_j1_m = -666;

    m_j2_ntrk = -666;
    m_j2_btag = -666;
    m_j2_pt = -666;
    m_j2_eta = -666;
    m_j2_phi = -666;
    m_j2_m = -666;

    m_j3_ntrk = -666;
    m_j3_btag = -666;
    m_j3_pt = -666;
    m_j3_eta = -666;
    m_j3_phi = -666;
    m_j3_m = -666;

    m_j4_ntrk = -666;
    m_j4_btag = -666;
    m_j4_pt = -666;
    m_j4_eta = -666;
    m_j4_phi = -666;
    m_j4_m = -666;

    m_j5_ntrk = -666;
    m_j5_btag = -666;
    m_j5_pt = -666;
    m_j5_eta = -666;
    m_j5_phi = -666;
    m_j5_m = -666;

    m_j6_ntrk = -666;
    m_j6_btag = -666;
    m_j6_pt = -666;
    m_j6_eta = -666;
    m_j6_phi = -666;
    m_j6_m = -666;

    m_met_phi = -666;
    m_met = -666;

    _v_ElectronTight->clear();
    _v_MuonTight->clear();
    _v_JetTight->clear();

    for(int i=0;i<_v_Electron->size();i++)
    {
        if( !_v_Electron->at(i).passPtEta() ) continue;
        if( _v_Electron->at(i).isLoose() )
            _v_ElectronTight->push_back(_v_Electron->at(i));
    }       

    for(int i=0;i<_v_Muon->size();i++)
    {
        //if( !_v_Muon->at(i).passPtEta() ) continue;
        if( _v_Muon->at(i).isLoose() )
            _v_MuonTight->push_back(_v_Muon->at(i));
    }      

    for(int i=0;i<_v_Jet->size();i++)
    {
        if( _v_Jet->at(i).isTight() && _v_Jet->at(i).pt() > 25. &&
                fabs(_v_Jet->at(i).eta()) < 2.4 )
            _v_JetTight->push_back(_v_Jet->at(i));
    }      

    std::vector<int> res = filterPt(0,0,0,
            _v_ElectronTight,
            _v_MuonTight);

    if( res[0] >= 0 )
    {
        int id = _v_Event->at(0).id();
        int run = _v_Event->at(0).run();
        int lumi = _v_Event->at(0).lumi();

        float metpt = _v_Event->at(0).metpt();
        float metphi = _v_Event->at(0).metphi();

        m_channel = _v_Event->at(0).tth_channel();

        m_met = metpt;
        m_met_phi = metphi;

        m_weight = _v_Event->at(0).mc_weight();

        int njets = _v_JetTight->size();
        int nbjets = 0;
        for(int ij=0;ij<njets;ij++)
        {
            if( _v_JetTight->at(ij).CSVv2() >= 0.244 ) nbjets++;
        }	

        int lidx1 = res[2];
        int lidx2 = res[3];
        int l1ise = res[4];
        int l2ise = res[5];

        int lidx3 = res[6];
        int lidx4 = res[7];
        int l3ise = res[8];
        int l4ise = res[9];

        int l1id = 0;
        int l1charge = 0;
        float l1pt = 0., l1eta = 0., l1phi = 0., l1m = 0.;
        int l2id = 0;
        int l2charge = 0;
        float l2pt = 0., l2eta = 0., l2phi = 0., l2m = 0.;
        bool l1_isLoose = 0;
        bool l2_isLoose = 0;
        bool l1_isTight = 0;
        bool l2_isTight = 0;
        bool l1_isTightMVA = 0;
        bool l2_isTightMVA = 0;
        bool l1_passCF = 0;
        bool l2_passCF = 0;

        int l3id = 0;
        int l3charge = 0;
        float l3pt = 0., l3eta = 0., l3phi = 0., l3m = 0.;

        int l4id = 0;
        int l4charge = 0;
        float l4pt = 0., l4eta = 0., l4phi = 0., l4m = 0.;

        if( l1ise == 1 ) 
        {
            l1id = _v_ElectronTight->at(lidx1).id();
            l1charge = _v_ElectronTight->at(lidx1).charge();
            l1pt = _v_ElectronTight->at(lidx1).pt();
            l1eta = _v_ElectronTight->at(lidx1).eta();
            l1phi = _v_ElectronTight->at(lidx1).phi();
            l1m = _v_ElectronTight->at(lidx1).m();
            l1_isLoose = _v_ElectronTight->at(lidx1).isLoose();
            l1_isTight = _v_ElectronTight->at(lidx1).isTight();
            l1_isTightMVA = _v_ElectronTight->at(lidx1).isTightMVA();
            //l1_passCF = _v_ElectronTight->at(lidx1).passChargeFlip();
            if( l2ise == 1 ) 
            {
                l2id = _v_ElectronTight->at(lidx2).id();
                l2charge = _v_ElectronTight->at(lidx2).charge();
                l2pt = _v_ElectronTight->at(lidx2).pt();
                l2eta = _v_ElectronTight->at(lidx2).eta();
                l2phi = _v_ElectronTight->at(lidx2).phi();
                l2m = _v_ElectronTight->at(lidx2).m();
                l2_isLoose = _v_ElectronTight->at(lidx2).isLoose();
                l2_isTight = _v_ElectronTight->at(lidx2).isTight();
                l2_isTightMVA = _v_ElectronTight->at(lidx2).isTightMVA();
                //l2_passCF = _v_ElectronTight->at(lidx2).passChargeFlip();
            }
            else
            {
                l2id = _v_MuonTight->at(lidx2).id();
                l2charge = _v_MuonTight->at(lidx2).charge();
                l2pt = _v_MuonTight->at(lidx2).pt();
                l2eta = _v_MuonTight->at(lidx2).eta();
                l2phi = _v_MuonTight->at(lidx2).phi();
                l2m = _v_MuonTight->at(lidx2).m();
                l2_isLoose = _v_MuonTight->at(lidx2).isLoose();
                l2_isTight = _v_MuonTight->at(lidx2).isTight();
                //l2_isTightMVA = _v_MuonTight->at(lidx2).isTightMVA();
                //l2_passCF = _v_MuonTight->at(lidx2).passChargeFlip();
            }	     
        }	
        else 
        {
            l1id = _v_MuonTight->at(lidx1).id();
            l1charge = _v_MuonTight->at(lidx1).charge();
            l1pt = _v_MuonTight->at(lidx1).pt();
            l1eta = _v_MuonTight->at(lidx1).eta();
            l1phi = _v_MuonTight->at(lidx1).phi();
            l1m = _v_MuonTight->at(lidx1).m();
            l1_isLoose = _v_MuonTight->at(lidx1).isLoose();
            l1_isTight = _v_MuonTight->at(lidx1).isTight();
            //l1_isTightMVA = _v_MuonTight->at(lidx1).isTightMVA();
            //l1_passCF = _v_MuonTight->at(lidx1).passChargeFlip();
            if( l2ise == 1 ) 
            {
                l2id = _v_ElectronTight->at(lidx2).id();
                l2charge = _v_ElectronTight->at(lidx2).charge();
                l2pt = _v_ElectronTight->at(lidx2).pt();
                l2eta = _v_ElectronTight->at(lidx2).eta();
                l2phi = _v_ElectronTight->at(lidx2).phi();
                l2m = _v_ElectronTight->at(lidx2).m();
                l2_isLoose = _v_ElectronTight->at(lidx2).isLoose();
                l2_isTight = _v_ElectronTight->at(lidx2).isTight();
                //l2_isTightMVA = _v_ElectronTight->at(lidx2).isTightMVA();
                //l2_passCF = _v_ElectronTight->at(lidx2).passChargeFlip();
            }
            else
            {
                l2id = _v_MuonTight->at(lidx2).id();
                l2charge = _v_MuonTight->at(lidx2).charge();
                l2pt = _v_MuonTight->at(lidx2).pt();
                l2eta = _v_MuonTight->at(lidx2).eta();
                l2phi = _v_MuonTight->at(lidx2).phi();
                l2m = _v_MuonTight->at(lidx2).m();
                l2_isLoose = _v_MuonTight->at(lidx2).isLoose();
                l2_isTight = _v_MuonTight->at(lidx2).isTight();
                //l2_isTightMVA = _v_MuonTight->at(lidx2).isTightMVA();
                //l2_passCF = _v_MuonTight->at(lidx2).passChargeFlip();
            }
        }

        if( lidx3 != -1 && l3ise == 1 )
        {
            l3id = _v_ElectronTight->at(lidx3).id();
            l3charge = _v_ElectronTight->at(lidx3).charge();
            l3pt = _v_ElectronTight->at(lidx3).pt();
            l3eta = _v_ElectronTight->at(lidx3).eta();
            l3phi = _v_ElectronTight->at(lidx3).phi();
            l3m = _v_ElectronTight->at(lidx3).m();
        }
        if( lidx3 != -1 && l3ise == 0 )
        {
            l3id = _v_MuonTight->at(lidx3).id();
            l3charge = _v_MuonTight->at(lidx3).charge();
            l3pt = _v_MuonTight->at(lidx3).pt();
            l3eta = _v_MuonTight->at(lidx3).eta();
            l3phi = _v_MuonTight->at(lidx3).phi();
            l3m = _v_MuonTight->at(lidx3).m();
        }	

        if( lidx4 != -1 && l4ise == 1 )
        {
            l4id = _v_ElectronTight->at(lidx4).id();
            l4charge = _v_ElectronTight->at(lidx4).charge();
            l4pt = _v_ElectronTight->at(lidx4).pt();
            l4eta = _v_ElectronTight->at(lidx4).eta();
            l4phi = _v_ElectronTight->at(lidx4).phi();
            l4m = _v_ElectronTight->at(lidx4).m();
        }
        if( lidx4 != -1 && l4ise == 0 )
        {
            l4id = _v_MuonTight->at(lidx4).id();
            l4charge = _v_MuonTight->at(lidx4).charge();
            l4pt = _v_MuonTight->at(lidx4).pt();
            l4eta = _v_MuonTight->at(lidx4).eta();
            l4phi = _v_MuonTight->at(lidx4).phi();
            l4m = _v_MuonTight->at(lidx4).m();
        }	

        // leptons
        if( l1id != 0 )
        {
            m_l1_type = (abs(l1id) == 11) ? 1 : 2;
            m_l1_charge = l1charge;
            m_l1_pt = l1pt;
            m_l1_eta = l1eta;
            m_l1_phi = l1phi;
            m_l1_m = l1m;
        }	
        if( l2id != 0 )
        {
            m_l2_type = (abs(l2id) == 11) ? 1 : 2;
            m_l2_charge = l2charge;
            m_l2_pt = l2pt;
            m_l2_eta = l2eta;
            m_l2_phi = l2phi;
            m_l2_m = l2m;
        }
        if( l3id != 0 )
        {
            m_l3_type = (abs(l3id) == 11) ? 1 : 2;
            m_l3_charge = l3charge;
            m_l3_pt = l3pt;
            m_l3_eta = l3eta;
            m_l3_phi = l3phi;
            m_l3_m = l3m;
        }
        if( l4id != 0 )
        {
            m_l4_type = (abs(l4id) == 11) ? 1 : 2;
            m_l4_charge = l4charge;
            m_l4_pt = l4pt;
            m_l4_eta = l4eta;
            m_l4_phi = l4phi;
            m_l4_m = l4m;
        }

        // jets
        for(int ij=0;ij<njets;ij++)
        {
            if( ij == 0 )
            {		  
                m_j1_ntrk = _v_JetTight->at(ij).ntrk();
                m_j1_btag = _v_JetTight->at(ij).CSVv2();
                m_j1_pt = _v_JetTight->at(ij).pt();
                m_j1_eta = _v_JetTight->at(ij).eta();
                m_j1_phi = _v_JetTight->at(ij).phi();
                m_j1_m = _v_JetTight->at(ij).m();
            }	     
            if( ij == 1 )
            {		  
                m_j2_ntrk = _v_JetTight->at(ij).ntrk();
                m_j2_btag = _v_JetTight->at(ij).CSVv2();
                m_j2_pt = _v_JetTight->at(ij).pt();
                m_j2_eta = _v_JetTight->at(ij).eta();
                m_j2_phi = _v_JetTight->at(ij).phi();
                m_j2_m = _v_JetTight->at(ij).m();
            }	     
            if( ij == 2 )
            {		  
                m_j3_ntrk = _v_JetTight->at(ij).ntrk();
                m_j3_btag = _v_JetTight->at(ij).CSVv2();
                m_j3_pt = _v_JetTight->at(ij).pt();
                m_j3_eta = _v_JetTight->at(ij).eta();
                m_j3_phi = _v_JetTight->at(ij).phi();
                m_j3_m = _v_JetTight->at(ij).m();
            }	     
            if( ij == 3 )
            {		  
                m_j4_ntrk = _v_JetTight->at(ij).ntrk();
                m_j4_btag = _v_JetTight->at(ij).CSVv2();
                m_j4_pt = _v_JetTight->at(ij).pt();
                m_j4_eta = _v_JetTight->at(ij).eta();
                m_j4_phi = _v_JetTight->at(ij).phi();
                m_j4_m = _v_JetTight->at(ij).m();
            }	     
            if( ij == 4 )
            {		  
                m_j5_ntrk = _v_JetTight->at(ij).ntrk();
                m_j5_btag = _v_JetTight->at(ij).CSVv2();
                m_j5_pt = _v_JetTight->at(ij).pt();
                m_j5_eta = _v_JetTight->at(ij).eta();
                m_j5_phi = _v_JetTight->at(ij).phi();
                m_j5_m = _v_JetTight->at(ij).m();
            }	     
            if( ij == 5 )
            {		  
                m_j6_ntrk = _v_JetTight->at(ij).ntrk();
                m_j6_btag = _v_JetTight->at(ij).CSVv2();
                m_j6_pt = _v_JetTight->at(ij).pt();
                m_j6_eta = _v_JetTight->at(ij).eta();
                m_j6_phi = _v_JetTight->at(ij).phi();
                m_j6_m = _v_JetTight->at(ij).m();
            }	     
        }       

        passSel = 0x0;
        passSel |= 1   << 0; 
        //	passSel |= (l1_isLoose && l2_isLoose) << 0;
        bool isEESS = (res[1] == 0 && l1ise && l2ise && l1_passCF && l2_passCF);
        bool isMMSS = (res[1] == 0 && !l1ise && !l2ise && l1_passCF && l2_passCF);
        passSel |= isEESS   << 1;
        passSel |= isMMSS   << 2;
        bool pass_lepPt = (l1pt >= 20. && l2pt >= 20.);
        passSel |= pass_lepPt   << 3;
        bool pass_isTight = (l1_isTight && l2_isTight);
        passSel |= pass_isTight   << 4;
        bool pass_tightMVA = (l1_isTightMVA && l2_isTightMVA);
        passSel |= pass_tightMVA   << 5;
        bool pass_nbjets = (nbjets >= 2);
        passSel |= pass_nbjets   << 6;
        bool pass_njets = (njets >= 4);
        passSel |= pass_njets   << 7;

        if( CHECK_BIT(passSel,0) ) m_sel = 0; // >=2lep
        if( CHECK_BIT(passSel,1) )
        {	     
            m_sel = 10; // eeSS
            if( CHECK_BIT(passSel,3) )
            {
                m_sel = 11; // eeSS_pt20
                if( CHECK_BIT(passSel,6) )
                {
                    m_sel = 12; // eeSS_pt20_bj2
                    if( CHECK_BIT(passSel,7) )
                        m_sel = 13; // eeSS_pt20_bj2_j4
                }		  
            }
        }
        if( CHECK_BIT(passSel,2) )
        {	     
            m_sel = 20; // mmSS
            if( CHECK_BIT(passSel,3) )
            {
                m_sel = 21; // mmSS_pt20
                if( CHECK_BIT(passSel,6) )
                {
                    m_sel = 22; // mmSS_pt20_bj2
                    if( CHECK_BIT(passSel,7) )
                        m_sel = 23; // mmSS_pt20_bj2_j4
                }		  
            }
        }

        if( CHECK_BIT(passSel,0) &&
                CHECK_BIT(passSel,2) &&
                CHECK_BIT(passSel,3) &&
                //	    CHECK_BIT(passSel,4) &&
                CHECK_BIT(passSel,5)
                ////	    CHECK_BIT(passSel,5) &&
                //	    CHECK_BIT(passSel,6) &&
                //	    CHECK_BIT(passSel,7)
          )
        {	     
            fprintf(_fevc,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %+2d  %6.2f %+4.2f %+4.2f    %6.1f  %+4.2f    %d \n",
                    run, lumi, id,
                    l1id, l1pt, l1eta, l1phi,
                    l2id, l2pt, l2eta, l2phi,
                    metpt, metphi,
                    njets);

            _trout->Fill();
        }	
    }   

    return 1;
}

std::vector<int> Hist::filterPt(double Pt1, double Pt2, double Pt3,
        std::vector<Electron>* ntElectron,
        std::vector<Muon>* ntMuon)
{
    int llChan = -1;
    int llChar = -1;

    bool foundPt1 = false;

    std::vector<std::pair<int,double> > el_pt;
    std::vector<std::pair<int,double> > mu_pt;
    std::vector<std::pair<int,double> > ll_pt;

    // descending
    std::sort(ntElectron->begin(),ntElectron->end(),Electron::sortPtPredicate);
    std::sort(ntMuon->begin(),ntMuon->end(),Muon::sortPtPredicate);

    for(int i=0;i<ntElectron->size();i++)
    {	
        el_pt.push_back(std::make_pair(i,ntElectron->at(i).pt()));
    }   
    for(int i=0;i<ntMuon->size();i++)
    {
        mu_pt.push_back(std::make_pair(i+1000,ntMuon->at(i).pt()));
    }

    ll_pt.reserve( el_pt.size() + mu_pt.size() );
    ll_pt.insert( ll_pt.end(), el_pt.begin(), el_pt.end() );
    ll_pt.insert( ll_pt.end(), mu_pt.begin(), mu_pt.end() );

    // descending
    std::sort( ll_pt.begin(), ll_pt.end(),
            boost::bind(&std::pair<int, double>::second, _1) >
            boost::bind(&std::pair<int, double>::second, _2));

    bool found_pt1 = false;
    bool found_pt2 = false;

    std::vector<int> id_remove;

    for(int i=0;i<ll_pt.size();i++)
    {
        if( found_pt1 && found_pt2 )
        {
            if( ll_pt.at(i).second < Pt3 )
                id_remove.push_back(i);
        }	

        if( ! found_pt2 )
        {	     
            if( ll_pt.at(i).second > Pt1 )
                found_pt1 = true;

            if( (found_pt1 && ll_pt.at(i).second < Pt2) || !found_pt1 )
                id_remove.push_back(i);
            else
                found_pt2 = true;
        }	
    }   

    for(int j=id_remove.size()-1;j>=0;j--)
    {
        int index = ll_pt.at(id_remove.at(j)).first;

        bool is_electron = true;

        if( index >= 1000 ) 
        {
            index = index - 1000;
            is_electron = false;
        }	     	

        ll_pt.erase(ll_pt.begin() + id_remove.at(j));

        if( is_electron )
            ntElectron->erase(ntElectron->begin() + index);
        else
            ntMuon->erase(ntMuon->begin() + index);
    }   

    bool l1_e = false;
    bool l1_m = false;

    bool l2_e = false;
    bool l2_m = false;

    bool l3_e = false;
    bool l3_m = false;

    bool l4_e = false;
    bool l4_m = false;

    int idx1 = -1;
    int idx2 = -1;
    int idx3 = -1;
    int idx4 = -1;

    // descending
    std::sort( ll_pt.begin(), ll_pt.end(),
            boost::bind(&std::pair<int, double>::second, _1) >
            boost::bind(&std::pair<int, double>::second, _2));

    if( ll_pt.size() > 1 )
    {
        if( ll_pt.at(0).first < 1000 ) l1_e = true;
        if( ll_pt.at(0).first >= 1000 ) l1_m = true;

        if( ll_pt.at(1).first < 1000 ) l2_e = true;
        if( ll_pt.at(1).first >= 1000 ) l2_m = true;

        if( ll_pt.size() > 2 )
        {
            if( ll_pt.at(2).first < 1000 ) l3_e = true;
            if( ll_pt.at(2).first >= 1000 ) l3_m = true;	     
        }	

        if( ll_pt.size() > 3 )
        {
            if( ll_pt.at(3).first < 1000 ) l4_e = true;
            if( ll_pt.at(3).first >= 1000 ) l4_m = true;	     
        }	

        // -1 - < 2 leptons
        //  0 - ee
        //  1 - mm
        //  2 - em

        // -1 - < 2 leptons
        //  0 - SS
        //  1 - OS

        if( l1_e && l2_e ) llChan = 0;
        if( l1_m && l2_m ) llChan = 1;
        if( (l1_e && l2_m) || 
                (l1_m && l2_e) ) llChan = 2;

        idx1 = (ll_pt.at(0).first < 1000) ? ll_pt.at(0).first : ll_pt.at(0).first - 1000;
        idx2 = (ll_pt.at(1).first < 1000) ? ll_pt.at(1).first : ll_pt.at(1).first - 1000;
        if( ll_pt.size() > 2 ) idx3 = (ll_pt.at(2).first < 1000) ? ll_pt.at(2).first : ll_pt.at(2).first - 1000;
        if( ll_pt.size() > 3 ) idx4 = (ll_pt.at(3).first < 1000) ? ll_pt.at(3).first : ll_pt.at(3).first - 1000;

        int cha1 = -1;
        int cha2 = -1;

        if( idx1 >= 0 )
        {	     
            if( l1_e ) cha1 = ntElectron->at(idx1).charge();
            if( l1_m ) cha1 = ntMuon->at(idx1).charge();
        }
        if( idx2 >= 0 )
        {	     
            if( l2_e ) cha2 = ntElectron->at(idx2).charge();
            if( l2_m ) cha2 = ntMuon->at(idx2).charge();
        }	

        if( idx1 >= 0 && idx2 >= 0 )
        {
            llChar = (cha1 * cha2 == 1) ? 0 : 1;
        }
    }

    std::vector<int> res;
    res.clear();
    res.push_back(llChan);
    res.push_back(llChar);
    res.push_back(idx1);
    res.push_back(idx2);
    res.push_back(l1_e);
    res.push_back(l2_e);

    res.push_back(idx3);
    res.push_back(idx4);
    res.push_back(l3_e);
    res.push_back(l4_e);

    return res;
}
