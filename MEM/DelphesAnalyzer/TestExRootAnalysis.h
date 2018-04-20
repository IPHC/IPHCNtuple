#ifndef TESTEXROOTANALYSIS_H
#define TESTEXROOTANALYSIS_H

#include <assert.h>
#include "MultiLeptonTree.h"

class Lepton {
public:
  Lepton();
  Lepton(Electron* elec):PT(elec->PT),Eta(elec->Eta),Phi(elec->Phi),Charge(elec->Charge),Id((elec->Charge<0)?11:-11){};
  Lepton(Muon* muon):PT(muon->PT),Eta(muon->Eta),Phi(muon->Phi),Charge(muon->Charge),Id((muon->Charge<0)?13:-13){};
  ~Lepton();
  Float_t PT;
  Float_t Eta;
  Float_t Phi;
  Int_t Charge;
  Int_t Id;
};

bool SortingLeptonPt( Lepton* l1, Lepton* l2)
{
    //cout << "l1PT="<<l1->PT<<" l2PT="<<l2->PT<<endl;
    if( l1->PT > l2->PT ) return true;
    else return false;
}

void selectBjets(std::vector<Jet*>& vSelectedJets, std::string BjetSel, int* ibsel1, int* ibsel2, bool doSelectOnlyBjets){

    //Selects the two highest b-tag jets. If only one b-tag select just this one.
    int ib1=-1, ib2=-1;

    if (BjetSel=="HighestBtagDiscrim"){
        float btag_max=-999, btag_max2=-999;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            Jet *jet = vSelectedJets.at(ib);	  
            //if (doSelectOnlyBjets && (vSelectedJets.at(ib).CSVv2()<0.5426)) continue;
            if (doSelectOnlyBjets && jet->BTag==0) continue;
/*
            if (vSelectedJets.at(ib).CSVv2()>btag_max){
                btag_max2 = btag_max;
                ib2 = ib1;
                btag_max = vSelectedJets.at(ib).CSVv2();
                ib1 = ib;
            }
            if (vSelectedJets.at(ib).CSVv2()<btag_max && vSelectedJets.at(ib).CSVv2()>btag_max2){
                btag_max2 = vSelectedJets.at(ib).CSVv2();
                ib2 = ib;
            }
*/
        }
    }
    if (BjetSel=="BtagHighestPt"){
        float pt_max=0, pt_max2=0;
        for (unsigned int ib=0; ib<vSelectedJets.size(); ib++){
            Jet *jet = vSelectedJets.at(ib);
            if (doSelectOnlyBjets && jet->BTag==0) continue;
            if (jet->PT>pt_max){
                pt_max2 = pt_max;
                ib2 = ib1;
                pt_max = jet->PT;
                ib1 = ib;
            }
            if (jet->PT<pt_max && jet->PT>pt_max2){
                pt_max2 = jet->PT;
                ib2 = ib;
            }
        }
    }

    *ibsel1 = ib1;
    *ibsel2 = ib2;

}

int GetTransferFunctionEta(float Reco_Eta){

  int iEta = -1;
  for (int iBin=0; iBin<3; iBin++){
    if (fabs(Reco_Eta) > EtaRange[iBin] && fabs(Reco_Eta) < EtaRange[iBin+1]) iEta=iBin;
  }
  return iEta;
}

int GetTransferFunctionEnergy(float Reco_E){

  int iEnergy = -1;
  for (int iBin=0; iBin<6; iBin++){
    if (Reco_E > EnergyRange[iBin] && Reco_E < EnergyRange[iBin+1]) iEnergy=iBin;
  }
  return iEnergy;
}

int GetTransferFunctionMet(float Reco_MET, float Reco_SumET){

  int iMet = -1;
  for (int iBin=0; iBin<2; iBin++){ //MET
    if (fabs(Reco_MET) > MetRange[iBin] && fabs(Reco_MET)<MetRange[iBin+1]) iMet = iBin;
  }
  for (int iBin=0; iBin<3; iBin++){ //SumET
    if (Reco_SumET > MetSumRange[iBin] && Reco_SumET<MetSumRange[iBin+1]) iMet = iMet*3 + iBin;
  }

  return iMet;
}

#endif

