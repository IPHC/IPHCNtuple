#ifndef PARTICLESELECTOR_H
#define PARTICLESELECTOR_H

#include "MultiLepton.h"
#include "MEPhaseSpace.h"

class ParticleSelector{

public:

  ParticleSelector();
  ~ParticleSelector();

  int nActualLeptons;
  int nActualJets;
  int nActualBjets;
 
  int nIntrinsicLeptons;  
  int nIntrinsicJets;
  int nIntrinsicBjets;

  bool doSelectOnlyBtaggedJets;
  double csvthresh;

  void SetIntrinsicNumOfParticles(MultiLepton*, int);

  void FillAllJets(MultiLepton*);

  int SelectJetsAndBjetsFromAllJet(MultiLepton*);
  void SelectBjets(MultiLepton*, string, int*, int*, bool);
  //void SelectLeptons(MultiLepton*);

  //void SetCatJets(MultiLepton*);

private:

};

ParticleSelector::ParticleSelector(){

  nActualLeptons = -1;
  nActualJets = -1;
  nActualBjets = -1;
  nIntrinsicLeptons = -1;
  nIntrinsicJets = -1;
  nIntrinsicBjets = -1;
  doSelectOnlyBtaggedJets = true;
  csvthresh = 0.5426;
}

ParticleSelector::~ParticleSelector(){

}

void ParticleSelector::FillAllJets(MultiLepton* multilepton_){

  (*multilepton_).AllJets.clear();

  for (unsigned int i=0; i<multilepton_->Bjets.size(); i++){
    for (unsigned int j=0; j<multilepton_->AllJets.size(); j++){
      if (multilepton_->Bjets.at(i).P4.Pt() != multilepton_->AllJets.at(j).P4.Pt())
        (*multilepton_).FillParticle("alljet", multilepton_->Bjets.at(i));
    }
  }

  for (unsigned int i=0; i<multilepton_->JetsHighestPt.size(); i++){
    for (unsigned int j=0; j<multilepton_->AllJets.size(); j++){
      if (multilepton_->JetsHighestPt.at(i).P4.Pt() != multilepton_->AllJets.at(j).P4.Pt())
        (*multilepton_).FillParticle("alljet", multilepton_->JetsHighestPt.at(i));
    }
  }

  for (unsigned int i=0; i<multilepton_->JetsClosestMw.size(); i++){
    for (unsigned int j=0; j<multilepton_->AllJets.size(); j++){
      if (multilepton_->JetsClosestMw.at(i).P4.Pt() != multilepton_->AllJets.at(j).P4.Pt())
        (*multilepton_).FillParticle("alljet", multilepton_->JetsClosestMw.at(i));
    }
  }

  for (unsigned int i=0; i<multilepton_->JetsLowestMjj.size(); i++){
    for (unsigned int j=0; j<multilepton_->AllJets.size(); j++){
      if (multilepton_->JetsLowestMjj.at(i).P4.Pt() != multilepton_->AllJets.at(j).P4.Pt())
        (*multilepton_).FillParticle("alljet", multilepton_->JetsLowestMjj.at(i));
    }
  }

  for (unsigned int i=0; i<multilepton_->JetsHighestEta.size(); i++){
    for (unsigned int j=0; j<multilepton_->AllJets.size(); j++){
      if (multilepton_->JetsHighestEta.at(i).P4.Pt() != multilepton_->AllJets.at(j).P4.Pt())
        (*multilepton_).FillParticle("alljet", multilepton_->JetsHighestEta.at(i));
    }
  }

  return;
}

int ParticleSelector::SelectJetsAndBjetsFromAllJet(MultiLepton* multilepton_){

  int res = -1;
  bool is2lss=false, is3l=false, is4l=false;

  if (multilepton_->AllJets.size()<2) return res;

  //B-jet cat
  int ib1=-1, ib2=-1;
  SelectBjets(multilepton_, "HighestBtagDiscrim", &ib1, &ib2, doSelectOnlyBtaggedJets); 
  (*multilepton_).Bjets.clear();
  if (ib1!=-1) (*multilepton_).FillParticle("bjet", multilepton_->AllJets.at(ib1));
  if (ib2!=-1) (*multilepton_).FillParticle("bjet", multilepton_->AllJets.at(ib2));

  //Lepton cat
  int tot_charge = 0;
  int tot_id = 0;
  if (multilepton_->Leptons.size()>=4) {
      for (unsigned int i=0; i<4; i++) {
          tot_charge += multilepton_->Leptons.at(i).Id/abs(multilepton_->Leptons.at(i).Id);
          tot_id += multilepton_->Leptons.at(i).Id;
      }
  }
  if (multilepton_->Leptons.size()>=4 && tot_charge==0 && (tot_id==0 || abs(tot_id)==2)) is4l = true;
  else if (multilepton_->Leptons.size()>=3) is3l = true;
  else if (multilepton_->Leptons.size()>=2 && multilepton_->Leptons.at(0).Id/abs(multilepton_->Leptons.at(0).Id)==-multilepton_->Leptons.at(1).Id/abs(multilepton_->Leptons.at(1).Id)) is2lss=true;
  else return res;

  res = 1;

  //2lss
  if (is2lss && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2>=4) (*multilepton_).kCatJets = kCat_2lss_2b_4j;
  else if (is2lss && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1>=4) (*multilepton_).kCatJets = kCat_2lss_1b_4j;
  else if (is2lss && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2==3) (*multilepton_).kCatJets = kCat_2lss_2b_3j;
  else if (is2lss && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1==3) (*multilepton_).kCatJets = kCat_2lss_1b_3j;
  else if (is2lss && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2==2) (*multilepton_).kCatJets = kCat_2lss_2b_2j;
  else if (is2lss && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2==1) (*multilepton_).kCatJets = kCat_2lss_2b_1j;
  else if (is2lss && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1==2) (*multilepton_).kCatJets = kCat_2lss_1b_2j;
  else if (is2lss && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1==1) (*multilepton_).kCatJets = kCat_2lss_1b_1j;
  //4l 
  else if (is4l && ib1!=-1 && ib2!=-1) (*multilepton_).kCatJets = kCat_4l_2b;
  else if (is4l && ib1!=-1 && ib2==-1) (*multilepton_).kCatJets = kCat_4l_1b;
  //3l
  else if (is3l && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2>=2) (*multilepton_).kCatJets = kCat_3l_2b_2j;
  else if (is3l && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1>=2) (*multilepton_).kCatJets = kCat_3l_1b_2j;
  else if (is3l && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2==1) (*multilepton_).kCatJets = kCat_3l_2b_1j;
  else if (is3l && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1==1) (*multilepton_).kCatJets = kCat_3l_1b_1j;
  else if (is3l && ib1!=-1 && ib2!=-1 && multilepton_->AllJets.size()-2==0) (*multilepton_).kCatJets = kCat_3l_2b_0j;
  else if (is3l && ib1!=-1 && ib2==-1 && multilepton_->AllJets.size()-1==0) (*multilepton_).kCatJets = kCat_3l_1b_0j;
  else if (is3l && ib1==-1 && ib2==-1 && multilepton_->AllJets.size()==1) (*multilepton_).kCatJets = kCat_3l_0b_1j;
  else if (is3l && ib1==-1 && ib2==-1 && multilepton_->AllJets.size()==0) (*multilepton_).kCatJets = kCat_3l_0b_0j;
  else  (*multilepton_).kCatJets = -1;

  //Choose 2 jets
    float pt_max=0, pt_max2=0, eta_max=0, eta_max2=0; int ij1=-1, ij2=-1;
    float diffmass_min = 10000, mass_min = 10000, mass=10000;
    int ik1=-1, ik2=-1, il1=-1, il2=-1, ie1=-1, ie2=-1;
    for (unsigned int ij=0; ij<multilepton_->AllJets.size(); ij++){
        if ((int)ij==ib1 || (int)ij==ib2) continue;
        if (multilepton_->AllJets.at(ij).P4.Pt() > pt_max ) {
            pt_max2 = pt_max;
            ij2 = ij1;
            pt_max = multilepton_->AllJets.at(ij).P4.Pt();
            ij1 = ij;
        } 
        if (multilepton_->AllJets.at(ij).P4.Pt() < pt_max && multilepton_->AllJets.at(ij).P4.Pt() > pt_max2){
            pt_max2 = multilepton_->AllJets.at(ij).P4.Pt(); 
            ij2 = ij; 
        }
	if (fabs(multilepton_->AllJets.at(ij).P4.Eta()) > fabs(eta_max)) {
	    eta_max2 = eta_max;
            ie2 = ie1;
            eta_max = multilepton_->AllJets.at(ij).P4.Eta();
            ie1 = ij;
	}
	if (fabs(multilepton_->AllJets.at(ij).P4.Eta()) < fabs(eta_max) && fabs(multilepton_->AllJets.at(ij).P4.Eta()) > fabs(eta_max2)){
	    eta_max2 = multilepton_->AllJets.at(ij).P4.Eta();
            ie2 = ij;
	}
        for (unsigned int ik=0; ik<multilepton_->AllJets.size(); ik++){
            if (ik==ij) continue;
            if ((int)ik==ib1 || (int)ik==ib2) continue;
	    mass = (multilepton_->AllJets.at(ij).P4+multilepton_->AllJets.at(ik).P4).M();
            if (TMath::Abs(mass-80.419)<diffmass_min){
                ik1=ij;
                ik2=ik;
                diffmass_min = TMath::Abs(mass-80.419);
            } 
            if (mass<mass_min){
                il1=ij;
                il2=ik;
                mass_min = mass;
            } 
        } 
    }  

 //Choose 2 more jets
    int io1=-1, io2=-1, ip1=-1, ip2=-1, im1=-1, im2=-1;
    diffmass_min = 10000, mass_min = 10000, pt_max2 = 0, pt_max = 0;
    for (unsigned int im=0; im<multilepton_->AllJets.size(); im++){
        if ((int)im==ib1 || (int)im==ib2 || (int)im==ik1 || (int)im==ik2) continue;
        if (multilepton_->AllJets.at(im).P4.Pt() > pt_max ) {
            pt_max2 = pt_max;
            im2 = im1;
            pt_max = multilepton_->AllJets.at(im).P4.Pt();
            im1 = im;
        }
        if (multilepton_->AllJets.at(im).P4.Pt() < pt_max && multilepton_->AllJets.at(im).P4.Pt() > pt_max2){
            pt_max2 = multilepton_->AllJets.at(im).P4.Pt();
            im2 = im;
        }
        for (unsigned int in=0; in<multilepton_->AllJets.size(); in++){
            if ((int)in==ib1 || (int)in==ib2 || (int)in==ik1 || (int)in==ik2 || in==im) continue;
            mass = (multilepton_->AllJets.at(im).P4+multilepton_->AllJets.at(in).P4).M();
            if (TMath::Abs(mass-80.419)<diffmass_min){
                io1=im;
                io2=in;
                diffmass_min = TMath::Abs(mass-80.419);
            }
            if (mass<mass_min){
                ip1=im;
                ip2=in;
                mass_min = mass;
            }
        }
    }

  (*multilepton_).JetsHighestPt.clear();
  (*multilepton_).JetsClosestMw.clear();
  (*multilepton_).JetsLowestMjj.clear();
  (*multilepton_).JetsHighestEta.clear(); 

  //First pair
  if (ij1!=-1 && ij2==-1) (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(ij1));
  if (ij1!=-1 && ij2!=-1) {
     (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(ij1));
     (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(ij2));
  }
  if (ik1!=-1 && ik2!=-1){
     (*multilepton_).FillParticle("jetClosestMw", multilepton_->AllJets.at(ik1));
     (*multilepton_).FillParticle("jetClosestMw", multilepton_->AllJets.at(ik2));
  }
  if (il1!=-1 && il2!=-1){
     (*multilepton_).FillParticle("jetLowestMjj", multilepton_->AllJets.at(il1));
     (*multilepton_).FillParticle("jetLowestMjj", multilepton_->AllJets.at(il2));
  }
  if (ie1!=-1 && ie2!=-1){
     (*multilepton_).FillParticle("jetHighestEta", multilepton_->AllJets.at(ie1));
     (*multilepton_).FillParticle("jetHighestEta", multilepton_->AllJets.at(ie2));
  }

  //Second pair (first one: closest to Mw)
  if (is2lss && ij1!=-1 && ij2!=-1){
    if (im1!=-1 && im2==-1){
      (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(im1));
    }
    if (im1!=-1 && im2!=-1){
      (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(im1));
      (*multilepton_).FillParticle("jetHighestPt", multilepton_->AllJets.at(im2));
    }
    if (io1!=-1 && io2!=-1){
      (*multilepton_).FillParticle("jetClosestMw", multilepton_->AllJets.at(io1));
      (*multilepton_).FillParticle("jetClosestMw", multilepton_->AllJets.at(io2));
    }
    if (ip1!=-1 && ip2!=-1){
      (*multilepton_).FillParticle("jetLowestMjj", multilepton_->AllJets.at(ip1));
      (*multilepton_).FillParticle("jetLowestMjj", multilepton_->AllJets.at(ip2));
    }
  }

 return res;
}

void ParticleSelector::SelectBjets(MultiLepton* multilepton_, std::string BjetSel, int* ibsel1, int* ibsel2, bool doSelectOnlyBjets){

  int ib1=-1, ib2=-1;

  if (BjetSel=="HighestBtagDiscrim"){
    float btag_max=-999, btag_max2=-999;
    for (unsigned int ib=0; ib<multilepton_->AllJets.size(); ib++){
      if (doSelectOnlyBjets && multilepton_->AllJets.at(ib).CSV < csvthresh) continue;
      if (multilepton_->AllJets.at(ib).CSV > btag_max){
        btag_max2 = btag_max;
        ib2 = ib1;
        btag_max = multilepton_->AllJets.at(ib).CSV;
        ib1 = ib;
      }
      if (multilepton_->AllJets.at(ib).CSV < btag_max && multilepton_->AllJets.at(ib).CSV > btag_max2){
	  btag_max2 = multilepton_->AllJets.at(ib).CSV;
	  ib2 = ib;
      }
    }
  }

  *ibsel1 = ib1;
  *ibsel2 = ib2;

}

void ParticleSelector::SetIntrinsicNumOfParticles(MultiLepton* multilepton_, int hyp){

  //To be  called once Jets collection is filled (in Permutations.h)

  nActualLeptons = multilepton_->Leptons.size(); 
  //nActualJets = multilepton_->AllJets.size(); //or a maximum value of 4 ?
  nActualJets = multilepton_->Jets.size();
  nActualBjets = multilepton_->Bjets.size(); //or a maximum value of 2 ?

  if (hyp==kMEM_TTH_TopAntitopHiggsDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 2;
    }
    if (nActualLeptons==4) {
      nIntrinsicLeptons=4;
      nIntrinsicJets = 0;
    }
  }

  if (hyp==kMEM_TTH_TopAntitopHiggsSemiLepDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==2) {
      nIntrinsicLeptons = 2;
      nIntrinsicJets = 4;
    }
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 2;
    }
  }

  if (hyp==kMEM_TTLL_TopAntitopDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 2;
    }
    if (nActualLeptons==4) {
      nIntrinsicLeptons=4;
      nIntrinsicJets = 0;
    }
  }

  if (hyp==kMEM_TTW_TopAntitopDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==2) {
      nIntrinsicLeptons = 2;
      nIntrinsicJets = 2;
    }
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 0;
    }
  }


  if (hyp==kMEM_TTWJJ_TopAntitopDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==2) {
      nIntrinsicLeptons = 2;
      nIntrinsicJets = 4;
    }
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 2;
    }
  }

  if (hyp==kMEM_TTbar_TopAntitopFullyLepDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 2;
      nIntrinsicJets = 0;
    }
  }

  if (hyp==kMEM_TTbar_TopAntitopSemiLepDecay){
    nIntrinsicBjets = 2;
    if (nActualLeptons==2) {
      nIntrinsicLeptons = 1;
      nIntrinsicJets = 2;
    }
    if (nActualLeptons==3){
      nIntrinsicLeptons = 1;
      nIntrinsicJets = 2;
    }
  }

  if (hyp==kMEM_TLLJ_TopLepDecay){
    nIntrinsicBjets = 1;
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 1;
    }
  }

  if (hyp==kMEM_WZJJ_LepDecay){
    nIntrinsicBjets = 0;
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 2;
    }
  }

  if (hyp==kMEM_THJ_TopLepDecay){
    nIntrinsicBjets = 1;
    if (nActualLeptons==2) {
      nIntrinsicLeptons = 2;
      nIntrinsicJets = 3;
    }
    if (nActualLeptons==3) {
      nIntrinsicLeptons = 3;
      nIntrinsicJets = 1;
    }
  }
}




#endif
