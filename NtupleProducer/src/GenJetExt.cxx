#include "include/NtupleProducer.h"

ClassImp(GenJetExt)

GenJetExt::GenJetExt()
{
}

GenJetExt::~GenJetExt()
{
}

void GenJetExt::read()
{
   ID = idx;
   
   genJet_pt = ntP->genJet_pt->at(idx);
   genJet_eta = ntP->genJet_eta->at(idx);
   genJet_phi = ntP->genJet_phi->at(idx);
   genJet_E = ntP->genJet_E->at(idx);   
}


void GenJetExt::init()
{
   genJet_pt  = -100;
   genJet_eta = -100;
   genJet_phi = -100;
   genJet_E   = -100;   
}


bool GenJetExt::sel()
{    
   bool isSel = false;
   
   if (genJet_pt > 25 && fabs(genJet_eta) < 2.5) isSel = true;
   
   return isSel;    			 
}
