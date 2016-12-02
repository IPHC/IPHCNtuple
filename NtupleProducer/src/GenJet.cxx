#include "include/NtupleProducer.h"

ClassImp(GenJet)

GenJet::GenJet()
{
}

GenJet::~GenJet()
{
}

void GenJet::read()
{
    _ID = idx;

    _genJet_pt = ntP->genJet_pt->at(idx);
    _genJet_eta = ntP->genJet_eta->at(idx);
    _genJet_phi = ntP->genJet_phi->at(idx);
    _genJet_E = ntP->genJet_E->at(idx);
          
}


void GenJet::init()
{

    _genJet_pt  = -100;
    _genJet_eta = -100;
    _genJet_phi = -100;
    _genJet_E   = -100;
    			 
}


bool GenJet::sel()
{
    
    bool isSel = false;
    
    if (_genJet_pt > 25 && fabs(_genJet_eta) < 2.5) isSel = true;
    
    return isSel;
    			 
}
