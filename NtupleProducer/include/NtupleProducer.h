#ifndef NTUPLEPRODUCER_H
#define NTUPLEPRODUCER_H

#include "Tree.h"
#include "Event.h"
#include "Electron.h"
#include "Muon.h"
#include "Tau.h"
#include "Jet.h"
#include "Truth.h"
#include "GenJet.h"
#include "TriggerObj.h"
#include "Ntuple.h"
// JES
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

extern Tree             *ntP;
extern TChain           *ch;
extern Ntuple           *nt;
extern std::vector<int> *evdebug;

#endif
