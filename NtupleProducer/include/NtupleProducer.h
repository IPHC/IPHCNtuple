#ifndef NTUPLEPRODUCER_H
#define NTUPLEPRODUCER_H

#include "Tree.h"
#include "EventExt.h"
#include "ElectronExt.h"
#include "MuonExt.h"
#include "TauExt.h"
#include "JetExt.h"
#include "TruthExt.h"
#include "GenJetExt.h"
#include "TriggerObjExt.h"
#include "Ntuple.h"
#include "Sync.h"


// JES
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "TRandom3.h"

extern Tree *ntP;
extern TChain *ch;
extern Ntuple *nt;
extern Sync *sc;

#endif
