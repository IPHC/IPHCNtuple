#ifndef ANALYZER
#define ANALYZER

#include <vector>
#include <iostream>
#include <fstream>

#include "Electron.h"
#include "Muon.h"
#include "Jet.h"
#include "Event.h"
#include "Truth.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TRandom3.h"

#define CHECK_BIT(var,pos) (((var) & (1<<(pos))) != 0 )

#endif
