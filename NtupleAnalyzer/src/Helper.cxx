#include "../include/Helper.h"
#include <assert.h>
#include "TH2F.h"
#include "TFile.h"

# define M_PI           3.14159265358979323846

bool SortingLeptonPt( Lepton l1, Lepton l2)
{
    if( l1.pt > l2.pt ) return true;
    else return false;
}

float DeltaRLeptonJet( Lepton l1, Jet j1)
{
    float dEta = l1.eta - j1.eta;
    float dPhi = l1.phi - j1.phi;
    while(dPhi >= M_PI) dPhi -= 2*M_PI;
    while(dPhi < -M_PI) dPhi += 2*M_PI;
    float deltaR = sqrt( dEta*dEta + dPhi*dPhi);
    return deltaR;
}

float DeltaRJets( Jet j1, Jet j2)
{
    float dEta = j1.eta - j2.eta;
    float dPhi = j1.phi - j2.phi;
    while(dPhi >= M_PI) dPhi -= 2*M_PI;
    while(dPhi < -M_PI) dPhi += 2*M_PI;
    float deltaR = sqrt( dEta*dEta + dPhi*dPhi);
    return deltaR;
}
