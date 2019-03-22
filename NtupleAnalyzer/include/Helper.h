#ifndef HELPER_H
#define HELPER_H

bool SortingJetPt(Jet j1, Jet j2);
bool SortingLeptonPt(Lepton l1, Lepton l2);
bool SortingLeptonConePt(Lepton l1, Lepton l2);
float DeltaRLeptonJet(Lepton l1, Jet j1);
float DeltaRJets(Jet j1, Jet j2);

float fwdjet_eventWeight_2017_option3_modified(float eta);

#endif
