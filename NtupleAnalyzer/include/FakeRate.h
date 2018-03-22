#ifndef FAKERATE_H
#define FAKERATE_H


void fillFRhistos(TFile *fileFR);


double get_FR_weight(std::vector<double> leptonsPts, std::vector<double> leptonsEtas, std::vector<int> leptonsIds);

TH2D* h_FR_wgt_el;
TH2D* h_FR_wgt_mu;

#endif
