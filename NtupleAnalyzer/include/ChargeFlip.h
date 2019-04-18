#ifndef CHARGEFLIP_H
#define CHARGEFLIP_H


void fillQFhistos(TFile *fileQF);

double get_QF_weight( std::vector<double> leptonsPts, std::vector<double> leptonsEtas);

TH2D* h_QF_wgt;

#endif
