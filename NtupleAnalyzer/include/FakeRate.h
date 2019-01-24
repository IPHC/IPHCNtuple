#ifndef FAKERATE_H
#define FAKERATE_H

//--------------------------------------------

void Fill_FR_Histograms(TFile *fileFR);
float Get_FR_Weight(std::vector<double>, std::vector<double>, std::vector<int>, TString, bool, bool);
float Compute_FakeRate_Weight(std::vector<double>, std::vector<double>, std::vector<int>, TH2D*, TH2D*);

//--------------------------------------------

//Histograms
// TH2D* h_FR_wgt_el;
// TH2D* h_FR_wgt_mu;

TH2D* h_FR_nominal_el = 0;
TH2D* h_FR_nominal_mu = 0;
TH2D* h_FR_normUp_el = 0;
TH2D* h_FR_normDown_el = 0;
TH2D* h_FR_normUp_mu = 0;
TH2D* h_FR_normDown_mu = 0;
TH2D* h_FR_pt1_el = 0;
TH2D* h_FR_pt2_el = 0;
TH2D* h_FR_pt1_mu = 0;
TH2D* h_FR_pt2_mu = 0;
TH2D* h_FR_be1_el = 0;
TH2D* h_FR_be2_el = 0;
TH2D* h_FR_be1_mu = 0;
TH2D* h_FR_be2_mu = 0;

//--------------------------------------------

#endif
