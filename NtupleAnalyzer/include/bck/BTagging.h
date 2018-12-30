#ifndef BTAGGING_H
#define BTAGGING_H


void fillCSVhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( std::vector<double> jetPts,
            std::vector<double> jetEtas,
            std::vector<double> jetCSVs,
            std::vector<int> jetFlavors,
            int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

TH1D* h_csv_wgt_hf[9][5];
TH1D* c_csv_wgt_hf[9][5];
TH1D* h_csv_wgt_lf[9][4][3];
TH2D* h_eff_btagging_b;
TH2D* h_eff_btagging_c;
TH2D* h_eff_btagging_l;

#endif
