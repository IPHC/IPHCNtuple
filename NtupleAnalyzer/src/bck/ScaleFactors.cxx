#include <assert.h>
#include "TH2F.h"
#include "TFile.h"

TFile *_file_recoToLoose_leptonSF_mu = NULL;
TH2F *_histo_recoToLoose_leptonSF_mu = NULL;
TFile *_file_recoToLoose_leptonSF_el = NULL;
TH2F *_histo_recoToLoose_leptonSF_el1 = NULL;
TH2F *_histo_recoToLoose_leptonSF_el2 = NULL;

float _get_recoToLoose_leptonSF_ttH(int pdgid, float pt, float eta, int nlep, float var){

    if (var!=0) assert(0); // NOT IMPLEMENTED

    if (!_histo_recoToLoose_leptonSF_mu) {
        _file_recoToLoose_leptonSF_mu = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/mu_eff_recoToLoose_ttH.root","read");
        _histo_recoToLoose_leptonSF_mu = (TH2F*)(_file_recoToLoose_leptonSF_mu->Get("FINAL"));
    }
    if (!_histo_recoToLoose_leptonSF_el1) {
        _file_recoToLoose_leptonSF_el = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/kinematicBinSFele.root","read");
        _histo_recoToLoose_leptonSF_el1 = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("MVAVLooseFO_and_IDEmu_and_TightIP2D"));
        _histo_recoToLoose_leptonSF_el2 = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("MiniIso0p4_vs_AbsEta"));
    }

    if (abs(pdgid)==13){
        TH2F *hist = _histo_recoToLoose_leptonSF_mu;
        int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(eta)));
        int ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pt)));
        return hist->GetBinContent(etabin,ptbin);
    }
    if (abs(pdgid)==11){
        TH2F *hist = _histo_recoToLoose_leptonSF_el1;
        int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
        int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
        float out = hist->GetBinContent(ptbin,etabin);
        hist = _histo_recoToLoose_leptonSF_el2;
        ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
        etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
        out *= hist->GetBinContent(ptbin,etabin);
        return out;
    }

    assert(0);
    return -999;

}

TFile *_file_looseToTight_leptonSF_mu_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_mu_2lss = NULL;
TFile *_file_looseToTight_leptonSF_el_2lss = NULL;
TH2F *_histo_looseToTight_leptonSF_el_2lss = NULL;
TFile *_file_looseToTight_leptonSF_mu_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_mu_3l = NULL;
TFile *_file_looseToTight_leptonSF_el_3l = NULL;
TH2F *_histo_looseToTight_leptonSF_el_3l = NULL;

float _get_looseToTight_leptonSF_ttH(int pdgid, float _pt, float eta, int nlep, float var){

    float pt = std::min(float(79.9),_pt);

    if (var!=0) assert(0); // NOT IMPLEMENTED

    if (!_histo_looseToTight_leptonSF_mu_2lss) {
        _file_looseToTight_leptonSF_mu_2lss = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/lepMVAEffSF_m_2lss.root","read");
        _histo_looseToTight_leptonSF_mu_2lss = (TH2F*)(_file_looseToTight_leptonSF_mu_2lss->Get("sf"));
    }
    if (!_histo_looseToTight_leptonSF_el_2lss) {
        _file_looseToTight_leptonSF_el_2lss = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/lepMVAEffSF_e_2lss.root","read");
        _histo_looseToTight_leptonSF_el_2lss = (TH2F*)(_file_looseToTight_leptonSF_el_2lss->Get("sf"));
    }
    if (!_histo_looseToTight_leptonSF_mu_3l) {
        _file_looseToTight_leptonSF_mu_3l = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/lepMVAEffSF_m_3l.root","read");
        _histo_looseToTight_leptonSF_mu_3l = (TH2F*)(_file_looseToTight_leptonSF_mu_3l->Get("sf"));
    }
    if (!_histo_looseToTight_leptonSF_el_3l) {
        _file_looseToTight_leptonSF_el_3l = new TFile("/afs/cern.ch/user/p/peruzzi/work/tthtrees/cms_utility_files/lepMVAEffSF_e_3l.root","read");
        _histo_looseToTight_leptonSF_el_3l = (TH2F*)(_file_looseToTight_leptonSF_el_3l->Get("sf"));
    }

    TH2F *hist = 0;
    if (abs(pdgid)==13) hist = (nlep>2) ? _histo_looseToTight_leptonSF_mu_3l : _histo_looseToTight_leptonSF_mu_2lss;
    else if (abs(pdgid)==11) hist = (nlep>2) ? _histo_looseToTight_leptonSF_el_3l : _histo_looseToTight_leptonSF_el_2lss;
    assert(hist);
    int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(fabs(eta))));
    return hist->GetBinContent(ptbin,etabin);

}

float leptonSF_ttH(int pdgid, float pt, float eta, int nlep, float var=0){

    float recoToLoose = _get_recoToLoose_leptonSF_ttH(pdgid,pt,eta,nlep,var);
    float looseToTight = _get_looseToTight_leptonSF_ttH(pdgid,pt,eta,nlep,var);
    float res = recoToLoose*looseToTight;
    assert (res>0);
    return res;

}

float triggerSF_ttH(int pdgid1, float pt1, int pdgid2, float pt2, int nlep, float var_ee=0){
    if (var_ee!=0) assert(0); // NOT IMPLEMENTED
    if (nlep>2) return 1;
    if (abs(pdgid1)==11 && abs(pdgid2)==11){
        if (std::max(pt1,pt2)<40) return 0.95;
        else return 0.99;
    }
    else if (abs(pdgid1)==13 && abs(pdgid2)==13) {
        return 1.;
    }
    else return 0.98;
}
