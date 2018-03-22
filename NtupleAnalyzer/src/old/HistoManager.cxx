#include "../include/HistoManager.h"


//------------------------------------------------------
//initialize the historams for the analysis
//------------------------------------------------------

HistoManager::HistoManager()
{
    numb_histo1D = 0;
    numb_histo2D = 0;
}

HistoManager::~HistoManager()
{

}


//------------------------------------------------------
//add an histogram to the collection of histograms
//------------------------------------------------------
void HistoManager::addHisto(TString var, TString selstep, TString channel, TString sample, int nbins, float min, float max)
{
    TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
    TH1F * thehisto = new TH1F(histoname,histoname,nbins,min,max);
    thehisto->Sumw2();
    histo1D_list.push_back(thehisto);
    histo1D_map[histoname.Data()] = numb_histo1D;
    numb_histo1D++;
}


void HistoManager::addHisto2D( TString var, TString selstep, TString channel, TString sample, int nbins1, float min1, float max1, int nbins2, float min2, float max2)
{
    TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
    TH2F * thehisto = new TH2F(histoname,histoname,  nbins1, min1, max1, nbins2, min2, max2  );
    thehisto->Sumw2();
    histo2D_list.push_back(thehisto);
    histo2D_map[histoname.Data()] = numb_histo2D;
    numb_histo2D++;
}

// in combine datacard format
void HistoManager::addHisto(TString process, TString syst, int nbins, float min, float max)
{  
    TString histoname = process+"_"+syst;
    TH1F * thehisto = new TH1F(histoname,histoname,nbins,min,max);
    thehisto->Sumw2();
    histo1D_list.push_back(thehisto);
    histo1D_map[histoname.Data()] = numb_histo1D;
    numb_histo1D++;
}

//------------------------------------------------------
//add an histogram to the collection of histograms
//------------------------------------------------------
void HistoManager::fillHisto(TString var, TString selstep, TString channel, TString sample, float val, float weight)
{
    TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
    histo1D_list[histo1D_map[histoname.Data()]]->Fill(val, weight);
}


void HistoManager::fillHisto2D(TString var, TString selstep, TString channel, TString sample,float val1, float val2, float weight)
{
    TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
    histo2D_list[histo2D_map[histoname.Data()]]->Fill(val1, val2, weight);
}

// in combine datacard format
void HistoManager::fillHisto(TString process, TString syst, float val, float weight)
{
    TString histoname = process+"_"+syst;
    histo1D_list[histo1D_map[histoname.Data()]]->Fill(val, weight);
}



void HistoManager::writeHisto()
{
}


void HistoManager::addHisto1D(TH1F* h1)
{  
  histo1D_list.push_back(h1);
}


void HistoManager::addHisto2D(TH2F* h2)
{
  histo2D_list.push_back(h2);
}


TH1F* HistoManager::getHisto1D(TString var, TString selstep, TString channel, TString sample)
{  
  TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
  return histo1D_list[histo1D_map[histoname.Data()]];
}


TH2F* HistoManager::getHisto2D(TString var, TString selstep, TString channel, TString sample)
{
  TString histoname = var+"_"+selstep+"_"+channel+"__"+sample;
  return histo2D_list[histo2D_map[histoname.Data()]];
}


TH1F* HistoManager::getHisto1D(TString process, TString syst)
{  
  TString histoname = process+"_"+syst;
  return histo1D_list[histo1D_map[histoname.Data()]];
}
