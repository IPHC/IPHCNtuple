// NB : using "max(1, min(Nbins, FindBin) )" ensures that we always read a bin between 1 and N (included)

#include "../include/FakeRate.h"

// fill the histograms (done once)
void fillFRhistos(TFile* fileFR)
{
    //2016
    //h_FR_wgt_el = (TH2D*)fileFR->Get("FR_mva075_el_data_comb");
    //h_FR_wgt_mu = (TH2D*)fileFR->Get("FR_mva075_mu_data_comb");
    
    //2017
    h_FR_wgt_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC");
    h_FR_wgt_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb");

    return;
}


//1 fake -->  + f / (1-f)
//2 fakes --> - f1*f2 / ((1-f1) * (1-f2))
//3 fakes --> + f1*f2*f3 / ((1-f1) * (1-f2) * (1-f3))
double get_FR_weight( std::vector<double> leptonsPts, std::vector<double> leptonsEtas, std::vector<int> leptonsIds)
{
    double weight    = -1;
    
    if(!leptonsPts.size() ) {std::cout<<"Problem in FakeRate.cxx : all leptons are tight !"<<std::endl; return 0;}

    for(int i=0; i<leptonsPts.size(); i++)
    {   
        double leptonPt   = leptonsPts[i];
        double leptonEta  = fabs( leptonsEtas[i] ); //Use absolute eta
        int    leptonId   = abs(  leptonsIds[i]   );
	
	double f = 0;

        if(leptonId == 11)
        {
	    //Find corresponding bin
	    int x = std::max(1, std::min(h_FR_wgt_el->GetNbinsX(), h_FR_wgt_el->GetXaxis()->FindBin(leptonPt) ) );
	    int y = std::max(1, std::min(h_FR_wgt_el->GetNbinsY(), h_FR_wgt_el->GetYaxis()->FindBin(leptonEta) ) );
	
            f = h_FR_wgt_el->GetBinContent(x,y);
        }
        else if(leptonId == 13)
        {
	    //Find corresponding bin
	    int x = std::max(1, std::min(h_FR_wgt_mu->GetNbinsX(), h_FR_wgt_mu->GetXaxis()->FindBin(leptonPt) ) );
	    int y = std::max(1, std::min(h_FR_wgt_mu->GetNbinsY(), h_FR_wgt_mu->GetYaxis()->FindBin(leptonEta) ) );
	    
            f = h_FR_wgt_mu->GetBinContent(x,y);
        }
	else
	{
	    std::cout<<"Problem in FakeRate.cxx : wrong lepton flavour !"<<std::endl;
	    return 0;
	}
	
	if(f==0 || f ==1) {cout<<"Problem in FakeRate.cxx : factor = "<<f<<endl;}

        weight*= -f / (1-f);
    }
    
    return weight;
}




