// NB : using "max(1, min(Nbins, FindBin) )" ensures that we always read a bin between 1 and N (included)

#include "../include/FakeRate.h"

using namespace std;

//--------------------------------------------
// ########  ########    ###    ########     ##     ## ####  ######  ########  #######   ######
// ##     ## ##         ## ##   ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##
// ##     ## ##        ##   ##  ##     ##    ##     ##  ##  ##          ##    ##     ## ##
// ########  ######   ##     ## ##     ##    #########  ##   ######     ##    ##     ##  ######
// ##   ##   ##       ######### ##     ##    ##     ##  ##        ##    ##    ##     ##       ##
// ##    ##  ##       ##     ## ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##
// ##     ## ######## ##     ## ########     ##     ## ####  ######     ##     #######   ######
//--------------------------------------------

//Get the histograms from files
void Fill_FR_Histograms(TFile* fileFR)
{
    //--- 2016
    //h_FR_el = (TH2D*)fileFR->Get("FR_mva075_el_data_comb");
    //h_FR_mu = (TH2D*)fileFR->Get("FR_mva075_mu_data_comb");

    //--- 2017

    //Nominal
    h_FR_nominal_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC");
    h_FR_nominal_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb");

    //Norm
    h_FR_normUp_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_up");
    h_FR_normDown_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_down");
    h_FR_normUp_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_up");
    h_FR_normDown_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_down");

    //pt
    h_FR_pt1_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_pt1");
    h_FR_pt2_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_pt2");

    h_FR_pt1_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_pt1");
    h_FR_pt2_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_pt2");

    //be (?)
    h_FR_be1_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_be1");
    h_FR_be2_el = (TH2D*)fileFR->Get("FR_mva090_el_data_comb_NC_be2");
    h_FR_be1_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_be1");
    h_FR_be2_mu = (TH2D*)fileFR->Get("FR_mva090_mu_data_comb_be2");

    return;
}




//--------------------------------------------
//  ######     ###    ##       ##          ######## ##     ## ##    ##  ######
// ##    ##   ## ##   ##       ##          ##       ##     ## ###   ## ##    ##
// ##        ##   ##  ##       ##          ##       ##     ## ####  ## ##
// ##       ##     ## ##       ##          ######   ##     ## ## ## ## ##
// ##       ######### ##       ##          ##       ##     ## ##  #### ##
// ##    ## ##     ## ##       ##          ##       ##     ## ##   ### ##    ##
//  ######  ##     ## ######## ########    ##        #######  ##    ##  ######
//--------------------------------------------

/**
 * Call the function computing the FR, and pass it the relevant histograms
 */
float Get_FR_Weight(std::vector<double> leptonsPts, std::vector<double> leptonsEtas, std::vector<int> leptonsIds, TString FR_type, bool isChannel_withoutEl=false, bool isChannel_withoutMu=false)
{
    //Pointers to the member histograms -- depends on FR_type, pass it to function
    TH2D* h_FR_el = 0;
    TH2D* h_FR_mu = 0;

    float FR_weight; //FR weight to retun

    if(FR_type == "" || FR_type == "nominal")
    {
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_norm_elUp")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_normUp_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_norm_elDown")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_normDown_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_norm_muUp")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_normUp_mu;
    }
    else if(FR_type == "FR_norm_muDown")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_normDown_mu;
    }
    else if(FR_type == "FR_pt_elUp")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_pt1_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_pt_elDown")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_pt2_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_pt_muUp")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_pt1_mu;
    }
    else if(FR_type == "FR_pt_muDown")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_pt2_mu;
    }
    else if(FR_type == "FR_be_elUp")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_be1_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_be_elDown")
    {
        if(isChannel_withoutEl) {return 0;} //set to nominal
        h_FR_el = h_FR_be2_el;
        h_FR_mu = h_FR_nominal_mu;
    }
    else if(FR_type == "FR_be_muUp")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_be1_mu;
    }
    else if(FR_type == "FR_be_muDown")
    {
        if(isChannel_withoutMu) {return 0;} //set to nominal
        h_FR_el = h_FR_nominal_el;
        h_FR_mu = h_FR_be2_mu;
    }
    else
    {
        cout<<"[FakeRate.cxx, l."<<__LINE__<<"] Error: wrong FR_type value !"<<endl;
        return 0;
    }

    FR_weight = Compute_FakeRate_Weight(leptonsPts, leptonsEtas, leptonsIds, h_FR_el, h_FR_mu);

    return FR_weight;
}







//--------------------------------------------
//  ######   ######## ########    ######## ########
// ##    ##  ##          ##       ##       ##     ##
// ##        ##          ##       ##       ##     ##
// ##   #### ######      ##       ######   ########
// ##    ##  ##          ##       ##       ##   ##
// ##    ##  ##          ##       ##       ##    ##
//  ######   ########    ##       ##       ##     ##
//--------------------------------------------


//1 fake -->  + f / (1-f)
//2 fakes --> - f1*f2 / ((1-f1) * (1-f2))
//3 fakes --> + f1*f2*f3 / ((1-f1) * (1-f2) * (1-f3))

/**
 * Compute the FR for the event
 */
float Compute_FakeRate_Weight( std::vector<double> leptonsPts, std::vector<double> leptonsEtas, std::vector<int> leptonsIds, TH2D* h_FR_el, TH2D* h_FR_mu)
{
    double weight    = -1;

    if(!leptonsPts.size() ) {std::cout<<"Problem in FakeRate.cxx : all leptons are tight !"<<std::endl; return 0;}

    for(int i=0; i<leptonsPts.size(); i++)
    {
        double leptonPt   = leptonsPts[i];
        double leptonEta  = fabs( leptonsEtas[i] );
        int    leptonId   = abs(  leptonsIds[i]   );

	    double f = 0;

        if(leptonId == 11)
        {
	    //Find corresponding bin
	    int x = std::max(1, std::min(h_FR_el->GetNbinsX(), h_FR_el->GetXaxis()->FindBin(leptonPt) ) );
	    int y = std::max(1, std::min(h_FR_el->GetNbinsY(), h_FR_el->GetYaxis()->FindBin(leptonEta) ) );

            f = h_FR_el->GetBinContent(x,y);
        }
        else if(leptonId == 13)
        {
	    //Find corresponding bin
	    int x = std::max(1, std::min(h_FR_mu->GetNbinsX(), h_FR_mu->GetXaxis()->FindBin(leptonPt) ) );
	    int y = std::max(1, std::min(h_FR_mu->GetNbinsY(), h_FR_mu->GetYaxis()->FindBin(leptonEta) ) );

            f = h_FR_mu->GetBinContent(x,y);
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
