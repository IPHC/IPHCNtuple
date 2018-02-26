#include "tHq_analysis.h"

using namespace std;

int main(int argc, char **argv) //Can choose region (tZq/WZ/ttZ) at execution
{
    if(argc >= 2) {cout<<BOLD(FRED("Error -- no arguments"))<<endl; return 1;}

//---------------------------------



//---------------------------------------------------------------------------
//  #######  ########  ######## ####  #######  ##    ##  ######
// ##     ## ##     ##    ##     ##  ##     ## ###   ## ##    ##
// ##     ## ##     ##    ##     ##  ##     ## ####  ## ##
// ##     ## ########     ##     ##  ##     ## ## ## ##  ######
// ##     ## ##           ##     ##  ##     ## ##  ####       ##
// ##     ## ##           ##     ##  ##     ## ##   ### ##    ##
//  #######  ##           ##    ####  #######  ##    ##  ######
//---------------------------------------------------------------------------

    //Matrix Element Method ==> TRUE
    bool include_MEM_variables = false;

    double set_luminosity = 35.862; //Moriond 2017 - 35.862fb
    TString plot_extension = ".png";

    //Templates options
    int nofbin_templates = 10; //Templates binning ==> 10 bins

    TString template_name = "ttV";

//---------------------------------------------------------------------------
// ########  #### ########            ########     ###    ######## ##     ##
// ##     ##  ##  ##     ##           ##     ##   ## ##      ##    ##     ##
// ##     ##  ##  ##     ##           ##     ##  ##   ##     ##    ##     ##
// ##     ##  ##  ########            ########  ##     ##    ##    #########
// ##     ##  ##  ##   ##             ##        #########    ##    ##     ##
// ##     ##  ##  ##    ##  ###       ##        ##     ##    ##    ##     ##
// ########  #### ##     ## ###       ##        ##     ##    ##    ##     ##
//---------------------------------------------------------------------------
//Set here the path of the directory containing all the Ntuples.

    TString dir_ntuples; TString t_name;
    dir_ntuples="input_ntuples";
    t_name = "Tree";

//---------------------------------------------------------------------------
//  ######  ##     ##    ###    ##    ## ##    ## ######## ##        ######
// ##    ## ##     ##   ## ##   ###   ## ###   ## ##       ##       ##    ##
// ##       ##     ##  ##   ##  ####  ## ####  ## ##       ##       ##
// ##       ######### ##     ## ## ## ## ## ## ## ######   ##        ######
// ##       ##     ## ######### ##  #### ##  #### ##       ##             ##
// ##    ## ##     ## ##     ## ##   ### ##   ### ##       ##       ##    ##
//  ######  ##     ## ##     ## ##    ## ##    ## ######## ########  ######
//---------------------------------------------------------------------------
    std::vector<TString > thechannellist;

    thechannellist.push_back("");

//---------------------------------------------------------------------------
//  ######     ###    ##     ## ########  ##       ########  ######
// ##    ##   ## ##   ###   ### ##     ## ##       ##       ##    ##
// ##        ##   ##  #### #### ##     ## ##       ##       ##
//  ######  ##     ## ## ### ## ########  ##       ######    ######
//       ## ######### ##     ## ##        ##       ##             ##
// ##    ## ##     ## ##     ## ##        ##       ##       ##    ##
//  ######  ##     ## ##     ## ##        ######## ########  ######
//---------------------------------------------------------------------------
    std::vector<TString> thesamplelist;
    std::vector<int> v_color; //sample <-> color

//-------------------
    //DATA --- 1 sample in first position
    // thesamplelist.push_back("Data");

    //Signal --- must be placed before backgrounds
    thesamplelist.push_back("tHq");        v_color.push_back(kGreen+1);
    thesamplelist.push_back("tHW");       v_color.push_back(920); //grey

    //BKG
    thesamplelist.push_back("WZ");        v_color.push_back(kMagenta); //grey
    thesamplelist.push_back("ZZ");        v_color.push_back(kYellow+1);
    thesamplelist.push_back("ttZ");       v_color.push_back(kRed-4);
    thesamplelist.push_back("ttW");       v_color.push_back(kRed+3);//Keep ttW & ttH samples together (coloring)
    thesamplelist.push_back("ttH");       v_color.push_back(kRed+3);
    thesamplelist.push_back("tZq");       v_color.push_back(kOrange+1);

    //Fakes
    // thesamplelist.push_back("TTbar");       v_color.push_back(kOrange+1);
    // thesamplelist.push_back("DY");       v_color.push_back(kOrange+1);

    //-- Use custom color palette
    bool use_custom_colorPalette = false;
//-------------------


//---------------------------------------------------------------------------
// ########  ########  ########       ##     ##    ###    ########   ######
// ##     ## ##     ##    ##          ##     ##   ## ##   ##     ## ##    ##
// ##     ## ##     ##    ##          ##     ##  ##   ##  ##     ## ##
// ########  ##     ##    ##          ##     ## ##     ## ########   ######
// ##     ## ##     ##    ##           ##   ##  ######### ##   ##         ##
// ##     ## ##     ##    ##            ## ##   ##     ## ##    ##  ##    ##
// ########  ########     ##             ###    ##     ## ##     ##  ######
//---------------------------------------------------------------------------
    std::vector<TString > thevarlist; //Variables used in BDT

//------------------------ for tZq
    thevarlist.push_back("nJet25");
    thevarlist.push_back("maxEtaJet25");
    thevarlist.push_back("lepCharge");
    thevarlist.push_back("nJetEta1");
    thevarlist.push_back("dEtaFwdJetBJet");
    thevarlist.push_back("dEtaFwdJet2BJet");
    thevarlist.push_back("dEtaFwdJetClosestLep");
    thevarlist.push_back("dPhiHighestPtSSPair");
    thevarlist.push_back("minDRll");
    thevarlist.push_back("Lep3Pt");


//---------------------------------------------------------------------------
//  ######  ##    ##  ######  ######## ######## ##     ##    ###    ######## ####  ######   ######
// ##    ##  ##  ##  ##    ##    ##    ##       ###   ###   ## ##      ##     ##  ##    ## ##    ##
// ##         ####   ##          ##    ##       #### ####  ##   ##     ##     ##  ##       ##
//  ######     ##     ######     ##    ######   ## ### ## ##     ##    ##     ##  ##        ######
//       ##    ##          ##    ##    ##       ##     ## #########    ##     ##  ##             ##
// ##    ##    ##    ##    ##    ##    ##       ##     ## ##     ##    ##     ##  ##    ## ##    ##
//  ######     ##     ######     ##    ######## ##     ## ##     ##    ##    ####  ######   ######
//---------------------------------------------------------------------------

    bool use_systematics = true;

    vector<TString> thesystlist;
    thesystlist.push_back("");

//---------------------------------------------------------------------------
// ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##        ######     ###    ##       ##        ######
// ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ##       ##    ##   ## ##   ##       ##       ##    ##
// ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ##       ##        ##   ##  ##       ##       ##
// ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##       ##       ##     ## ##       ##        ######
// ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##       ######### ##       ##             ##
// ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ###       ##    ## ##     ## ##       ##       ##    ##
// ##        #######  ##    ##  ######     ##    ####  #######  ##    ##        ######  ##     ## ######## ########  ######
//---------------------------------------------------------------------------

    tHq_analysis* theAnalysis = new tHq_analysis(thesamplelist, thesystlist, thechannellist, thevarlist, v_color, plot_extension, dir_ntuples, t_name, nofbin_templates, use_custom_colorPalette);
    theAnalysis->Set_Luminosity(set_luminosity); if(theAnalysis->stop_program) {return 1;}

    // theAnalysis->Train_BDT("", template_name, false);

    theAnalysis->Produce_Templates(template_name);

    // theAnalysis->Produce_Control_Histograms(template_name);

    // theAnalysis->Draw_Control_Plots(""); //If "tt" or "ttV", will use only these as bkg

    delete theAnalysis;
}
