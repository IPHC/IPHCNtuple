#include "../include/TTbarHiggsBTagEff.h"
#include "TSystem.h"

TTbarHiggsBTagEff::TTbarHiggsBTagEff() 
{

}

TTbarHiggsBTagEff::~TTbarHiggsBTagEff() 
{

    delete theHistoManager;

}


TTbarHiggsBTagEff::TTbarHiggsBTagEff(TString inputFileName, TChain *tree, TString sampleName, TString treeName, TString outputFileName, 
        //bool isdata, float xsec, float lumi, int nowe, 
        int nmax)
{    

    //
    //_isdata = isdata;
    //_xsec = xsec;
    //_lumi = lumi;
    //_nowe = nowe;
    _nmax = nmax;
    _outputFileName = outputFileName;
    //_sampleName = sampleName;


    tree = new TChain(treeName.Data());

    std::ifstream infile;
    infile.open(inputFileName.Data());
    std::string ifile = "";
    while( getline(infile, ifile) )
    {
        std::string fnameStr = std::string(ifile);

        tree->Add(fnameStr.c_str());

        std::cout << "file: " << fnameStr << std::endl;
    }   
    infile.close();

    Init(tree);

    theHistoManager = new HistoManager();


    TString outputFileNameRoot = _outputFileName+".root";
    _outputFile = new TFile(outputFileNameRoot.Data(), "recreate");  

}

void TTbarHiggsBTagEff::createHistograms()
{    

    _outputFile->cd();

    //
    theHistoManager->addHisto2D("pTvsEta_b",       "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaLoose_b",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaMedium_b", "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaTight_b",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);

    //  
    theHistoManager->addHisto2D("pTvsEta_c",       "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaLoose_c",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaMedium_c", "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaTight_c",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);

    //
    theHistoManager->addHisto2D("pTvsEta_l",       "", "", "",    20,    -2.5,    2.5,    20,    0,    200);   
    theHistoManager->addHisto2D("pTvsEtaLoose_l",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaMedium_l", "", "", "",    20,    -2.5,    2.5,    20,    0,    200);
    theHistoManager->addHisto2D("pTvsEtaTight_l",  "", "", "",    20,    -2.5,    2.5,    20,    0,    200);


}  


void TTbarHiggsBTagEff::writeHistograms()
{  
    _outputFile->cd();

    //
    TH2F* EffpTvsEtaLoose_b = (TH2F*)theHistoManager->getHisto2D("pTvsEtaLoose_b", "", "", "")->Clone("EffpTvsEtaLoose_b");
    EffpTvsEtaLoose_b->SetTitle("pTvsEtaLoose_b");
    EffpTvsEtaLoose_b->Divide(theHistoManager->getHisto2D("pTvsEta_b", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaLoose_b);

    TH2F* EffpTvsEtaMedium_b = (TH2F*)theHistoManager->getHisto2D("pTvsEtaMedium_b", "", "", "")->Clone("EffpTvsEtaMedium_b");
    EffpTvsEtaMedium_b->SetTitle("pTvsEtaMedium_b");
    EffpTvsEtaMedium_b->Divide(theHistoManager->getHisto2D("pTvsEta_b", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaMedium_b);

    //TH2F* EffpTvsEtaTight_b = (TH2F*)theHistoManager->getHisto2D("pTvsEtaTight_b", "", "", "")->Clone("EffpTvsEtaTight_b");
    //EffpTvsEtaTight_b->SetTitle("pTvsEtaTight_b");
    //EffpTvsEtaTight_b->Divide(theHistoManager->getHisto2D("pTvsEta_b", "", "", ""));
    //theHistoManager->addHisto2D(EffpTvsEtaTight_b);

    //
    TH2F* EffpTvsEtaLoose_c = (TH2F*)theHistoManager->getHisto2D("pTvsEtaLoose_c", "", "", "")->Clone("EffpTvsEtaLoose_c");
    EffpTvsEtaLoose_c->SetTitle("pTvsEtaLoose_c");
    EffpTvsEtaLoose_c->Divide(theHistoManager->getHisto2D("pTvsEta_c", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaLoose_c);

    TH2F* EffpTvsEtaMedium_c = (TH2F*)theHistoManager->getHisto2D("pTvsEtaMedium_c", "", "", "")->Clone("EffpTvsEtaMedium_c");
    EffpTvsEtaMedium_c->SetTitle("pTvsEtaMedium_c");
    EffpTvsEtaMedium_c->Divide(theHistoManager->getHisto2D("pTvsEta_c", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaMedium_c);

    //TH2F* EffpTvsEtaTight_c = (TH2F*)theHistoManager->getHisto2D("pTvsEtaTight_c", "", "", "")->Clone("EffpTvsEtaTight_c");
    //EffpTvsEtaTight_c->SetTitle("pTvsEtaTight_c");
    //EffpTvsEtaTight_c->Divide(theHistoManager->getHisto2D("pTvsEta_c", "", "", ""));
    //theHistoManager->addHisto2D(EffpTvsEtaTight_c);

    //
    TH2F* EffpTvsEtaLoose_l = (TH2F*)theHistoManager->getHisto2D("pTvsEtaLoose_l", "", "", "")->Clone("EffpTvsEtaLoose_l");
    EffpTvsEtaLoose_l->SetTitle("pTvsEtaLoose_l");
    EffpTvsEtaLoose_l->Divide(theHistoManager->getHisto2D("pTvsEta_l", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaLoose_l);

    TH2F* EffpTvsEtaMedium_l = (TH2F*)theHistoManager->getHisto2D("pTvsEtaMedium_l", "", "", "")->Clone("EffpTvsEtaMedium_l");
    EffpTvsEtaMedium_l->SetTitle("pTvsEtaMedium_l");
    EffpTvsEtaMedium_l->Divide(theHistoManager->getHisto2D("pTvsEta_l", "", "", ""));
    theHistoManager->addHisto2D(EffpTvsEtaMedium_l);

    //TH2F* EffpTvsEtaTight_l = (TH2F*)theHistoManager->getHisto2D("pTvsEtaTight_l", "", "", "")->Clone("EffpTvsEtaTight_l");
    //EffpTvsEtaTight_l->SetTitle("pTvsEtaTight_l");
    //EffpTvsEtaTight_l->Divide(theHistoManager->getHisto2D("pTvsEta_l", "", "", ""));
    //theHistoManager->addHisto2D(EffpTvsEtaTight_l);

    //
    std::vector<TH2F*> the2DHisto =  theHistoManager->getHisto2D_list();
    for(unsigned int i=0; i<the2DHisto.size(); i++)  the2DHisto[i]->Write();



    _outputFile->Close();


}


void TTbarHiggsBTagEff::Init(TChain *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;

    vJet = new std::vector<Jet>();

    fChain->SetBranchAddress("Jet",  &vJet);

}

void TTbarHiggsBTagEff::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();

    int nentries_max = nentries;
    if ( _nmax != -1 && _nmax < nentries ) nentries_max = _nmax;

    std::cout << "Number of input events = " << nentries_max << std::endl;

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries_max; jentry++) 
    {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;

        if(jentry%10000 == 0) std::cout << "number of processed events " << jentry << std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        ////////////////////////////////////////////////
        //
        ////////////////////////////////////////////////



        for(unsigned int ijet=0; ijet < vJet->size() ; ijet++)
        {             
            if ( vJet->at(ijet).pt() > 20. && fabs(vJet->at(ijet).eta()) < 2.5 )
            { 
                std::cout << "flavour "<< vJet->at(ijet).jet_hadronFlavour() << std::endl;

                if (vJet->at(ijet).jet_hadronFlavour()==5) 
                {
                    theHistoManager->fillHisto2D("pTvsEta_b", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    if (vJet->at(ijet).CSVv2() >= 0.5426) theHistoManager->fillHisto2D("pTvsEtaLoose_b", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    if (vJet->at(ijet).CSVv2() >= 0.8484) theHistoManager->fillHisto2D("pTvsEtaMedium_b", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    //if (vJet->at(ijet).CSVv2()) theHistoManager->fillHisto2D("pTvsEtaTight_b", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                }

                if (vJet->at(ijet).jet_hadronFlavour()==4) 
                { 
                    theHistoManager->fillHisto2D("pTvsEta_c", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1);
                    if (vJet->at(ijet).CSVv2() >= 0.5426) theHistoManager->fillHisto2D("pTvsEtaLoose_c", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    if (vJet->at(ijet).CSVv2() >= 0.8484) theHistoManager->fillHisto2D("pTvsEtaMedium_c", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    //if (vJet->at(ijet).CSVv2()) theHistoManager->fillHisto2D("pTvsEtaTight_c", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                }

                if (vJet->at(ijet).jet_hadronFlavour()==0)
                {  

                    std::cout << "light"<<std::endl;

                    theHistoManager->fillHisto2D("pTvsEta_l", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1);
                    if (vJet->at(ijet).CSVv2() >= 0.5426) theHistoManager->fillHisto2D("pTvsEtaLoose_l", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    if (vJet->at(ijet).CSVv2() >= 0.8484) theHistoManager->fillHisto2D("pTvsEtaMedium_l", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                    //if (vJet->at(ijet).jet_hadronFlavour()) theHistoManager->fillHisto2D("pTvsEtaTight_l", "", "", "", vJet->at(ijet).eta(), vJet->at(ijet).pt(),1); 
                }


            }//pT/eta	     

        }


    }


}
