#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Config.h"


using namespace std;

void Train_BDT()
{
    TMVA::Tools::Instance();

	(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 300;

	TFile* output_file = TFile::Open( "test.root", "RECREATE" );

    TString t_name = "Tree";

	TMVA::DataLoader *dataloader = new TMVA::DataLoader();
    dataloader->AddVariable("nJet25", 'F');
    dataloader->AddVariable("nJetEta1", 'F');

    TString inputfile = "~/tHq.root";
    TFile* file_input = TFile::Open(inputfile);
    TTree* tree = (TTree*) file_input->Get(t_name);

    dataloader->AddSignalTree(tree, 1);

    inputfile = "~/ttZ_LO.root";
    file_input = TFile::Open(inputfile);
    tree = (TTree*) file_input->Get(t_name);

    dataloader->AddBackgroundTree(tree, 1);

    dataloader->PrepareTrainingAndTestTree("", "", "nTrain_Signal=100:nTrain_Background=100:nTest_Signal=100:nTest_Background=100:SplitMode=Random:NormMode=NumEvents:!V");

    TMVA::Factory *factory = new TMVA::Factory("test", output_file, "V:!Silent:Color:DrawProgressBar:Correlations=True:AnalysisType=Classification");

    TString method_options= "!H:!V:NTrees=200:BoostType=Grad";

    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", method_options);

    // Train MVAs using the set of training events
    factory->TrainAllMethods();
    factory->TestAllMethods(); // ---- Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // ----- Evaluate and compare performance of all configured MVAs

    delete dataloader; dataloader = NULL;
	delete factory; factory = NULL;
    output_file->Close(); output_file = NULL;
    file_input->Close(); file_input = NULL;
}




int main()
{
    Train_BDT();
}
