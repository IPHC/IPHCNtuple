//by Nicolas Tonon (IPHC)

#include "tHq_analysis.h"

using namespace std;

//---------------------------------------------------------------------------
// ####    ##    ##    ####    ########
//  ##     ###   ##     ##        ##
//  ##     ####  ##     ##        ##
//  ##     ## ## ##     ##        ##
//  ##     ##  ####     ##        ##
//  ##     ##   ###     ##        ##
// ####    ##    ##    ####       ##
//---------------------------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

//Overloaded constructor
tHq_analysis::tHq_analysis(vector<TString> thesamplelist, vector<TString> thesystlist, vector<TString> thechannellist, vector<TString> thevarlist, vector<int> thecolorlist, TString theplotextension, TString dirntuplename, TString treename, int nofbin_templates, bool use_custom_colors)
{
	TH1::SetDefaultSumw2();
	gStyle->SetErrorX(0.);

	mkdir("outputs", 0777);
	mkdir("plots", 0777);

	stop_program = false;

	dir_ntuples = dirntuplename;
	t_name = treename;

	luminosity_rescale = 1;

	use_custom_colorPalette = use_custom_colors;

	plot_extension = theplotextension;

	nbin = nofbin_templates;

	sample_list.resize(thesamplelist.size());
	for(int i=0; i<thesamplelist.size(); i++)
	{
		sample_list[i] = thesamplelist[i];
	}

	syst_list.resize(thesystlist.size());
	for(int i=0; i<thesystlist.size(); i++)
	{
		syst_list[i] = thesystlist[i];
	}

	channel_list.resize(thechannellist.size());
	for(int i=0; i<thechannellist.size(); i++)
	{
		channel_list[i] = thechannellist[i];
	}

	var_list.resize(thevarlist.size());
	var_list_floats.resize(thevarlist.size());
	for(int i=0; i<thevarlist.size(); i++)
	{
		var_list[i] = thevarlist[i];
		var_list_floats[i] = 0;
	}

	color_list.resize(thecolorlist.size());
	for(int i=0; i<thecolorlist.size(); i++)
	{
		color_list[i] = thecolorlist[i];
	}

	if(use_custom_colorPalette) {Set_Custom_ColorPalette(v_custom_colors, color_list);}


	if(plot_extension != ".png" && plot_extension != ".pdf") {cout<<FRED("Wrong plot_extension value ! Abort !")<<endl; stop_program=true;}
	if(nofbin_templates < 0) {cout<<FRED("Wrong nbins value ! Abort !")<<endl; stop_program=true;}
}

tHq_analysis::~tHq_analysis()
{
	if(use_custom_colorPalette)
	{
		for(int icolor=0; icolor<color_list.size(); icolor++)
		{
			delete v_custom_colors[icolor];
		}
	}
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/**
 * Compute the luminosity re-scaling factor (MC),  to be used thoughout the code
 * @param desired_luminosity [Value of the desired lumi in fb-1]
 */
void tHq_analysis::Set_Luminosity(double desired_luminosity)
{
	double current_luminosity = 35.862; //Moriond 2017 //CHANGED
	this->luminosity_rescale = desired_luminosity / current_luminosity;

	if(luminosity_rescale !=1 )
	{
		cout<<endl<<BOLD(FBLU("##################################"))<<endl;
		cout<<"--- Using luminosity scale factor : "<<desired_luminosity<<" / "<<current_luminosity<<" = "<<luminosity_rescale<<" ! ---"<<endl;
		cout<<BOLD(FBLU("##################################"))<<endl<<endl;
	}
}







//---------------------------------------------------------------------------
// ########    ########        ###       ####    ##    ##    ####    ##    ##     ######
//    ##       ##     ##      ## ##       ##     ###   ##     ##     ###   ##    ##    ##
//    ##       ##     ##     ##   ##      ##     ####  ##     ##     ####  ##    ##
//    ##       ########     ##     ##     ##     ## ## ##     ##     ## ## ##    ##   ####
//    ##       ##   ##      #########     ##     ##  ####     ##     ##  ####    ##    ##
//    ##       ##    ##     ##     ##     ##     ##   ###     ##     ##   ###    ##    ##
//    ##       ##     ##    ##     ##    ####    ##    ##    ####    ##    ##     ######
//---------------------------------------------------------------------------


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

/*
* Train, test and evaluate the BDT with signal and bkg MC
* (Use a DataLoader ; seems that it is required in newer root versions)
* @param channel
* @param type           [tt or ttV]
* @param write_ranking_info [Save variable ranking info in file or not (must set to false if optimizing BDT variable list)]
 */
void tHq_analysis::Train_BDT(TString channel, TString type, bool write_ranking_info)
{
	cout<<endl<<BOLD(FYEL("##################################"))<<endl;
	cout<<FYEL("---TRAINING ---")<<endl;
	cout<<BOLD(FYEL("##################################"))<<endl<<endl;

	// cout<<endl<<BOLD(FGRN("-- USE ONLY ttZ / WZ / ZZ / tZq FOR TRAINING!"))<<endl;


	//---------------------------------------------------------------
	//--- Could modify here the name of dir. containing the BDT weights (default = "weights")
	TMVA::gConfig().GetIONames().fWeightFileDir = type;

    // This loads the TMVA libraries
    TMVA::Tools::Instance();

	//Allows to bypass a protection in TMVA::TransformationHandler, cf. description in source file:
	// if there are too many input variables, the creation of correlations plots blows up memory and basically kills the TMVA execution --> avoid above critical number (which can be user defined)
	(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 300;

	TString output_file_name = "outputs/BDT";
	if(channel != "") {output_file_name+= "_" + channel;}
	// output_file_name+= this->filename_suffix;
	output_file_name+= "_training.root";

	TFile* output_file = TFile::Open( output_file_name, "RECREATE" );

	// Create the factory object
	// TMVA::Factory* factory = new TMVA::Factory(type.Data(), output_file, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );
	TMVA::DataLoader *dataloader = new TMVA::DataLoader("weights_tHq"); //If no TString given in arg, will store weights in : default/weights/...

	// Define the input variables that shall be used for the MVA training
	for(int i=0; i<var_list.size(); i++)
	{
		dataloader->AddVariable(var_list[i].Data(), 'F');
	}

	//Choose if the cut variables are used in BDT or not
	//NOTE : spectator vars are not used for training/evalution, but possible to check their correlations, etc.
	// for(int i=0; i<v_cut_name.size(); i++)
	{
		// cout<<"Is "<<v_cut_name[i]<<" used ? "<<(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("=="))<<endl;

		//if we ask "var == x", all the selected events will have the value x, so can't use it as discriminant variable !
		// if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {dataloader->AddVariable(v_cut_name[i].Data(), 'F');}
		// else {dataloader->AddSpectator(v_cut_name[i].Data(), v_cut_name[i].Data(), 'F');}
	}
	for(int i=0; i<v_add_var_names.size(); i++)
	{
		dataloader->AddSpectator(v_add_var_names[i].Data(), v_add_var_names[i].Data(), 'F');
	}

    TFile *f(0);
 	std::vector<TFile *> files_to_close;

	for(int isample=0; isample<sample_list.size(); isample++)
    {
		if(type == "tt")
		{
			if(!sample_list[isample].Contains("tHq") && !sample_list[isample].Contains("tHW") && !sample_list[isample].Contains("ttbar") ) {continue;}
		}
		else if(type == "ttV")
		{
			if(!sample_list[isample].Contains("tHq") && !sample_list[isample].Contains("tHW") && !sample_list[isample].Contains("ttZ") && !sample_list[isample].Contains("ttW") ) {continue;}
		}


        // Read training and test data
        // --- Register the training and test trees
		TString inputfile;
		inputfile = dir_ntuples + "/" + sample_list[isample] + ".root";

	    TFile* file_input = 0;
		file_input = TFile::Open(inputfile.Data() );
		if(!file_input) {cout<<BOLD(FRED(<<inputfile.Data()<<" not found!"))<<endl; continue;}
		files_to_close.push_back(file_input);

		TTree* tree = 0;
		tree = (TTree*) file_input->Get(t_name.Data());
		if(tree==0)
		{
			cout<<BOLD(FRED("ERROR :"))<<" file "<<inputfile.Data()<<" --> *tree = 0 !"<<endl; continue;
		}


        // global event weights per tree (see below for setting event-wise weights)
        Double_t signalWeight     = 1.0;
        Double_t backgroundWeight = 1.0;

    //-- Choose between absolute/relative weights for training

		if(sample_list[isample].Contains("tHq") || sample_list[isample].Contains("tHW") )
		{
			dataloader->AddSignalTree(tree, signalWeight);
			dataloader->SetSignalWeightExpression("weight");
			// dataloader->SetSignalWeightExpression("fabs(weight)");
		}
		else
		{
			dataloader->AddBackgroundTree(tree, backgroundWeight);
			dataloader->SetBackgroundWeightExpression("weight");
			// dataloader->SetBackgroundWeightExpression("fabs(weight)");
		}
		cout<<FYEL("-- Using *RELATIVE* weights --")<<endl;
		// cout<<FYEL("-- Using *ABSOLUTE* weights --")<<endl;
    }


	// Apply additional cuts on the signal and background samples (can be different)
	TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

	if(channel != "all" && channel != "")
    {
        // if(channel == "uuu")         		{mycuts = "Channel==0"; mycutb = "Channel==0";}
        // else if(channel == "uue" )     		{mycuts = "Channel==1"; mycutb = "Channel==1";}
        // else if(channel == "eeu"  )     	{mycuts = "Channel==2"; mycutb = "Channel==2";}
		// else if(channel == "eee"   )     	{mycuts = "Channel==3"; mycutb = "Channel==3";}
		// else 								{cout << "WARNING : wrong channel name while training " << endl;}
    }

//--------------------------------
//--- Apply cuts during training
	TString tmp = "";
/*
	for(int ivar=0; ivar<v_cut_name.size(); ivar++)
	{
		if(v_cut_def[ivar] != "")
		{
			if(!v_cut_def[ivar].Contains("&&") && !v_cut_def[ivar].Contains("||")) {tmp+= v_cut_name[ivar] + v_cut_def[ivar];} //If cut contains only 1 condition
			else if(v_cut_def[ivar].Contains("&&") && v_cut_def[ivar].Contains("||")) {cout<<BOLD(FRED("ERROR ! Wrong cut definition !"))<<endl;}
			else if(v_cut_def[ivar].Contains("&&") )//If '&&' in the cut, break it in 2
			{
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
				tmp+= " && ";
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
			}
			else if(v_cut_def[ivar].Contains("||") )//If '||' in the cut, break it in 2
			{
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
				tmp+= " || ";
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
			}
		}

		//Complicated way of concatenating the TStrings
		if((ivar+1) < v_cut_name.size() && v_cut_def[ivar+1] != "")
		{
			for(int i=0; i<ivar+1; i++)
			{
				if(v_cut_def[i] != "") {tmp += " && "; break;}
			}
		}

		//cout<<"tmp = "<<tmp<<endl;
	}

	// cout<<"Total cut chain : "<<tmp<<endl;

	if(tmp != "") {mycuts+= tmp; mycutb+= tmp;}
*/

//--------------------------------
	if(mycuts != mycutb) {cout<<__LINE__<<FRED("PROBLEM : cuts are different for signal and background ! If this is normal, modify code -- Abort")<<endl; return;}

    // Tell the factory how to use the training and testing events    //
    // If no numbers of events are given, half of the events in the tree are used for training, and the other half for testing:
	dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

	// TString method_title = channel + this->filename_suffix; //So that the output weights are labelled differently for each channel
	TString method_title = "all"; //So that the output weights are labelled differently for each channel

	//--- Boosted Decision Trees -- Choose method
	TMVA::Factory *factory = new TMVA::Factory("BDT"+type, output_file, "V:!Silent:Color:DrawProgressBar::AnalysisType=Classification");


	//--- BDT METHOD
    factory->BookMethod(dataloader, TMVA::Types::kBDT, method_title,"!H:!V:NTrees=100:nCuts=100:MaxDepth=2:BoostType=Grad:Shrinkage=0.4");

	// TMVA Multilayer Perceptron -- Example from TMVA
	// factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

//--------------------------------------
	output_file->cd();


	mkdir("outputs/Rankings", 0777); //Dir. containing variable ranking infos


	TString ranking_file_path = "outputs/rank.txt";

	if(write_ranking_info) cout<<endl<<endl<<endl<<FBLU("NB : Temporarily redirecting standard output to file '"<<ranking_file_path<<"' in order to save Ranking Info !!")<<endl<<endl<<endl;

	std::ofstream out("ranking_info_tmp.txt"); //Temporary name
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	if(write_ranking_info) std::cout.rdbuf(out.rdbuf()); //redirect std::cout to text file --> Ranking info will be saved !

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

	if(write_ranking_info) std::cout.rdbuf(coutbuf); //reset to standard output again

	//-- NB : Test & Evaluation recap in the output files
    factory->TestAllMethods(); // ---- Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // ----- Evaluate and compare performance of all configured MVAs


    // --------------------------------------------------------------
    // Save the output
    output_file->Close();
    std::cout << "==> Wrote root file: " << output_file->GetName() << std::endl;
    std::cout << "==> TMVA is done!" << std::endl;

	for(unsigned int i=0; i<files_to_close.size(); i++) {files_to_close[i]->Close(); delete files_to_close[i];}


	if(write_ranking_info)
	{
		MoveFile("./ranking_info_tmp.txt", ranking_file_path);
		Extract_Ranking_Info(ranking_file_path, type, channel); //Extract only ranking info from TMVA output
	}
	else {system("rm ./ranking_info_tmp.txt");} //Else remove the temporary ranking file

	delete dataloader; dataloader = NULL;
	delete factory; factory = NULL;
	output_file->Close(); output_file = NULL;


	return;
}













//---------------------------------------------------------------------------
//  ######  ########  ########    ###    ######## ########       ######## ######## ##     ## ########  ##          ###    ######## ########  ######
// ##    ## ##     ## ##         ## ##      ##    ##                ##    ##       ###   ### ##     ## ##         ## ##      ##    ##       ##    ##
// ##       ##     ## ##        ##   ##     ##    ##                ##    ##       #### #### ##     ## ##        ##   ##     ##    ##       ##
// ##       ########  ######   ##     ##    ##    ######            ##    ######   ## ### ## ########  ##       ##     ##    ##    ######    ######
// ##       ##   ##   ##       #########    ##    ##                ##    ##       ##     ## ##        ##       #########    ##    ##             ##
// ##    ## ##    ##  ##       ##     ##    ##    ##                ##    ##       ##     ## ##        ##       ##     ##    ##    ##       ##    ##
//  ######  ##     ## ######## ##     ##    ##    ########          ##    ######## ##     ## ##        ######## ##     ##    ##    ########  ######
//---------------------------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/**
 * Uses output from training (weights, ...) and input ntuples to create distributions of the BDT discriminants
 * @param  template_name         BDT, BDTttZ or mTW
 * @param  fakes_from_data       If true, use fakes from data sample
 * @param  real_data             If true, use real data ; else, looks for pseudodata
 * @param  cut_on_BDT            Cut value on BDT (-> cran create templates in CR)
 */
void tHq_analysis::Produce_Templates(TString template_name)
{
	cout<<endl<<BOLD(FYEL("##################################"))<<endl;
	if(template_name == "tt" || template_name == "ttV") {cout<<FYEL("--- Producing "<<template_name<<" Templates ---")<<endl;}
	else {cout<<BOLD(FRED("--- ERROR : invalid template_name value ! Exit !"))<<endl; cout<<"Valid names are : BDT/BDTttZ/mTW !"<<endl; return;}
	cout<<BOLD(FYEL("##################################"))<<endl<<endl;

	TString output_file_name = "outputs/Templates_" + template_name;
	output_file_name+= ".root";

	TFile* file_output = TFile::Open( output_file_name, "RECREATE" );

	reader = new TMVA::Reader( "!Color:!Silent" );

	// Name & adress of local variables which carry the updated input values during the event loop
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used -- same order
	for(int i=0; i<var_list.size(); i++)
	{
		reader->AddVariable(var_list[i].Data(), &var_list_floats[i]); //cout<<"Added variable "<<var_list[i]<<endl;
	}

	/*
	for(int i=0; i<v_cut_name.size(); i++)
	{
		if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {reader->AddVariable(v_cut_name[i].Data(), &v_cut_float[i]);}
		else
		{
			// cout<<v_cut_name[i]<<endl;
			reader->AddSpectator(v_cut_name[i].Data(), &v_cut_float[i]);
		}
	}
	for(int i=0; i<v_add_var_names.size(); i++)
	{
		// cout<<v_add_var_names[i]<<endl;
		reader->AddSpectator(v_add_var_names[i].Data(), &v_add_var_floats[i]);
	}*/


	// --- Book the MVA methods
	TString dir    = "weights_tHq/" + template_name + "/";
	// TString dir    = "weights_2016/";

	TString MVA_method_name = "";
	TString weightfile = "";

	for(int ichan=0; ichan<channel_list.size(); ichan++) //Book the method for each channel (separate BDTs)
	{
		TString template_name_MVA = template_name + "_all"; //---

		MVA_method_name = template_name_MVA + " method";
		weightfile = dir + "BDT" + template_name_MVA + ".weights.xml";
		// weightfile = dir + "thq_vs_ttv_3l_BDTG.weights.xml";

		if(!Check_File_Existence(weightfile) ) {cout<<BOLD(FRED("Weight file "<<weightfile<<" not found ! Abort"))<<endl; return;}

		reader->BookMVA(MVA_method_name, weightfile);
	}


	TFile* file_input;
	TTree* tree(0);

	// TH1F *hist_uuu = 0, *hist_uue = 0, *hist_eeu = 0, *hist_eee = 0;
	TH1F *hist_all = 0;

	// --- Systematics loop
	for(int isyst=0; isyst<syst_list.size(); isyst++)
	{
		cout<<endl<<endl<<FGRN("Systematic "<<syst_list[isyst]<<" ("<<isyst+1<<"/"<<syst_list.size()<<") :")<<endl;

		//Loop on samples, syst., events ---> Fill histogram/channel --> Write()
		for(int isample=0; isample<sample_list.size(); isample++)
		{

			TString inputfile = dir_ntuples + "/" + sample_list[isample] + ".root";


			if(!Check_File_Existence(inputfile) ) {cout<<inputfile.Data()<<" not found!"<<endl; continue;}
			file_input = TFile::Open( inputfile.Data() );


			// Book output histograms
			if (template_name.Contains("tt")) //create histogram for each channel (-1 = bkg, +1 = signal)
			{
				hist_all     = new TH1F( template_name+"_all",           template_name+"_all",           nbin, -1, 1 );

				// hist_uuu     = new TH1F( (template_name+"_uuu").Data(),           (template_name+"_uuu").Data(),           nbin, -1, 1 );
				// hist_uue     = new TH1F( (template_name+"_uue").Data(),           (template_name+"_uue").Data(),           nbin, -1, 1 );
				// hist_eeu     = new TH1F( (template_name+"_eeu").Data(),           (template_name+"_eeu").Data(),           nbin, -1, 1 );
				// hist_eee     = new TH1F( (template_name+"_eee").Data(),           (template_name+"_eee").Data(),           nbin, -1, 1 );
			}
			else {cout<<FRED("ERROR : wrong template name !")<<endl; return;}


			tree = 0;
			tree = (TTree*) file_input->Get(t_name.Data());

			if(!tree && syst_list[isyst]=="") {cout<<BOLD(FRED("ERROR : nominal tree not found ! Abort"))<<endl; return;}
			else if(!tree) {cout<<BOLD(FRED("ERROR : tree "<<syst_list[isyst]<<" not found in "<<inputfile<< "! Skip this systematic"))<<endl; break;}

//--- Prepare the event tree -- Set Branch Addresses
//WARNING : the last SetBranchAddress overrides the previous ones !! Be careful not to associate branches twice !
			/*
			for(int i=0; i<v_add_var_names.size(); i++)
			{
				tree->SetBranchAddress(v_add_var_names[i].Data(), &v_add_var_floats[i]);
			}
			*/
			for(int i=0; i<var_list.size(); i++)
			{
				tree->SetBranchAddress(var_list[i].Data(), &var_list_floats[i]);
			}
			/*
			for(int i=0; i<v_cut_name.size(); i++)
			{
				tree->SetBranchAddress(v_cut_name[i].Data(), &v_cut_float[i]);
			}
			*/

			//Dedicated variables, easier to access in event loop
			float theVar = -666;
			// tree->SetBranchAddress("", &theVar);
			// float i_channel = 9; tree->SetBranchAddress("Channel", &i_channel);


			float weight;
			//For all other systematics, only the events weights change
			tree->SetBranchAddress("weight", &weight);


			cout<<endl<< "--- "<<sample_list[isample]<<" : Processing: " << tree->GetEntries() << " events" << std::endl;
//------------------------------------------------------------
// --- START EVENT LOOP ---

			// int n_entries = 100;
			int n_entries = tree->GetEntries();

			for(int ievt=0; ievt<n_entries; ievt++)
			{
				if(ievt%10000==0) cout<<ievt<<" / "<<n_entries<<endl;

				weight = 0;
				theVar=-666;
				// i_channel = 9;

				tree->GetEntry(ievt);

				// NOTE -- Can simulate signal strength here
				// if(sample_list[isample].Contains("tZq")) {weight*=3;}

				bool isChannelToKeep = true;
				if(!isChannelToKeep) continue;

//------------------------------------------------------------
//------------------------------------------------------------
//---- APPLY CUTS HERE  ----
				float cut_tmp = 0; bool pass_all_cuts = true;
				/*
				for(int ivar=0; ivar<v_cut_name.size(); ivar++)
				{
					if(v_cut_def[ivar] == "") {continue;}

					else if((sample_list[isample].Contains("Fakes") || sample_list[isample].Contains("DY") || sample_list[isample].Contains("TT") || sample_list[isample].Contains("WW")) && (v_cut_name[ivar] == "AdditionalMuonIso" || v_cut_name[ivar] == "AdditionalEleIso") && v_cut_def[ivar].Contains("<") )
					{
						cout<<endl<<endl<<BOLD(FYEL("Not applying '<' cuts on isolation for Fakes samples !"))<<endl<<endl<<endl;

						continue;
					}

					//Can't set Branch address of same branch to 2 different variables (only last one will be filled)
					//Since I define myself a 'mTW' variable, it is a problem if one of the cuts is on mTW (because the associated float won't be filled)
					//---> If v_cut_name[i] == "mTW", need to make sure that we use the variable which is filled !
					if(v_cut_name[ivar] == "mTW") {v_cut_float[ivar] = mTW;}
					//Idem for channel value
					if(v_cut_name[ivar] == "Channel") {v_cut_float[ivar] = i_channel;}

					// cout<<v_cut_name[ivar]<<" "<<v_cut_float[ivar]<<endl;

					if(!v_cut_def[ivar].Contains("&&") && !v_cut_def[ivar].Contains("||")) //If cut contains only 1 condition
					{
						cut_tmp = Find_Number_In_TString(v_cut_def[ivar]);
						if(v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] < cut_tmp)		 {pass_all_cuts = false; break;}
						else if(v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] > cut_tmp)	 {pass_all_cuts = false; break;}
						else if(v_cut_def[ivar].Contains(">") && !v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] <= cut_tmp)	 {pass_all_cuts = false; break;}
						else if(v_cut_def[ivar].Contains("<") && !v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] >= cut_tmp) 	 {pass_all_cuts = false; break;}
						else if(v_cut_def[ivar].Contains("==") && v_cut_float[ivar] != cut_tmp)  {pass_all_cuts = false; break;}
						else if(v_cut_def[ivar].Contains("!=") && v_cut_float[ivar] == cut_tmp)  {pass_all_cuts = false; break;}
					}
					else if(v_cut_def[ivar].Contains("&&") && v_cut_def[ivar].Contains("||")) {cout<<BOLD(FRED("ERROR ! Wrong cut definition !"))<<endl;}
					else if(v_cut_def[ivar].Contains("&&") )//If '&&' in the cut def, break it in 2
					{
						TString cut1 = Break_Cuts_In_Two(v_cut_def[ivar]).first;
						TString cut2 = Break_Cuts_In_Two(v_cut_def[ivar]).second;
						//CUT 1
						cut_tmp = Find_Number_In_TString(cut1);
						if(cut1.Contains(">=") && v_cut_float[ivar] < cut_tmp)			 {pass_all_cuts = false; break;}
						else if(cut1.Contains("<=") && v_cut_float[ivar] > cut_tmp)		 {pass_all_cuts = false; break;}
						else if(cut1.Contains(">") && !v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] <= cut_tmp)		 {pass_all_cuts = false; break;}
						else if(cut1.Contains("<") && !v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] >= cut_tmp) 	 {pass_all_cuts = false; break;}
						else if(cut1.Contains("==") && v_cut_float[ivar] != cut_tmp) 	 {pass_all_cuts = false; break;}
						else if(cut1.Contains("!=") && v_cut_float[ivar] == cut_tmp) 	 {pass_all_cuts = false; break;}
						//CUT 2
						cut_tmp = Find_Number_In_TString(cut2);
						if(cut2.Contains(">=") && v_cut_float[ivar] < cut_tmp)			 {pass_all_cuts = false; break;}
						else if(cut2.Contains("<=") && v_cut_float[ivar] > cut_tmp)		 {pass_all_cuts = false; break;}
						else if(cut2.Contains(">") && !v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] <= cut_tmp)		 {pass_all_cuts = false; break;}
						else if(cut2.Contains("<") && !v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] >= cut_tmp) 	 {pass_all_cuts = false; break;}
						else if(cut2.Contains("==") && v_cut_float[ivar] != cut_tmp) 	 {pass_all_cuts = false; break;}
						else if(cut2.Contains("!=") && v_cut_float[ivar] == cut_tmp) 	 {pass_all_cuts = false; break;}
					}
					else if(v_cut_def[ivar].Contains("||") )//If '||' in the cut def, break it in 2
					{
						TString cut1 = Break_Cuts_In_Two(v_cut_def[ivar]).first;
						TString cut2 = Break_Cuts_In_Two(v_cut_def[ivar]).second;

						bool pass_cut1 = true; bool pass_cut2 = true; //Need to pass at least 1 cut

						//CUT 1
						cut_tmp = Find_Number_In_TString(cut1);
						if(cut1.Contains(">=") && v_cut_float[ivar] < cut_tmp)			 {pass_cut1 = false;}
						else if(cut1.Contains("<=") && v_cut_float[ivar] > cut_tmp)		 {pass_cut1 = false;}
						else if(cut1.Contains(">") && !v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] <= cut_tmp)		 {pass_cut1 = false;}
						else if(cut1.Contains("<") && !v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] >= cut_tmp) 	 {pass_cut1 = false;}
						else if(cut1.Contains("==") && v_cut_float[ivar] != cut_tmp) 	 {pass_cut1 = false;}
						else if(cut1.Contains("!=") && v_cut_float[ivar] == cut_tmp) 	 {pass_cut1 = false;}
						//CUT 2
						cut_tmp = Find_Number_In_TString(cut2);
						if(cut2.Contains(">=") && v_cut_float[ivar] < cut_tmp)			 {pass_cut2 = false;}
						else if(cut2.Contains("<=") && v_cut_float[ivar] > cut_tmp)		 {pass_cut2 = false;}
						else if(cut2.Contains(">") && !v_cut_def[ivar].Contains(">=") && v_cut_float[ivar] <= cut_tmp)		 {pass_cut2 = false;}
						else if(cut2.Contains("<") && !v_cut_def[ivar].Contains("<=") && v_cut_float[ivar] >= cut_tmp) 	 {pass_cut2 = false;}
						else if(cut2.Contains("==") && v_cut_float[ivar] != cut_tmp) 	 {pass_cut2 = false;}
						else if(cut2.Contains("!=") && v_cut_float[ivar] == cut_tmp) 	 {pass_cut2 = false;}

						if(!pass_cut1 && !pass_cut2) {pass_all_cuts = false; break;}
					}
				}*/


				if(!pass_all_cuts) {continue;}

				TString chan_name = "all";
				// if(i_channel == 0) chan_name = "uuu";
				// else if(i_channel == 1) chan_name = "uue";
				// else if(i_channel == 2) chan_name = "eeu";
				// else if(i_channel == 3) chan_name = "eee";


				double xmax_h = hist_all->GetXaxis()->GetXmax();
				int lastbin_h = hist_all->GetNbinsX();
				double mva_value = -999;

				mva_value = reader->EvaluateMVA(MVA_method_name);
				if(mva_value < xmax_h ) {hist_all->Fill(mva_value, weight);}
				else {Fill_Last_Bin_TH1F(hist_all, weight);} //Put overflow in last bin (no info lost)
			} //end entries loop


//---------

			//Re-scale to desired luminosity, unless it's data
			if(sample_list[isample] != "Data")
			{
				hist_all->Scale(luminosity_rescale);

				// hist_uuu->Scale(luminosity_rescale); hist_uue->Scale(luminosity_rescale); hist_eeu->Scale(luminosity_rescale); hist_eee->Scale(luminosity_rescale);
			}


			// --- Write histograms
			file_output->cd();

			//NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)]
			TString output_histo_name = "";
			TString syst_name = "";

			TString sample_name = sample_list[isample];
			// if(real_data && sample_list[isample] == "Data") {sample_name = "DATA";} //THETA CONVENTION
			if(sample_list[isample] == "Data")
			{
				sample_name = "data_obs";
			}

			TString template_name_output = template_name;
			// template_name_output = "BDT2l";

			output_histo_name = template_name_output+"_all__" + sample_name + syst_name;
			hist_all->Write(output_histo_name);

			// output_histo_name = template_name_output+"_uuu__" + sample_name + syst_name;
			// hist_uuu->Write(output_histo_name.Data());
			// output_histo_name = template_name_output+"_uue__" + sample_name + syst_name;
			// hist_uue->Write(output_histo_name.Data());
			// output_histo_name = template_name_output+"_eeu__" + sample_name + syst_name;
			// hist_eeu->Write(output_histo_name.Data());
			// output_histo_name = template_name_output+"_eee__" + sample_name + syst_name;
			// hist_eee->Write(output_histo_name.Data());

			cout<<"Done with "<<sample_list[isample]<<" sample"<<endl;

			delete tree; tree = NULL;
			delete hist_all; hist_all = NULL;
			// delete hist_uuu; delete hist_uue; delete hist_eeu; delete hist_eee; //Free memory
			file_input->Close(); file_input = NULL;
		} //end sample loop

		cout<<"Done with syst : "<<syst_list[isyst]<<endl;
	} 	//end syst loop


	cout<<endl<<FYEL("==> Created root file: "<<file_output->GetName()<<" containing the output histograms")<<endl<<endl;
	delete file_output; file_output = NULL;
	delete reader; reader = NULL;


	return;
}












//---------------------------------------------------------------------------
//  ######      ########            ##     ## ####  ######  ########  ######
// ##    ##     ##     ##           ##     ##  ##  ##    ##    ##    ##    ##
// ##           ##     ##           ##     ##  ##  ##          ##    ##
// ##           ########            #########  ##   ######     ##     ######
// ##           ##   ##             ##     ##  ##        ##    ##          ##
// ##    ## ### ##    ##  ###       ##     ##  ##  ##    ##    ##    ##    ##
//  ######  ### ##     ## ###       ##     ## ####  ######     ##     ######
//---------------------------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void tHq_analysis::Produce_Control_Histograms()
{
	cout<<endl<<BOLD(FYEL("##################################"))<<endl;
	cout<<FYEL("--- Create Control Histograms ---")<<endl;
	cout<<BOLD(FYEL("##################################"))<<endl<<endl;

	TString output_file_name = "outputs/Control_Histograms_tHq.root";
	TFile* f_output = TFile::Open( output_file_name, "RECREATE" );

	int binning_choice = 0;

	for(int isample=0; isample<sample_list.size(); isample++)
	{
		cout<<"--- Sample "<<sample_list[isample]<<endl;

		TString input_file_name = dir_ntuples + "/" + sample_list[isample] + ".root";
		if(!Check_File_Existence(input_file_name) ) {cout<<FRED("Error ! File "<<input_file_name<<" not found !")<<endl;}
		TFile* f_input = TFile::Open(input_file_name);
		TTree* tree = (TTree*) f_input->Get("Tree");

		for(int ivar=0; ivar<var_list.size(); ivar++)
		{
			for(int ichan=0; ichan<channel_list.size(); ichan++)
			{
				for(int isyst=0; isyst<syst_list.size(); isyst++)
				{
					TString hname = var_list[ivar];
					// if( !f_input->GetListOfKeys()->Contains(hname) ) {cout<<FRED("Error ! Leaf "<<hname<<" not found in TTree!")<<endl;}
					TH1F* h_tmp = (TH1F*) f_input->Get(hname);

					if(!binning_choice) //Different binning if summing all bkgs (TMVA-like histos)
					{
						if(var_list[ivar] == "Lep3Pt")						{h_tmp = new TH1F( "","", 20, 0, 80);}
						else if(var_list[ivar] == "dEtaFwdJetBJet")			{h_tmp = new TH1F( "","", 18, 0, 6);}
						else if(var_list[ivar] == "dEtaFwdJet2BJet")		{h_tmp = new TH1F( "","", 18, 0, 7);}
						else if(var_list[ivar] == "dEtaFwdJetClosestLep")	{h_tmp = new TH1F( "","", 18, 0, 6);}
						else if(var_list[ivar] == "nJet25")					{h_tmp = new TH1F( "","", 6, 1.5, 7.5);}
						else if(var_list[ivar] == "minDRll")				{h_tmp = new TH1F( "","", 18, 0, 4);}
						else if(var_list[ivar] == "maxEtaJet25")			{h_tmp = new TH1F( "","", 18, 0, 5);}
						else if(var_list[ivar] == "dPhiHighestPtSSPair")	{h_tmp = new TH1F( "","", 18, 0, 3.2);}
						else if(var_list[ivar] == "nJetEta1")				{h_tmp = new TH1F( "","", 7, 0.5, 7.5);}
						else if(var_list[ivar] == "lepCharge")				{h_tmp = new TH1F( "","", 3, -1.5, 1.5);}
						else {cout<<FRED("Unknown variable name : "<<var_list[ivar]<<"! Continue !")<<endl;}
					}
					else //Different binning if summing all bkgs (TMVA-like histos)
					{
						if(var_list[ivar] == "Lep3Pt")						{h_tmp = new TH1F( "","", 40, 0, 145);}
						else if(var_list[ivar] == "dEtaFwdJetBJet")			{h_tmp = new TH1F( "","", 40, 0, 7);}
						else if(var_list[ivar] == "dEtaFwdJet2BJet")		{h_tmp = new TH1F( "","", 40, 0, 7);}
						else if(var_list[ivar] == "dEtaFwdJetClosestLep")	{h_tmp = new TH1F( "","", 40, 0, 6.5);}
						else if(var_list[ivar] == "nJet25")					{h_tmp = new TH1F( "","", 9, 2, 11);}
						else if(var_list[ivar] == "minDRll")				{h_tmp = new TH1F( "","", 40, 0, 3.7);}
						else if(var_list[ivar] == "maxEtaJet25")			{h_tmp = new TH1F( "","", 40, 0, 4.5);}
						else if(var_list[ivar] == "dPhiHighestPtSSPair")	{h_tmp = new TH1F( "","", 40, 0, 3.5);}
						else if(var_list[ivar] == "nJetEta1")				{h_tmp = new TH1F( "","", 10, 0, 10);}
						else if(var_list[ivar] == "lepCharge")				{h_tmp = new TH1F( "","", 7, -3, 4);}
						else {cout<<FRED("Unknown variable name : "<<var_list[ivar]<<"! Continue !")<<endl;}
					}

					float weight = 0, mc_weight = 0, tmp = 0;
					tree->SetBranchAddress(var_list[ivar], &tmp); //One variable at a time
					tree->SetBranchAddress("weight", &weight); //Each branch has its own "Weight" (different values among branches)
					tree->SetBranchAddress("mc_weight", &mc_weight); //Each branch has its own "Weight" (different values among branches)

					int nentries = tree->GetEntries();
					for(int ientry = 0; ientry<nentries; ientry++)
					{
						weight = 0; mc_weight = 0; tmp = 0;
						tree->GetEntry(ientry); //Read event

						// weight*= mc_weight;

						//Avoid to get overflow because of inappropriate binning --> Put it into last bin instead! (underflow in first bin)
						if(tmp < h_tmp->GetXaxis()->GetXmax() && tmp > h_tmp->GetXaxis()->GetXmin() ) {h_tmp->Fill(tmp, weight);}
						else if(tmp >= h_tmp->GetXaxis()->GetXmax() ) {Fill_Last_Bin_TH1F(h_tmp, weight);} //overflow
						else if(tmp <= h_tmp->GetXaxis()->GetXmin() ) {Fill_First_Bin_TH1F(h_tmp, weight);} //underflow
					}

					TString output_histo_name = var_list[ivar] + "_" + sample_list[isample];
					f_output->cd();
					h_tmp->Write(output_histo_name);

					delete h_tmp; h_tmp = NULL;
				}
			}
		}

		delete tree; tree = NULL;
		delete f_input; f_input = NULL;
	}

	f_output->Close(); f_output = NULL;

	return;
}







//-----------------------------------------------------------------------------
// ########  ########     ###    ##      ##        ######  ########        ########  ##        #######  ########  ######
// ##     ## ##     ##   ## ##   ##  ##  ##       ##    ## ##     ##       ##     ## ##       ##     ##    ##    ##    ##
// ##     ## ##     ##  ##   ##  ##  ##  ##       ##       ##     ##       ##     ## ##       ##     ##    ##    ##
// ##     ## ########  ##     ## ##  ##  ##       ##       ########        ########  ##       ##     ##    ##     ######
// ##     ## ##   ##   ######### ##  ##  ##       ##       ##   ##         ##        ##       ##     ##    ##          ##
// ##     ## ##    ##  ##     ## ##  ##  ##       ##    ## ##    ##        ##        ##       ##     ##    ##    ##    ##
// ########  ##     ## ##     ##  ###  ###         ######  ##     ##       ##        ########  #######     ##     ######
//-----------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


void tHq_analysis::Draw_Control_Plots(TString bkg_choice)
{
	cout<<endl<<BOLD(FYEL("##################################"))<<endl;
	cout<<FYEL("--- Draw Control Plots ---")<<endl;
	cout<<BOLD(FYEL("##################################"))<<endl<<endl;

	TString input_file_name = "outputs/Control_Histograms_tHq.root";
	if(!Check_File_Existence(input_file_name) ) {cout<<FRED("Error ! File "<<input_file_name<<" not found !")<<endl;}
	TFile* f_input = TFile::Open(input_file_name);

	vector<TString> MC_samples_legend; //List the MC samples which are actually used (to get correct legend)

	for(int ivar=0; ivar<var_list.size(); ivar++)
	{
		TCanvas* c1 = new TCanvas("c1","c1", 1000, 800);
		TLegend* qw = new TLegend(.79,.55,0.999,0.999);

		vector<TH1F*> v_MC_histo;
		MC_samples_legend.clear();

		if(bkg_choice.Contains("tt") )
		{
			v_MC_histo.resize(2);
			MC_samples_legend.resize(2);
		}

		for(int ichan=0; ichan<channel_list.size(); ichan++)
		{
			for(int isample=0; isample<sample_list.size(); isample++)
			{
				TString hname = var_list[ivar] + "_" + sample_list[isample];
				if( !f_input->GetListOfKeys()->Contains(hname) ) {cout<<FRED("Error ! Leaf "<<hname<<" not found in TTree!")<<endl;}
				TH1F* h_tmp = (TH1F*) f_input->Get(hname);

				if(!bkg_choice.Contains("tt")) {h_tmp->Scale(1./h_tmp->Integral() );}

				h_tmp->SetLineColor(color_list[isample]);
				h_tmp->SetLineWidth(2);

				if(bkg_choice.Contains("tt") )
				{
					if(sample_list[isample] == "tHq" || sample_list[isample] == "tHW")
					{
						if(!v_MC_histo[0]) v_MC_histo[0] = (TH1F*) h_tmp->Clone();
						else {v_MC_histo[0]->Add( (TH1F*) h_tmp->Clone() );}

						MC_samples_legend[0] = "Signal";
					}
					else if(bkg_choice == "tt" && sample_list[isample] == "TT")
					{
						if(!v_MC_histo[1]) v_MC_histo[1] = (TH1F*) h_tmp->Clone();
						else {v_MC_histo[1]->Add( (TH1F*) h_tmp->Clone() );}

						MC_samples_legend[1] = "ttbar";
					}
					else if(bkg_choice == "ttV" && (sample_list[isample] == "ttZ" || sample_list[isample] == "ttW") )
					{
						if(!v_MC_histo[1]) v_MC_histo[1] = (TH1F*) h_tmp->Clone();
						else {v_MC_histo[1]->Add( (TH1F*) h_tmp->Clone() );}

						MC_samples_legend[1] = "ttV";
					}
				}
				else
				{
					v_MC_histo.push_back( (TH1F*) h_tmp->Clone() );MC_samples_legend.push_back(sample_list[isample]);
				}

				delete h_tmp; h_tmp = NULL;
			}

		}

		if(bkg_choice.Contains("tt") )
		{
			for(int i=0; i<v_MC_histo.size(); i++)
			{
				v_MC_histo[i]->Scale(1./v_MC_histo[i]->Integral() );
			}
		}

		//--- Draw histograms
		double ymax = 0;
		for(int i=0; i<v_MC_histo.size(); i++)
		{
			if(v_MC_histo[i]->GetMaximum() > ymax) {ymax = v_MC_histo[i]->GetMaximum();}
		}
		for(int i=0; i<v_MC_histo.size(); i++)
		{
			v_MC_histo[i]->SetMaximum(ymax * 1.1);
			v_MC_histo[i]->Draw("same HIST");
		}
		for(int i=0; i<MC_samples_legend.size(); i++)
		{
			qw->AddEntry(v_MC_histo[i], MC_samples_legend[i], "L");
		}

		qw->Draw("same");

		//--- Format Canvas
		c1->SetTitle(var_list[ivar]);

		TString output_name = "plots/" + var_list[ivar] + plot_extension;
		c1->SaveAs(output_name);

		for(int i=0; i<v_MC_histo.size(); i++)
		{
			delete v_MC_histo[i]; v_MC_histo[i] = NULL;
		}

		delete qw; qw = NULL;
		delete c1; c1 = NULL;
	}

	return;
}
