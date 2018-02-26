#include "../include/Hist.h"
//#include "../include/TTbarDileptonAnalysis.h"
#include "../include/TTbarHiggsMultileptonAnalysis.h"
//#include "../include/TTbarHiggsTFAnalysis.h"
#include "../include/TTbarHiggsBTagEff.h"

#include "TChain.h"

#include <iostream>

int main(int argc, char *argv[])
{
    if( argc < 3 )                                                                                                          
    {
        std::cout << "NtupleAnalyzer usage:"         	  	   << std::endl;
        std::cout << "--file: input filename"        	  	   << std::endl;
        std::cout << "--tree: TTree name"            	  	   << std::endl; 
        std::cout << "--outfile: output path"             	   << std::endl;
        std::cout << "--nmax: max number of events"       	   << std::endl;
        std::cout << "--nowe: number of weighted events"  	   << std::endl;
        std::cout << "--xsec: cross section (pb)"    	  	   << std::endl;
        std::cout << "--lumi: luminosity in (pb-1)"  	  	   << std::endl;
        std::cout << "--isdata: is data or not"      	  	   << std::endl;
	std::cout << "--doSystCombine: do systematics for combine" << std::endl;
        std::cout << "--dataset: dataset name"      	           << std::endl;
        exit(1);
    }

    const char *fname_str  = "list.txt";
    const char *stream_str = "Nt";
    int nmax = -1;
    int nowe = -1;
    float lumi = -1.;
    float xsec = -1.;
    bool isdata = 0;
    bool doSystCombine = 0;
    const char *outfile_str = "./";
    const char *dataset_str = "";
 
    for(int i=0;i<argc;i++)
    {
        if( ! strcmp(argv[i],"--file") )          fname_str  = argv[i+1];
        if( ! strcmp(argv[i],"--tree") )          stream_str = argv[i+1];
        if( ! strcmp(argv[i],"--outfile") )       outfile_str = argv[i+1];
        if( ! strcmp(argv[i],"--nmax") )          nmax = atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--nowe") )          nowe = atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--lumi") )          lumi = atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--xsec") )          xsec = atof(argv[i+1]);
        if( ! strcmp(argv[i],"--isdata") )        isdata = (bool) atoi(argv[i+1]);
	if( ! strcmp(argv[i],"--doSystCombine") ) doSystCombine = (bool) atoi(argv[i+1]);
        if( ! strcmp(argv[i],"--dataset") )       dataset_str = argv[i+1];
    }

    const char *fname = fname_str;                                                                                         
    const char *stream = stream_str;                                                                                        
    const char *outfile = outfile_str;                                                                                        
    const char *dataset = dataset_str;                                                                                        

    std::cout << "--file=" << fname  << std::endl;                                                                          
    std::cout << "--tree=" << stream << std::endl;                                                                          
    std::cout << "--nmax=" << nmax   << std::endl;
    std::cout << "--nowe=" << nowe   << std::endl;
    std::cout << "--xsec=" << xsec   << std::endl;
    std::cout << "--lumi=" << lumi   << std::endl;
    std::cout << "--outfile=" << outfile << std::endl;
    std::cout << "--isdata="  << isdata  << std::endl;
    std::cout << "--doSystCombine=" << doSystCombine  << std::endl;
    std::cout << "--dataset=" << dataset  << std::endl;

    //TChain *thetree = 0;

    // TTH Transfer Function analysis

    // TTbarHiggsTFAnalysis *TTHTFanalysis = new TTbarHiggsTFAnalysis(fname,thetree,"TTbarHiggs",stream);
    // TTHTFanalysis->createHistograms();
    // TTHTFanalysis->Loop();
    // TTHTFanalysis->writeHistograms();

    // TTH Trigger Studies

    // TTbar...

    // TTH Charge Flip Determination

    // TTbar...

    // TTH Fake Rate Determination

    // TTbar...
  
      
    //tHq analysis
    TTbarHiggsMultileptonAnalysis* tHqAnalysis = new TTbarHiggsMultileptonAnalysis(fname,dataset,stream,outfile,isdata,doSystCombine,xsec,lumi,nowe,nmax);
    
    //tHqAnalysis->createHistograms();
    
    tHqAnalysis->Loop();
    
    //tHqAnalysis->writeHistograms();
    
    delete tHqAnalysis;
}
