#include "include/NtupleProducer.h"
#include <iostream>
#include <iomanip>

ClassImp(TriggerObjExt)

TriggerObjExt::TriggerObjExt()
{
}

TriggerObjExt::~TriggerObjExt()
{
}

void TriggerObjExt::read()
{
   
   // 
   ID = idx;
   
   pT  = ntP->triggerobject_pt->at(ID);
   eta = ntP->triggerobject_eta->at(ID); 
   phi = ntP->triggerobject_phi->at(ID);
   
   collection = ntP->triggerobject_collection->at(ID);
   
   //std::cout <<"pt, eta, phi, collection "<< std::endl;          
   //std::cout <<_pT<<" "<<_eta<<" "<<_phi<<" "<<_collection<< std::endl;
   
   //
   filterIds_n    =  ntP->triggerobject_filterIds_n->at(ID);
   filterLabels_n =  ntP->triggerobject_filterLabels_n->at(ID);
   pathNamesAll_n =  ntP->triggerobject_pathNamesAll_n->at(ID);
   
   //std::cout <<"Number of filters, labels, paths per trigger object "<< std::endl;
   //std::cout <<_filterIds_n  <<" "<< _filterLabels_n <<" "<< _pathNamesAll_n << std::endl;
   
   pathNamesAll_offset = 0;
   filterLabels_offset = 0;
   filterIds_offset    = 0;
  
   for (int i=0; i<ID; i++)
     {
	pathNamesAll_offset += ntP->triggerobject_pathNamesAll_n->at(i);
	filterIds_offset    += ntP->triggerobject_filterIds_n->at(i);
	filterLabels_offset += ntP->triggerobject_filterLabels_n->at(i);
     }   
     
   //std::cout <<"offsets "<< std::endl;
   //std::cout <<_filterIds_offset  <<" "<< _filterLabels_offset <<" "<< _pathNamesAll_offset << std::endl;
     
   //    
   std::vector<string> pathNamesAll      = *ntP->triggerobject_pathNamesAll;
   std::vector<bool> pathNamesAll_isL3   = *ntP->triggerobject_pathNamesAll_isL3;
   std::vector<bool> pathNamesAll_isLF   = *ntP->triggerobject_pathNamesAll_isLF;
   std::vector<bool> pathNamesAll_isBoth = *ntP->triggerobject_pathNamesAll_isBoth;
   std::vector<bool> pathNamesAll_isNone = *ntP->triggerobject_pathNamesAll_isNone;
   std::vector<int> filterIds            = *ntP->triggerobject_filterIds;
   std::vector<std::string> filterLabels = *ntP->triggerobject_filterLabels;

   
   pathNamesAll = std::vector<std::string>(pathNamesAll.begin() + pathNamesAll_offset, pathNamesAll.begin() + pathNamesAll_offset + pathNamesAll_n);
   
   pathNamesAll_isL3   = std::vector<bool>(pathNamesAll_isL3.begin() + pathNamesAll_offset, pathNamesAll_isL3.begin() + pathNamesAll_offset + pathNamesAll_n);
   pathNamesAll_isLF   = std::vector<bool>(pathNamesAll_isLF.begin() + pathNamesAll_offset, pathNamesAll_isLF.begin() + pathNamesAll_offset + pathNamesAll_n);
   pathNamesAll_isBoth = std::vector<bool>(pathNamesAll_isBoth.begin() + pathNamesAll_offset, pathNamesAll_isBoth.begin() + pathNamesAll_offset + pathNamesAll_n);
   pathNamesAll_isNone = std::vector<bool>(pathNamesAll_isNone.begin() + pathNamesAll_offset, pathNamesAll_isNone.begin() + pathNamesAll_offset + pathNamesAll_n);
    
   filterIds    = std::vector<int>(filterIds.begin()+filterIds_offset, filterIds.begin()+filterIds_offset+filterIds_n);
   filterLabels = std::vector<std::string>(filterLabels.begin()+filterLabels_offset, filterLabels.begin()+filterLabels_offset+filterLabels_n);

    
   /*
   
   std::cout<<"========= paths size " <<  _pathNamesAll.size() << std::endl;  
     
   for(int j=0;j<_pathNamesAll.size();j++)
   { 
    std::cout<<_pathNamesAll.at(j)<<" "<< _pathNamesAll_isL3.at(j)<<" "<< 
   	       _pathNamesAll_isLF.at(j)<<" "<<_pathNamesAll_isBoth.at(j)<<" "<< _pathNamesAll_isNone.at(j)  <<std::endl;}
  
   std::cout<<"========= filters size"<< std::endl;
   for(int j=0;j<_filterIds.size();j++)
   { 
    std::cout<< _filterIds.at(j)  <<std::endl;}
   
   std::cout<<"========= labels size "<< std::endl;
   for(int j=0;j<_filterLabels.size();j++)
   { 
    std::cout<< _filterLabels.at(j)  <<std::endl;}
   
   */     
}

void TriggerObjExt::init()
{  
   
   pT  = -999.; 
   eta = -999.; 
   phi = -999.; 
   collection = "empty";
   
   pathNamesAll_n = 0; 
   filterLabels_n = 0;
   filterIds_n = 0; 
   
   pathNamesAll_offset = 0; //to decode flattrees
   filterLabels_offset = 0; 
   filterIds_offset = 0; 
   
   filterIds.clear();
   filterLabels.clear();
    
   pathNamesAll.clear();    
   pathNamesAll_isL3.clear();
   pathNamesAll_isLF.clear();
   pathNamesAll_isBoth.clear();
   pathNamesAll_isNone.clear();   	
}

bool TriggerObjExt::sel()
{  
   bool isTrigger = false;
      
   //if (std::find(_filterIds.begin(),_filterIds.end(), 82)!= _filterIds.end() || std::find(_filterIds.begin(),_filterIds.end(), 83) != _filterIds.end()) isTrigger = true;
    
   for(int j=0;j<pathNamesAll_n;j++)
     { 
	//SL
	std::size_t ok0 = pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL");
	std::size_t ok1 = pathNamesAll.at(j).find("HLT_IsoMu20");
	std::size_t ok2 = pathNamesAll.at(j).find("HLT_IsoTkMu20");
	std::size_t ok3 = pathNamesAll.at(j).find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL");
	std::size_t ok4 = pathNamesAll.at(j).find("HLT_Ele23_WPLoose_Gsf");
	
	//Tri-lep
	std::size_t ok5 = pathNamesAll.at(j).find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL");
	std::size_t ok6 = pathNamesAll.at(j).find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL");
	std::size_t ok7 = pathNamesAll.at(j).find("HLT_TripleMu_12_10_5");
	std::size_t ok8 = pathNamesAll.at(j).find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL");
    
	//Di-lep
	std::size_t ok9  = pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");
	std::size_t ok10 = pathNamesAll.at(j).find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL");
	std::size_t ok11 = pathNamesAll.at(j).find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL");
	std::size_t ok12 = pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL");
	std::size_t ok13 = pathNamesAll.at(j).find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL");
    
	std::size_t ok22 = pathNamesAll.at(j).find("HLT_Ele22_eta2p1_WPLoose_Gsf");
	
	if ( ok0!=std::string::npos  || ok1!=std::string::npos  || ok2!=std::string::npos  || ok3!=std::string::npos  || 
	     ok4!=std::string::npos  ||
	     ok5!=std::string::npos  || ok6!=std::string::npos  || ok7!=std::string::npos  || ok8!=std::string::npos  ||
	     ok9!=std::string::npos  || ok10!=std::string::npos || ok11!=std::string::npos || ok12!=std::string::npos ||
	     ok13!=std::string::npos || ok22!=std::string::npos 
	   ) 
	  {
	     isTrigger = true;
	     break;
	  }
     }
   
   for(int j=0;j<filterLabels_n;j++)
     { 
	std::size_t ok14 = filterLabels.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
	std::size_t ok15 = filterLabels.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
	std::size_t ok16 = filterLabels.at(j).find("hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17");
	std::size_t ok17 = filterLabels.at(j).find("hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8");
	std::size_t ok18 = filterLabels.at(j).find("hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17");
	std::size_t ok19 = filterLabels.at(j).find("hltDiMuonGlbFiltered17TrkFiltered8");
	std::size_t ok20 = filterLabels.at(j).find("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4");
	std::size_t ok21 = filterLabels.at(j).find("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4");
	std::size_t ok23 = filterLabels.at(j).find("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2");
	std::size_t ok24 = filterLabels.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
	std::size_t ok25 = filterLabels.at(j).find("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
   
	if ( ok14!=std::string::npos || ok15!=std::string::npos || ok16!=std::string::npos ||
	     ok17!=std::string::npos || ok18!=std::string::npos || ok19!=std::string::npos || ok20!=std::string::npos ||
	     ok21!=std::string::npos || ok23!=std::string::npos || ok24!=std::string::npos || ok25!=std::string::npos
	   ) 
	  {
	     isTrigger = true;
	     break;
	  }	
     }          
 
   //std::cout <<"isTrigger " << isTrigger << std::endl;
   
   return isTrigger;
}
