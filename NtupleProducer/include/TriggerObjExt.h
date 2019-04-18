#ifndef TriggerObjExt_H
#define TriggerObjExt_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "TriggerObj.h"

class TriggerObjExt : public TriggerObj
{
 public:

   TriggerObjExt();
   virtual ~TriggerObjExt();
   
   void read();
   void init();
   bool sel();

 protected:
   
   int ID;
   
   float pT;
   float eta; 
   float phi;
   std::string collection;
   
   int filterIds_n;
   std::vector<int> filterIds;
   
   int filterLabels_n;
   std::vector<std::string> filterLabels;
   
   int pathNamesAll_n;
   std::vector<std::string> pathNamesAll;
   std::vector<bool> pathNamesAll_isL3;
   std::vector<bool> pathNamesAll_isLF;
   std::vector<bool> pathNamesAll_isBoth;
   std::vector<bool> pathNamesAll_isNone;
   
   int pathNamesAll_offset;
   int filterLabels_offset; 
   int filterIds_offset; 
   
   ClassDef(TriggerObjExt,1)
};

#endif
