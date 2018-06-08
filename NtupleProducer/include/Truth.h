#ifndef TRUTH_H
#define TRUTH_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Truth : public Base
{
 public:
   
   int mc_truth_n;
   std::vector<int>    mc_truth_id;
   std::vector<int>    mc_truth_label;
   std::vector<float>  mc_truth_pt;
   std::vector<float>  mc_truth_eta;
   std::vector<float>  mc_truth_phi;
   std::vector<float>  mc_truth_E;
   
   int gen_n;
   float gen_PVz;
   std::vector<float>  gen_pt;
   std::vector<float>  gen_eta;
   std::vector<float>  gen_phi;
   std::vector<float>  gen_m;
   std::vector<int>    gen_id;
   std::vector<int>    gen_status;
   std::vector<int>    gen_mother_id;
   
   Truth();
   virtual ~Truth();
   
   ClassDef(Truth,1)
};

#endif
