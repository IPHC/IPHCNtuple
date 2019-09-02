#ifndef BASE_H
#define BASE_H

#include "Helper.h"

extern unsigned int idx;

class Base : public Helper
{
 public:
   Base();
   virtual ~Base();
   
   template< class BR >
     
     bool CHECK(BR br)
       {
	  bool res = false;
	  if( br )
	    {
	       if( idx < br->size() )
		 res = true;
	    }
	  
	  return res;
       }
   
   float pt;
   float conept;
   float eta;
   float phi;
   float E;
   float m;
   int   iElec;
   int   iMuon;
   int   iTau;
   bool  isFakeableTTH;
   bool  isLooseTTH;
   bool  isMediumTTH;
   bool  isTightTTH;
   bool  hasMCMatch;
   bool  hasChargeMCMatch;
   bool  hasPhotonMCMatch; //for ele only
   int   charge;
   bool  tightCharge;
   float lepMVA;
   int matchedJetId;

   float gen_pt;
   float gen_eta;
   float gen_phi;
   float gen_m;
   float gen_E;
   int gen_status;
   int gen_id;
   float gen_dr;
   
   ClassDef(Base,1)
};

struct sort_by_pt
{   
   bool operator () (const Base& lhs, const Base& rhs)
     {	
	return lhs.pt > rhs.pt;
     }
};

struct sort_by_conept
{   
   bool operator () (const Base& lhs, const Base& rhs)
     {	
	return lhs.conept > rhs.conept;
     }
};

#endif
