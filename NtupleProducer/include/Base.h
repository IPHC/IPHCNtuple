#ifndef BASE_H
#define BASE_H

#include "TObject.h"
#include "TLorentzVector.h"

extern unsigned int idx;

class Base : public TObject
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
   
   float GetDPhi(float phi1,float phi2);
   float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
   
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
   int   charge;
   bool  tightCharge;
   float lepMVA;
   
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
