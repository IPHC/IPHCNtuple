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
   
   float _pt;
   float _eta;
   float _phi;
   float _E;
   int   iElec;
   int   iMuon;
   int   iTau;
   bool  _isTightTTH;
   bool  _hasMCMatch;
   int   _charge;
   bool  _tightCharge;
   float _lepMVA;
   
   ClassDef(Base,1)
};

struct sort_by_pt
{   
   bool operator () (const Base& lhs, const Base& rhs)
     {	
	return lhs._pt > rhs._pt;
     }
};

#endif
