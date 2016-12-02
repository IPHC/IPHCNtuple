#include "include/Base.h"
#include "include/NtupleProducer.h"

ClassImp(Base)

Base::Base()
{
}

Base::~Base()
{
}

float Base::GetDPhi(float phi1,float phi2)
{
   float deltaPhi=fabs(phi1-phi2);
   if (deltaPhi>M_PI) deltaPhi=2*M_PI-deltaPhi;
   return deltaPhi;
}

float Base::GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
      if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
