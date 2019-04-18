#ifndef LEPTON_H
#define LEPTON_H

#include "Electron.h"
#include "Muon.h"

#include "TLorentzVector.h"

class Lepton
{

    public:

        Lepton();
        virtual ~Lepton();
	
	float           pt;
	float           conept;
        //float           ptCor;
        //float           ptUnc;
        float           eta;
        float           phi;
        float           E;

        TLorentzVector  p4;
	
	int id;
        int             idx;

        bool            isFakeableTTH;
	bool            isLooseTTH;
        bool            isTightTTH;
        float           lepMVA;

        bool            cutEventSel;
	bool            passTightCharge;
        bool            noLostHits;

        int             charge;
	bool		hasMCMatch;
	bool		hasChargeMCMatch;
	bool		hasPhotonMCMatch;
	

        template <class T> void setLepton(T *lep, int idx, bool isE, bool isMu)
        {
            pt                 = lep->pt;
	    conept                 = lep->conept;
            //ptCor              = lep->ptCor();
            //ptUnc              = lep->ptUnc();
            eta                = lep->eta;
            phi                = lep->phi;
            E                  = lep->E;

            //p4.SetPtEtaPhiE(ptUnc,eta,phi,E);
	    p4.SetPtEtaPhiE(pt,eta,phi,E); //FIXME -- ptUnc ?

            idx                = idx;

            charge             = lep->charge;

    	    isFakeableTTH      = lep->isFakeableTTH;
	    isLooseTTH = lep->isLooseTTH;
	    
            //if(isE || isMu) {_isLooseTTH = lep->isLooseTTH;}
	    //else {_isLooseTTH = true;} //isLooseTTH not implemented for tau
	    
	    passTightCharge    = lep->tightCharge;
	    isTightTTH         = lep->isTightTTH;
            lepMVA         = lep->lepMVA;
	    
	    id = lep->id;
	    
	    hasMCMatch = lep->hasMCMatch;
	    hasChargeMCMatch = lep->hasChargeMCMatch;
	    if(isE) {hasPhotonMCMatch = lep->hasPhotonMCMatch;} //For ele only
        }
	
};

#endif
