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
        bool            noLostHits;

        int             charge;
	
	int Get_LeptonID(int iElec, int iMuon, int iTau, int charge)
	{
		int id = 0;
		if(iElec>=0) {id=11;}
		else if(iMuon>=0) {id=13;}
		else if(iTau>=0) {id=15;}	
		
		id*= charge;
		
		return id;
	}
	
	
        template <class T> void setLepton(T *lep, int idx, bool isE, bool isMu)
        {
            pt                 = lep->pt;
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
	    isTightTTH         = lep->isTightTTH;
            lepMVA         = lep->lepMVA;
	    
	    id = Get_LeptonID(lep->iElec, lep->iMuon, lep->iTau, lep->charge);


            //cutEventSel        = lep->cutEventSel();       // last set of cuts for electron, used at event selection only
            //noLostHits         = lep->noLostHits();
        }
	
};

#endif
