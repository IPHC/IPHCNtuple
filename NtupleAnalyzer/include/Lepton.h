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

        float pt()                  {return _pt;};
        float ptCor()               {return _ptCor;};
        float ptUnc()               {return _ptUnc;};
        float eta()                 {return _eta;};
        float phi()                 {return _phi;};
        float E()                   {return _E;};

        TLorentzVector p4()         {return _p4;};

        int     id()                {return _id;};

        int     idx()               {return _idx;};
        bool    isElectron()        {return _isElectron;};
        bool    isMuon()            {return _isMuon;};
        float   lepMVA_TTH()        {return _lepMVA_TTH;};

        bool    passTightCharge()   {return _passTightCharge;};
        bool    cutEventSel()       {return _cutEventSel;};
        bool    noLostHits()        {return _noLostHits;};

        bool isFakeableTTH()        {return _isFakeableTTH;};
	bool isLooseTTH()           {return _isLooseTTH;};
        bool isTightTTH()           {return _isTightTTH;};

        int charge()                {return _charge;};

        template <class T> void setLepton(T *lep, int idx, bool isE, bool isMu)
        {
            _pt                 = lep->pt();
            _ptCor              = lep->ptCor();
            _ptUnc              = lep->ptUnc();
            _eta                = lep->eta();
            _phi                = lep->phi();
            _E                  = lep->E();

            _p4.SetPtEtaPhiE(_ptUnc,_eta,_phi,_E);

            _id                 = lep->id();

            _idx                = idx;
            _isElectron         = isE;
    	    _isMuon             = isMu;

            _charge             = lep->charge();

    	    _isFakeableTTH      = lep->isFakeableTTH();
            if(isE || isMu) {_isLooseTTH = lep->isLooseTTH();}
	    else {_isLooseTTH = true;} //isLooseTTH not implemented for tau
	    _isTightTTH         = lep->isTightTTH();
            _lepMVA_TTH         = lep->lepMVA_TTH();

            _passTightCharge    = lep->passTightCharge();
            _cutEventSel        = lep->cutEventSel();       // last set of cuts for electron, used at event selection only
            _noLostHits         = lep->noLostHits();
        }

    protected:

        float           _pt;
        float           _ptCor;
        float           _ptUnc;
        float           _eta;
        float           _phi;
        float           _E;

        int             _id;

        TLorentzVector  _p4;

        int             _idx;
        bool            _isElectron;
        bool            _isMuon;

        bool            _isFakeableTTH;
	bool            _isLooseTTH;
        bool            _isTightTTH;
        float           _lepMVA_TTH;

        bool            _passTightCharge;
        bool            _cutEventSel;
        bool            _noLostHits;

        int             _charge;
};

#endif
