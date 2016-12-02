#ifndef Tau_H
#define Tau_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Tau : public Base
{
    public:

        Tau();
        virtual ~Tau();

        static bool sortPtPredicate(Tau lhs, Tau rhs)   {return (lhs.pt() > rhs.pt());};

        int ID()                                        {return _ID;};

        void setFakeType(int faketype)                  {_fakeType = faketype;};
        int  fakeType()                                 {return _fakeType;};

        bool sel();
        void read();
        void init();

        // kinematics
        float E()               {return _E;};
        float pt()              {return _pt;};
        float ptCor()           {return _ptCor;};
        float ptUnc()           {return _ptUnc;};
        float eta()             {return _eta;};
        float phi()             {return _phi;};
        float m()               {return _m;};
        int   charge()          {return _charge;};
        int   id()              {return _id;};

        // Id
        bool isLoose()          {return _isLoose;};
        bool isFakeableTTH()    {return _isFakeableTTH;};
        bool isTightTTH()       {return _isTightTTH;};
        float lepMVA_TTH()      {return _lepMVA_TTH;};
        bool passTightCharge()  {return _passTightCharge;};
        bool cutEventSel()      {return _cutEventSel;};
        bool noLostHits()       {return _noLostHits;};

        // Variables for Id
        float dxy()             {return _dxy;};
        float dz()              {return _dz;};

        // Other variables

    protected:

        int   _ID; //idx in flatree

        int   _fakeType;

        // General informations
        float _E;
        float _pt;
        float _ptCor;
        float _ptUnc;
        float _eta;
        float _phi;
        float _m;
        int   _charge;
        int   _id; //pdgid

        float _dxy;
        float _dz;

        bool  _isLoose;
        bool  _isFakeableTTH;
        bool  _isTightTTH;
        float _lepMVA_TTH;
        bool _passTightCharge;
        bool _cutEventSel;
        bool _noLostHits;

        //-------------------------------
        // http://kskovpen.web.cern.ch/kskovpen/IPHCFlatTree/table_MantaRay-patch7_20150829.html
        //-------------------------------

        int   _decayMode;
        bool  _hasLeadChargedHadrCand;
        float _leadingTrackPt;
        float _leadingTrackCharge;

        //-------------------------------
        // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        //-------------------------------


        float _decayModeFindingOldDMs;
        float _decayModeFindingNewDMs ;//(only for analysis involving high pT taus, like !SUSY searches)

        // isolation discriminators

        float _byLooseCombinedIsolationDeltaBetaCorr3Hits;
        float _byMediumCombinedIsolationDeltaBetaCorr3Hits;
        float _byTightCombinedIsolationDeltaBetaCorr3Hits;

        float _byLooseIsolationMVA3newDMwLT;  //placeholder, will it exist for Run II ??
        float _byMediumIsolationMVA3newDMwLT; //placeholder, will it exist for Run II ??
        float _byTightIsolationMVA3newDMwLT;  //placeholder, will it exist for Run II ??

        // isolation discriminators

        float _byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
        float _byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
        float _byTightCombinedIsolationDeltaBetaCorr3HitsdR03;

        float _byLooseIsolationMVArun2v1DBdR03oldDMwLT;
        float _byMediumIsolationMVArun2v1DBdR03oldDMwLT;
        float _byTightIsolationMVArun2v1DBdR03oldDMwLT;
        float _byVTightIsolationMVArun2v1DBdR03oldDMwLT;

        // raw values of the isolation

        float _byCombinedIsolationDeltaBetaCorrRaw3Hits;
        float _chargedIsoPtSum;
        float _neutralIsoPtSum;
        float _puCorrPtSum;

        // muon discriminators

        float _againstMuonLoose3;
        float _againstMuonTight3;

        // electron discriminators

        //AC8X
        float _againstElectronVLoose;
        float _againstElectronLoose;
        float _againstElectronMedium;
        float _againstElectronTight;

        float _pfEssential_jet_pt;
        float _pfEssential_jet_eta;
        float _pfEssential_jet_phi;
        float _pfEssential_jet_m;
        float _pfEssential_jetCorr_pt;
        float _pfEssential_jetCorr_eta;
        float _pfEssential_jetCorr_phi;
        float _pfEssential_jetCorr_m;
        float _pfEssential_hasSV;
        float _pfEssential_sv_x;
        float _pfEssential_sv_y;
        float _pfEssential_sv_z;
        float _pfEssential_flightLengthSig;
        float _pfEssential_dxy;
        float _pfEssential_dxy_error;
        float _pfEssential_dxy_Sig;

        ClassDef(Tau,1)
};

#endif
