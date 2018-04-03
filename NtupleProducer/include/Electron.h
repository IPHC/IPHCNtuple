#ifndef ELECTRON_H
#define ELECTRON_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Electron : public Base
{
    public:

        Electron();
        virtual ~Electron();

        static bool sortPtPredicate(Electron lhs, Electron rhs) {return (lhs.pt() > rhs.pt());};

        int ID()                             {return _ID;};

        void setFakeType(int faketype)       {_fakeType = faketype;};
        int fakeType()                       {return _fakeType;};
        bool sel();

        // General informations
        float E()                            {return _E;};
        float pt()                           {return _pt;};
        float ptCor()                        {return _ptCor;};
        float ptUnc()                        {return _ptUnc;};
        float eta()                          {return _eta;};
        float phi()                          {return _phi;};
        float m()                            {return _m;};
        int   charge()                       {return _charge;};
        int   id()                           {return _id;};

        // Id
        bool isLoose()                       {return _isLoose;};
        bool isTight()                       {return _isTight;};
        bool isLooseMVA()                    {return _isLooseMVA;};
        bool isTightMVA()                    {return _isTightMVA;};
        bool isLooseTTH()                    {return _isLooseTTH;};
        bool isFakeableTTH()                 {return _isFakeableTTH;};
        bool isTightTTH()                    {return _isTightTTH;};

        // Variables for Id
        float dxy()                             {return _dxy;};
        float dz()                              {return _dz;};
        bool  passPtEta()                       {return _passPtEta;};
        int   tightCharge()                     {return _tightCharge;};
        bool passTightCharge()                  {return _passTightCharge;};
        bool cutEventSel()                      {return _cutEventSel;};
        bool noLostHits()                       {return _noLostHits;};

        // MVA
        float lepMVA()                          {return _lepMVA;};
        float lepMVA_TTH()	                    {return _lepMVA_TTH;};

        float lepMVA_miniRelIsoCharged()        {return _lepMVA_miniRelIsoCharged;};
        float lepMVA_miniRelIsoNeutral()        {return _lepMVA_miniRelIsoNeutral;};
        float lepMVA_jetPtRelv2()               {return _lepMVA_jetPtRelv2;};
        //float lepMVA_jetDR()                 {return _lepMVA_jetDR;};
        float lepMVA_jetPtRatio()               {return _lepMVA_jetPtRatio;};
        float lepMVA_jetBTagCSV()               {return _lepMVA_jetBTagCSV;};
        float lepMVA_sip3d()                    {return _lepMVA_sip3d;};
        float lepMVA_dxy()                      {return _lepMVA_dxy;};
        float lepMVA_dz()                       {return _lepMVA_dz;};
        float lepMVA_mvaId()                    {return _lepMVA_mvaId;};
        float lepMVA_eta()		                {return _lepMVA_eta;};
        float lepMVA_jetNDauChargedMVASel()     {return _lepMVA_jetNDauChargedMVASel;};

        // Other variables
        bool passChargeFlip()                   {return _passChargeFlip;};
        bool hasMatchedConversion()             {return _hasMatchedConversion;};

        bool isGsfCtfScPixChargeConsistent()    {return _isGsfCtfScPixChargeConsistent;};

        //float dB3D()                           {return _dB3D;};
        //float edB3D()                          {return _edB3D;};

        float miniIso()                        {return _miniIso;};
        int   nlosthits()                      {return _nlosthits;};
        float sigmaIetaIeta()                  {return _sigmaIetaIeta;};
        float hadronicOverEm()                 {return _hadronicOverEm;};
        float deltaEtaSuperClusterTrackAtVtx() {return _deltaEtaSuperClusterTrackAtVtx;};
        float deltaPhiSuperClusterTrackAtVtx() {return _deltaPhiSuperClusterTrackAtVtx;};
        float see()                            {return _see;};
        float superCluster_eta()               {return _superCluster_eta;};
        float correctedEcalEnergy()            {return _correctedEcalEnergy;};
        float ecalEnergy()                     {return _ecalEnergy;};
        float eSuperClusterOverP()             {return _eSuperClusterOverP;};
        float trackMomentumError()             {return _trackMomentumError;};
        float mvaNonTrigV0()                   {return _mvaNonTrigV0;};

        bool  passCV()                         {return _passCV;};
        float ip3d()                           {return _ip3d;};
        float ip3dErr()                        {return _ip3dErr;};

        void read();
        void init();

        //NEW
        bool passMuOverlap() {return _passMuOverlap;}
        bool passConditions() {return _passConditions;}
        float eInvMinusPinv() {return _eInvMinusPinv;}


    protected:

        int     _ID;

        int     _fakeType;

        // General informations
        float   _E;
        float   _pt;
        float   _ptCor;
        float   _ptUnc;
        float   _eta;
        float   _phi;
        float   _m;
        int     _charge;
        int     _id;

        // Id
        bool    _isLooseCBId;
        bool    _isMediumCBId;
        bool    _isLoose;
        bool    _isMedium;
        bool    _isTight;
        bool    _isLooseMVA;
        bool    _isTightMVA;
        bool    _isLooseTTH;
        bool    _isFakeableTTH;
        bool    _isTightTTH;

        // Variables for Id
        float   _dxy;
        float   _dz;
        float   _miniIso;
        int     _nlosthits;
        bool    _passCV;
        bool    _isPCC;
        bool    _passPtEta;
        float   _ip3d;
        float   _ip3dErr;
        int     _tightCharge;
        bool    _passTightCharge;
        bool    _cutEventSel;
        bool    _noLostHits;

        // MVA
        float _lepMVA;
        float _lepMVA_TTH;

        float _lepMVA_miniRelIsoCharged;
        float _lepMVA_miniRelIsoNeutral;
        float _lepMVA_jetPtRelv2;
        //float _lepMVA_jetDR;
        float _lepMVA_jetPtRatio;
        float _lepMVA_jetBTagCSV;
        float _lepMVA_sip3d;
        float _lepMVA_dxy;
        float _lepMVA_dz;
        float _lepMVA_mvaId;

        float _lepMVA_eta;
        float _lepMVA_jetNDauChargedMVASel;

        // Other variables
        bool _passChargeFlip;
        bool _hasMatchedConversion;
        bool _isGsfCtfScPixChargeConsistent;

        //float _dB3D;
        //float _edB3D;

        float _sigmaIetaIeta;
        float _hadronicOverEm;
        float _correctedEcalEnergy;
        float _ecalEnergy;
        float _eSuperClusterOverP;
        float _deltaEtaSuperClusterTrackAtVtx;
        float _deltaPhiSuperClusterTrackAtVtx;
        float _see;
        float _superCluster_eta;

        float _trackMomentumError;
        float _mvaNonTrigV0;

        //NEw
        bool _passMuOverlap;
        bool _passConditions;
        float _eInvMinusPinv;

        ClassDef(Electron,1)
};

#endif
