#ifndef JET_H
#define JET_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Jet : public Base
{
    public:

        Jet();
        virtual ~Jet();

        static bool sortPtPredicate(Jet lhs, Jet rhs) {return (lhs.pt() > rhs.pt());};

        int   ID()       {return _ID;};

        bool  sel();
		void  sel_tZq(std::vector<Electron> &theelec, std::vector<Muon> &themuons);

        // kinematics
        float E()        		{return _E;};
        float pt()       		{return _pt;};
        float eta()      		{return _eta;};
        float phi()      		{return _phi;};
        float m()        		{return _m;};
		float pt_uncorrected()  {return _pt_uncorrected;};
		float uncert() 			{return _uncert;};
		bool passJetID() 		{return _passJetID;};
		bool passJetCleaning()	{return _passJetCleaning;};

        bool  isTight()  {return _isTight;};
        bool  isLoose()  {return _isLoose;};

        int   ntrk()     {return _ntrk;};

        float CSVv2()    {return _CSVv2;};

        float jet_partonFlavour() {return _jet_partonFlavour  ;};
        float jet_hadronFlavour() {return _jet_hadronFlavour  ;};

        float jet_genJet_pt()     {return  _jet_genJet_pt    ;};
		float jet_genJet_E()      {return  _jet_genJet_E   ;};
        /*float jet_genJet_eta()    {return  _jet_genJet_eta ;};
          float jet_genJet_phi()    {return  _jet_genJet_phi ;};
          float jet_genJet_m()      {return  _jet_genJet_m   ;};
          float jet_genJet_E()      {return  _jet_genJet_E   ;};
          float jet_genJet_status() {return  _jet_genJet_status ;};
          float jet_genJet_id()     {return  _jet_genJet_id  ;};*/

        float jet_genParton_pt()     {return  _jet_genParton_pt     ;};
		float jet_genParton_E()      {return  _jet_genParton_E      ;};
        /*float jet_genParton_eta()    {return  _jet_genParton_eta    ;};
          float jet_genParton_phi()    {return  _jet_genParton_phi    ;};
          float jet_genParton_m()      {return  _jet_genParton_m      ;};
          float jet_genParton_E()      {return  _jet_genParton_E      ;};
          float jet_genParton_status() {return  _jet_genParton_status ;};*/
          float jet_genParton_id()     {return  _jet_genParton_id     ;};

		float JES_uncert()	{return _JES_uncert;}
		float pt_JER()		{return _pt_JER;}
		float pt_JER_down()	{return _pt_JER_down;}
		float pt_JER_up()	{return _pt_JER_up;}

        void  read(bool isdata);
        void  init();
	void  setJESUncertainty(float JES_uncert);
        void  JECUncertainty();

    protected:

        int _ID;

        float _E;
        float _pt;
        float _eta;
        float _phi;
        float _m;

		float _uncert;

   		float _pt_uncorrected;
   		bool  _passJetID;
   		bool  _passJetCleaning;

        bool _isTight;
        bool _isLoose;

        int _ntrk;

        float _CSVv2;


        float _jet_partonFlavour  ;
        float _jet_hadronFlavour  ;

        float _jet_genJet_pt      ;
	    float _jet_genJet_E       ;
        /*
        float _jet_genJet_eta     ;
        float _jet_genJet_phi     ;
        float _jet_genJet_m       ;
        float _jet_genJet_E       ;
        float _jet_genJet_status  ;
        float _jet_genJet_id      ;
        */

        float _jet_genParton_pt     ;
	    float _jet_genParton_E      ;
        /*
        float _jet_genParton_eta    ;
        float _jet_genParton_phi    ;
        float _jet_genParton_m      ;
        float _jet_genParton_E      ;
        float _jet_genParton_status ;
        */
        float _jet_genParton_id     ;

        float _JES_uncert; // to be used with pt(1+/-uncert), E(1+/-uncert)

        float _pt_JER;
        float _pt_JER_down;
        float _pt_JER_up;

        ClassDef(Jet,1)
};

#endif
