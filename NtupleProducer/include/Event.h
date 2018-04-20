#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Event : public Base
{
    public:
        Event();
        virtual ~Event();

	double id()                       {return _id;};
        int run()                           {return _run;};
        int lumi()                          {return _lumi;};

        float rho()                         {return _rho;};

        float metpt()                       {return _metpt;};
        float metphi()                      {return _metphi;};
        float metpx()                       {return _metpt*cos(_metphi);};
        float metpy()                       {return _metpt*sin(_metphi);};
        float metsumet()                    {return _metsumet;};

        double metcov00()		            {return _metcov00;};
        double metcov01()		            {return _metcov01;};
        double metcov10()		            {return _metcov10;};
        double metcov11()		            {return _metcov11;};

        float metNoHF_pt()                  {return _metNoHF_pt;};
        float metNoHF_phi()                 {return _metNoHF_phi;};
        float metNoHF_sumet()               {return _metNoHF_sumet;};

        int pv_n()                          {return _pv_n;};
        float pv_z()                        {return _pv_z;};
        float pv_zError()                   {return _pv_zError;};

        float weight_scale_muF0p5()         {return _weight_scale_muF0p5;};
        float weight_scale_muF2()           {return _weight_scale_muF2;};
        float weight_scale_muR0p5()         {return _weight_scale_muR0p5;};
        float weight_scale_muR2()           {return _weight_scale_muR2;};
        std::vector<float> pdf_weights()    {return _pdf_weights;};
        std::vector<std::string> pdf_ids()  {return _pdf_ids;};

        float mc_weight()                   {return _mc_weight;};
        float mc_ptHat()                    {return _mc_ptHat;};
        int mc_pu_trueNumInt()              {return _mc_pu_trueNumInt;};

        bool  is_TRIGm()    {return  _TRIGm;   };
        bool  is_TRIGe()    {return  _TRIGe;   };
        bool  is_TRIGmTk()  {return  _TRIGmTk; };
        bool  is_TRIGee()   {return  _TRIGee;  };
        bool  is_TRIGmm()   {return  _TRIGmm;  };
        bool  is_TRIGme()   {return  _TRIGme;  };
        bool  is_TRIGem()   {return  _TRIGem;  };
        bool  is_TRIGmmTk() {return  _TRIGmmTk;};
        bool  is_TRIGeee()  {return  _TRIGeee; };
        bool  is_TRIGmme()  {return  _TRIGmme; };
        bool  is_TRIGeem()  {return  _TRIGeem; };
        bool  is_TRIGmmm()  {return  _TRIGmmm; };
        bool  is_TRIGee_noDz()  {return  _TRIGee_noDz;   };
        bool  is_TRIGmm_noDz()  {return  _TRIGmm_noDz;   };
        bool  is_TRIGmmTk_noDz(){return  _TRIGmmTk_noDz; };

        int tth_channel()                   {return _tth_channel;};

        float disc_TT()                     {return _disc_TT;};

        void read(bool isdata);
        void init();

    protected:

        long int   _id;
        int   _run;
        int   _lumi;

        float _rho;

        float _metpt;
        float _metphi;
        float _metsumet;
        double _metcov00;
        double _metcov01;
        double _metcov10;
        double _metcov11;

        float _metNoHF_pt;
        float _metNoHF_phi;
        float _metNoHF_sumet;

        int   _pv_n;
        float _pv_z;
        float _pv_zError;

        float _weight_scale_muF0p5;
        float _weight_scale_muF2;
        float _weight_scale_muR0p5;
        float _weight_scale_muR2;

        std::vector<float> _pdf_weights;
        std::vector<std::string> _pdf_ids;

        float _mc_weight;
        float _mc_ptHat;
        int   _mc_pu_trueNumInt;

       /* int   _trigger_pass;
        int   _trigger_pass_byname;
        int   _trigger_pass_byname_1;
        int   _trigger_pass_byname_1_noDz;*/

	    bool _TRIGm;
	    bool _TRIGe;
	    bool _TRIGmTk;
        bool _TRIGee;
	    bool _TRIGmm;
	    bool _TRIGme;
	    bool _TRIGem;
	    bool _TRIGmmTk;
        bool _TRIGeee;
	    bool _TRIGmme;
	    bool _TRIGeem;
	    bool _TRIGmmm;
        bool _TRIGee_noDz;
        bool _TRIGmm_noDz;
        bool _TRIGmmTk_noDz;

        float _disc_TT;

    public:
        int _tth_channel;

        ClassDef(Event,1)
};

#endif
