#ifndef TRUTH_H
#define TRUTH_H

#include "TObject.h"
#include "TLorentzVector.h"
#include "Base.h"

class Truth : public Base
{
    public:

        Truth();
        virtual ~Truth();

        int mc_truth_n()                        {return _mc_truth_n;};
        std::vector<int>   mc_truth_id()        {return _mc_truth_id;};
        std::vector<int>   mc_truth_label()     {return _mc_truth_label;};
        std::vector<float> mc_truth_pt()        {return _mc_truth_pt;};
        std::vector<float> mc_truth_eta()       {return _mc_truth_eta;};
        std::vector<float> mc_truth_phi()       {return _mc_truth_phi;};
        std::vector<float> mc_truth_E()         {return _mc_truth_E;};

        int gen_n()                             {return _gen_n;};
        std::vector<float> gen_pt()             {return _gen_pt;};
        std::vector<float> gen_eta()            {return _gen_eta;};
        std::vector<float> gen_phi()            {return _gen_phi;};
        std::vector<float> gen_m()              {return _gen_m;};
        std::vector<int>   gen_id()             {return _gen_id;};
        std::vector<int>   gen_status()         {return _gen_status;};
        std::vector<int>   gen_mother_id()      {return _gen_mother_id;};
        float gen_PVz()                         {return _gen_PVz;};

        float metGen_px()            	        {return _metGen_px;};
        float metGen_py()	     	            {return _metGen_py;};
        float metGen_pt()	     	            {return _metGen_pt;};
        float metGen_phi()	     	            {return _metGen_phi;};
        float metGen_sumet()	     	        {return _metGen_sumet;};
        float metGen_MuonEt()	     	        {return _metGen_MuonEt;};

        std::vector<int> Bjets_id()             {return _Bjets_id;};
        std::vector<int> Leptons_id()           {return _Leptons_id;};
        std::vector<int> Jets_id()              {return _Jets_id;};
        std::vector<int> AllJets_id()           {return _AllJets_id;};
        std::vector<int> JetsHighestPt_id()     {return _JetsHighestPt_id;};
        std::vector<int> JetsClosestMw_id()     {return _JetsClosestMw_id;};
        std::vector<int> JetsLowestMjj_id()     {return _JetsLowestMjj_id;};
        std::vector<int> QuarksFromWs_id()      {return _QuarksFromWs_id;};
        std::vector<int> JetsFromWs_id()        {return _JetsFromWs_id;};

        std::vector<float> Bjets_pt()           {return _Bjets_pt;};
        std::vector<float> Leptons_pt()         {return _Leptons_pt;};
        std::vector<float> Jets_pt()            {return _Jets_pt;};
        std::vector<float> AllJets_pt()         {return _AllJets_pt;};
        std::vector<float> JetsHighestPt_pt()   {return _JetsHighestPt_pt;};
        std::vector<float> JetsClosestMw_pt()   {return _JetsClosestMw_pt;};
        std::vector<float> JetsLowestMjj_pt()   {return _JetsLowestMjj_pt;};
        std::vector<float> QuarksFromWs_pt()    {return _QuarksFromWs_pt;};
        std::vector<float> JetsFromWs_pt()      {return _JetsFromWs_pt;};

        std::vector<float> Bjets_eta()          {return _Bjets_eta;};
        std::vector<float> Leptons_eta()        {return _Leptons_eta;};
        std::vector<float> Jets_eta()           {return _Jets_eta;};
        std::vector<float> AllJets_eta()        {return _AllJets_eta;};
        std::vector<float> JetsHighestPt_eta()  {return _JetsHighestPt_eta;};
        std::vector<float> JetsClosestMw_eta()  {return _JetsClosestMw_eta;};
        std::vector<float> JetsLowestMjj_eta()  {return _JetsLowestMjj_eta;};
        std::vector<float> QuarksFromWs_eta()   {return _QuarksFromWs_eta;};
        std::vector<float> JetsFromWs_eta()     {return _JetsFromWs_eta;};

        std::vector<float> Bjets_phi()          {return _Bjets_phi;};
        std::vector<float> Leptons_phi()        {return _Leptons_phi;};
        std::vector<float> Jets_phi()           {return _Jets_phi;};
        std::vector<float> AllJets_phi()        {return _AllJets_phi;};
        std::vector<float> JetsHighestPt_phi()  {return _JetsHighestPt_phi;};
        std::vector<float> JetsClosestMw_phi()  {return _JetsClosestMw_phi;};
        std::vector<float> JetsLowestMjj_phi()  {return _JetsLowestMjj_phi;};
        std::vector<float> QuarksFromWs_phi()   {return _QuarksFromWs_phi;};
        std::vector<float> JetsFromWs_phi()     {return _JetsFromWs_phi;};

        std::vector<float> Bjets_E()            {return _Bjets_E;};
        std::vector<float> Leptons_E()          {return _Leptons_E;};
        std::vector<float> Jets_E()             {return _Jets_E;};
        std::vector<float> AllJets_E()          {return _AllJets_E;};
        std::vector<float> JetsHighestPt_E()    {return _JetsHighestPt_E;};
        std::vector<float> JetsClosestMw_E()    {return _JetsClosestMw_E;};
        std::vector<float> JetsLowestMjj_E()    {return _JetsLowestMjj_E;};
        std::vector<float> QuarksFromWs_E()     {return _QuarksFromWs_E;};
        std::vector<float> JetsFromWs_E()       {return _JetsFromWs_E;};

        int boson_decay()  		                {return _boson_decay;};
        int ttbar_decay()  		                {return _ttbar_decay;};

        void read();
        void readMultiLepton();
        void init();


    protected:

        int _mc_truth_n;
        std::vector<int>    _mc_truth_id;
        std::vector<int>    _mc_truth_label;
        std::vector<float>  _mc_truth_pt;
        std::vector<float>  _mc_truth_eta;
        std::vector<float>  _mc_truth_phi;
        std::vector<float>  _mc_truth_E;

        int _gen_n;
        float _gen_PVz;
        std::vector<float>  _gen_pt;
        std::vector<float>  _gen_eta;
        std::vector<float>  _gen_phi;
        std::vector<float>  _gen_m;
        std::vector<int>    _gen_id;
        std::vector<int>    _gen_status;
        std::vector<int>    _gen_mother_id;

        float _metGen_px;
        float _metGen_py;
        float _metGen_pt;
        float _metGen_phi;
        float _metGen_sumet;
        float _metGen_MuonEt;

        std::vector<int> _Bjets_id;
        std::vector<int> _Leptons_id;
        std::vector<int> _Jets_id;
        std::vector<int> _AllJets_id;
        std::vector<int> _JetsHighestPt_id;
        std::vector<int> _JetsClosestMw_id;
        std::vector<int> _JetsLowestMjj_id;
        std::vector<int> _QuarksFromWs_id;
        std::vector<int> _JetsFromWs_id;

        std::vector<float> _Bjets_pt;
        std::vector<float> _Leptons_pt;
        std::vector<float> _Jets_pt;
        std::vector<float> _AllJets_pt;
        std::vector<float> _JetsHighestPt_pt;
        std::vector<float> _JetsClosestMw_pt;
        std::vector<float> _JetsLowestMjj_pt;
        std::vector<float> _QuarksFromWs_pt;
        std::vector<float> _JetsFromWs_pt;

        std::vector<float> _Bjets_eta;
        std::vector<float> _Leptons_eta;
        std::vector<float> _Jets_eta;
        std::vector<float> _AllJets_eta;
        std::vector<float> _JetsHighestPt_eta;
        std::vector<float> _JetsClosestMw_eta;
        std::vector<float> _JetsLowestMjj_eta;
        std::vector<float> _QuarksFromWs_eta;
        std::vector<float> _JetsFromWs_eta;

        std::vector<float> _Bjets_phi;
        std::vector<float> _Leptons_phi;
        std::vector<float> _Jets_phi;
        std::vector<float> _AllJets_phi;
        std::vector<float> _JetsHighestPt_phi;
        std::vector<float> _JetsClosestMw_phi;
        std::vector<float> _JetsLowestMjj_phi;
        std::vector<float> _QuarksFromWs_phi;
        std::vector<float> _JetsFromWs_phi;

        std::vector<float> _Bjets_E;
        std::vector<float> _Leptons_E;
        std::vector<float> _Jets_E;
        std::vector<float> _AllJets_E;
        std::vector<float> _JetsHighestPt_E;
        std::vector<float> _JetsClosestMw_E;
        std::vector<float> _JetsLowestMjj_E;
        std::vector<float> _QuarksFromWs_E;
        std::vector<float> _JetsFromWs_E;

        int _boson_decay;
        int _ttbar_decay;

        ClassDef(Truth,1)
};

#endif
