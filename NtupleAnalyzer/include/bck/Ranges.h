#ifndef SYSRANGES
#define SYSRANGES

namespace RANGE
{   
   double set_j_pt[3] = {25,0.,1.0};
   double set_j_JP[3] = {25,0.,2.5};
   double set_j_BJP[3] = {25,0.,12.};
   double set_j_ntrk[3] = {60,0.,60.};
   double set_j_eta[3] = {30,-3.,3.};
   double set_j_phi[3] = {30,-3.2,3.2};
   double set_j_mass[3] = {30,0.,100.};
   double set_j_mueta[3] = {30,-5.,5.};
   double set_j_muphi[3] = {30,-3.2,3.2};
   double set_j_muptrel[3] = {30,0.,30.};
   double set_j_svntrk[3] = {25,0.,25.};
   double set_j_svmass[3] = {30,0.,100.};
   double set_j_npv[3] = {25,0.,25.};
   double set_j_npu[3] = {50,0.,50.};
   double set_j_njet[3] = {12,0.,12.};
   double set_j_nmuon[3] = {5,0.,5.};
   double set_j_nsv[3] = {5,0.,5.};
   double set_j_nsvj[3] = {5,0.,5.};
   double set_j_nseltrk[3] = {60,0.,60.};
/*
   double set_j_JP[19][3] = {
      {25,0.,2.5}, // nosel	 
      {30,0.,1.4}, // pt20t30
      {30,0.,1.6}, // pt30t40
      {30,0.,1.8}, // pt40t50
      {30,0.,2.0}, // pt50t60
      {30,0.,2.1}, // pt60t70
      {30,0.,2.1}, // pt70t80
      {30,0.,2.2}, // pt80t100
      {30,0.,2.3}, // pt100t120
      {30,0.,2.4}, // pt120t160
      {30,0.,2.5}, // pt160t210
      {30,0.,2.5}, // pt210t260
      {30,0.,2.5}, // pt260t320
      {30,0.,2.5}, // pt320t400
      {30,0.,2.5}, // pt400t500
      {30,0.,2.5}, // pt500t600
      {30,0.,2.5}, // pt600t800
      {30,0.,2.5}, // pt800t1000
      {30,0.,2.5} // pt1000t1400
   };*/
   
   double set_j_mupt[19][3] = {
      {30,0.,110.}, // nosel	 
      {30,0.,25.}, // pt20t30
      {30,0.,25.}, // pt30t40
      {30,0.,30.}, // pt40t50
      {30,0.,30.}, // pt50t60
      {30,0.,35.}, // pt60t70
      {30,0.,40.}, // pt70t80
      {30,0.,45.}, // pt80t100
      {30,0.,50.}, // pt100t120
      {30,0.,55.}, // pt120t160
      {30,0.,70.}, // pt160t210
      {30,0.,90.}, // pt210t260
      {30,0.,100.}, // pt260t320
      {30,0.,110.}, // pt320t400
      {30,0.,110.}, // pt400t500
      {30,0.,110.}, // pt500t600
      {30,0.,110.}, // pt600t800
      {30,0.,110.}, // pt800t1000
      {30,0.,110.} // pt1000t1400
   };
   
   double set_j_2d_JP_vs_pt[6] = {
      25,0.,1.0,
      25,0.,2.5
   };
   
   double set_j_2d_JP_vs_mupt[6] = {
      50,0.,100.,
      25,0.,2.5
   };
   
   double set_j_2d_JP_vs_nseltrk[6] = {
      60,0.,60.,
      25,0.,2.5
   };

   double set_j_2d_JP_vs_njet[6] = {
      12,0.,12.,
      25,0.,2.5
   };
   
   double set_j_2d_JP_vs_jeta[6] = {
      25,-3.,3.,
      25,0.,2.5
   };

   double set_j_2d_JP_vs_nsv[6] = {
      5,0.,5.,
      25,0.,2.5
   };

   double set_j_2d_JP_vs_CSV[6] = {
      25,-0.1,1.1,
      25,0.,2.5
   };
   
   double set_j_2d_BJP_vs_CSV[6] = {
      25,-1.3,1.3,
      25,0.,12.
   };
   
   double set_j_2d_nseltrk_vs_njet[6] = {
      12,0.,12.,
      60,0.,60.
   };
   double set_j_2d_nseltrk_vs_ntrkgen[6] = {
      21,-1.,20.,
      60,0.,60.
   };

   double set_j_2d_nseltrk_vs_mupt[6] = {
      50,0.,100.,
      50,0.,60.
   };

   double set_j_2d_nseltrk_vs_jeta[6] = {
      25,-3.,3.,
      50,0.,60.
   };

   double set_j_2d_njet_vs_mupt[6] = {
      50,0.,100.,
      12,0.,12.
   };

   double set_j_2d_njet_vs_jeta[6] = {
      25,-3.,3.,
      12,0.,12.
   };

   double set_j_2d_mupt_vs_jeta[6] = {
      25,-3.,3.,
      50,0.,100.
   };

   double set_j_2d_nseltrk_vs_nsv[6] = {
      5,0.,5.,
      50,0.,60.
   };

   double set_j_2d_njet_vs_nsv[6] = {
      5,0.,5.,
      12,0.,12.
   };

   double set_j_2d_mupt_vs_nsv[6] = {
      5,0.,5.,
      50,0.,100.
   };

   double set_j_2d_jeta_vs_nsv[6] = {
      5,0.,5.,
      25,-3.,3.
   };

   double set_j_2d_pt_vs_nseltrk[6] = {
      50,0.,60.,
      25,0.,1.
   };
   
   double set_j_2d_pt_vs_njet[6] = {
      12,0.,12.,
      25,0.,1.
   };
   
   double set_j_2d_pt_vs_mupt[6] = {
      50,0.,100.,
      25,0.,1.
   };

   double set_j_2d_pt_vs_jeta[6] = {
      25,-3.,3.,
      25,0.,1.
   };

   double set_j_2d_pt_vs_nsv[6] = {
      5,0.,5.,
      25,0.,1.
   };
   
   // z y x
   double set_j_3d_ntrk_vs_pt_vs_mupt[9] = {
      30,0.,30.,
      25,0.,1.,
      25,0.,25.
   };
}

#endif
