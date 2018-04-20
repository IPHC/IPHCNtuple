
{
  gInterpreter->AddIncludePath("/afs/cern.ch/work/c/chanon/MG5_aMC_v2_5_5/Delphes/");
  gSystem->Load("libDelphes");
  gROOT->ProcessLine(".L TestExRootAnalysis.C+");
  TestExRootAnalysis();
}
