
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "Math/Integrator.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/GSLSimAnMinimizer.h"
#include "TPluginManager.h"

#include "Functions.h"
#include "SubGradient.h"
 
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, char *argv[])
{

  gPluginMgr->AddHandler("ROOT::Math::Minimizer", "subgradient", "SubGradientDescent", "SubGradient_h", "SubGradientDescent(const char*)");

  Functions* function = new Functions();
  function->par[0] = 1;
  function->par[1] = 1;
  function->par[2] = 0;

  ROOT::Math::Functor* toMinimize;
  toMinimize = new ROOT::Math::Functor(function, &Functions::Eval, 2);

  //SubGradientDescent* subgd = new SubGradientDescent();
  ROOT::Math::Minimizer* subgd = ROOT::Math::Factory::CreateMinimizer("subgradient","");
  subgd->SetFunction(*toMinimize);  

  subgd->SetLimitedVariable(0, "x", 0.5, 0.01, -10, 10);
  subgd->SetLimitedVariable(1, "y", 0.2, 0.01, -10, 10);

  subgd->SetMaxFunctionCalls(100000);
  subgd->Minimize();

  const double fBest = subgd->MinValue();
  const double* xBest = subgd->X();

  cout << "Minimum: "<<fBest<<endl;
  for (unsigned int ivar=0; ivar<subgd->NDim(); ivar++){
    cout << "Var"<<ivar<<": "<<xBest[ivar]<<endl;
  }

  return 0;
/*
  //ROOT::Minuit2::Minuit2Minimizer* minimizer = new ROOT::Minuit2::Minuit2Minimizer( ROOT::Minuit2::kMigrad );

  ROOT::Math::GSLSimAnMinimizer* minimizer = new ROOT::Math::GSLSimAnMinimizer();
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(100000);

  //minimizer->SetPrintLevel(0);
  minimizer->SetFunction(*toMinimize);

  minimizer->SetLimitedVariable(0, "x", 5, 0.01, -10, 10);

  minimizer->Minimize();

  const double *xs = minimizer->X();
  cout << function->Eval(xs) << endl;;

  return 0;
*/
}

