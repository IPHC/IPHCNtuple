#include "Math/IFunction.h"
//#include "Math/Functor.h"
#include "TRandom2.h"
#include "Math/Minimizer.h"
//#include "Fit/Fitter.h"

#include <iostream>
#include <ctime>

#define kFunctionValue 0

using namespace std;

class SubGradientDescent : public ROOT::Math::Minimizer //ROOT::Math::Minimizer
{
  public:


  int kStopCriterion;
  double tolerance;
  double alpha;
  bool foundNonNullGrad;
  int maxFunctionCalls;
  int verbose;

  SubGradientDescent();
  SubGradientDescent(const char*) :
    kStopCriterion(kFunctionValue),
    tolerance(1e-10),
    alpha(1.),
    foundNonNullGrad(false),
    maxFunctionCalls(10000),
    verbose(0)
  {}
  ~SubGradientDescent();


  //void SetFunction(ROOT::Math::Functor);
  virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func) ; 
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func)
   {
      SetFunction(static_cast<const ::ROOT::Math::IMultiGenFunction &> (func));
   }

  virtual bool SetVariable(unsigned int ivar, const std::string & name, double, double );
  virtual bool SetLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower, double upper); 
  void SetMaxFunctionCalls(unsigned int maxfcn);
  void SetPrintLevel(int);
  //unsigned int NDim();
  virtual unsigned int NDim() const ;
  virtual double MinValue() const ;

  void GetGradient();
  void NextPoint();
  bool CheckLimits();
  bool Minimize();

  TRandom2 rnd;
  bool CheckNullGrad();
  void RandomInit(); 
  void BackToBestPoint();

  bool CheckConvergence();
  //double* X();
  virtual const double *  X() const ; 

  void Printout();

  private:

  const ROOT::Math::IMultiGenFunction * 	fObjFunc;
  //ROOT::Math::Functor function;
  //const ROOT::Math::IMultiGenFunction function;
  //const ROOT::Math::IBaseFunctionMultiDim * fObjFunc;
  //const ROOT::Math::FitMethodFunction *function;
  double* x;
  double* xOld;
  double* xBest;
  double* xTmp;
  double* xprev;
  double* xnext;
  string* xName;
  double* xL;
  double* xU;
  double* steps;
  double* grad;
  //double* gradBest;

  //int verbose;
  int ntimesSamePoint;
  int nFunctionCalls;
  //int maxFunctionCalls;
  //bool foundNonNullGrad;
};

SubGradientDescent::SubGradientDescent(){
/*
  if (verbose>0) cout << "Instantiating SubGradientDescent"<<endl;
  kStopCriterion = kFunctionValue;
  tolerance = 1e-10;
  alpha = 1.;
  foundNonNullGrad = false;
  verbose = 1;
*/
}
/*
SubGradientDescent::SubGradientDescent(const char* a){
  //SubGradientDescent();
  if (verbose>0) cout << "Instantiating SubGradientDescent"<<endl;
  kStopCriterion = kFunctionValue;
  tolerance = 1e-10;
  alpha = 1.;
  foundNonNullGrad = false;
  verbose = 1;
}
*/
SubGradientDescent::~SubGradientDescent(){
}

void SubGradientDescent::SetFunction(const ROOT::Math::IMultiGenFunction & func){
//void SubGradientDescent::SetFunction(ROOT::Math::Functor functor){

  if (verbose>0) cout << "Setting function" <<endl;

  //function = ROOT::Math::Functor(func);
  //const ROOT::Math::FitMethodFunction * function_tmp = dynamic_cast<const ROOT::Math::FitMethodFunction *>(&func);
  //fObjFunc = new ROOT::Math::IMultiGenFunction(func);//dynamic_cast<const ROOT::Math::IMultiGenFunction &> (func);
  //ROOT::Math::Minimizer::SetFunction(static_cast<const ::ROOT::Math::IMultiGenFunction &> (func));
  //Minimizer::SetFunction(func);
  fObjFunc = static_cast<const ROOT::Math::IMultiGenFunction*> (&func);
  x = new double[func.NDim()];
  xOld = new double[func.NDim()];
  xBest = new double[func.NDim()];
  xTmp = new double[func.NDim()];
  xprev = new double[func.NDim()];
  xnext = new double[func.NDim()];
  xL = new double[func.NDim()];
  xU = new double[func.NDim()]; 
  steps = new double[func.NDim()];
  grad = new double[func.NDim()];
  xName = new string[func.NDim()];

  if (verbose>0) cout << "Function with "<<func.NDim()<<" variables"<<endl;

}

unsigned int SubGradientDescent::NDim() const{
  return fObjFunc->NDim();
}

bool SubGradientDescent::SetVariable(unsigned int ivar, const std::string & name, double val, double step){

  x[ivar] = val;
  xName[ivar] = name;
  steps[ivar] = step;

  return true;
}

bool SubGradientDescent::SetLimitedVariable(unsigned int ivar, const std::string & name, double val, double step, double lower, double uper){ 
//void SubGradientDescent::SetLimitedVariable(unsigned int ivar, string name, double val, double step, double valmin, double valmax){

  x[ivar] = val;
  xName[ivar] = name;
  steps[ivar] = step;
  xL[ivar] = lower;
  xU[ivar] = uper;

  return true;
}

void SubGradientDescent::SetMaxFunctionCalls(unsigned int maxfcn){

  maxFunctionCalls = maxfcn;
  if (verbose>0) cout << "SetMaxFunctionCalls maxFunctionCalls="<<maxFunctionCalls<<endl;

}

void SubGradientDescent::SetPrintLevel(int v){
  verbose = v;
}

void SubGradientDescent::GetGradient(){

  //double gradient=0;
  double subgrad_minus, subgrad_plus, grad_chosen=0;

  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++) {
    xprev[ivar] = x[ivar];
    xnext[ivar] = x[ivar];
  }

  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    xprev[ivar] = x[ivar] - steps[ivar];
    xnext[ivar] = x[ivar] + steps[ivar];
    //cout << "ivar="<<ivar<<" xprev="<<xprev[ivar]<<" f(xprev)="<<function(xprev)  <<" xnext="<<xnext[ivar] << " f(xnext)="<<function(xnext)<<endl;

    //gradient = ((*fObjFunc)(xnext)-(*fObjFunc)(xprev))/(2*steps[ivar]);
    subgrad_minus = ((*fObjFunc)(x)-(*fObjFunc)(xprev))/steps[ivar];
    subgrad_plus = ((*fObjFunc)(xnext)-(*fObjFunc)(x))/steps[ivar];
    //cout << "Subgradients of "<< xName[ivar] << " with step "<<steps[ivar]<<" : grad_minus="<<subgrad_minus<<" grad_plus="<<subgrad_plus<<" grad="<<gradient<<endl;
    
    if (fabs(subgrad_minus) < fabs(subgrad_plus)) grad_chosen = subgrad_minus;
    else if (fabs(subgrad_minus) >= fabs(subgrad_plus)) grad_chosen = subgrad_plus;

    //if (grad_chosen>1000) grad_chosen = 0;
 
    //if (subgrad_minus==0 && subgrad_plus!=0) grad_chosen = subgrad_plus;
    //else if (subgrad_minus!=0 && subgrad_plus==0) grad_chosen = subgrad_minus;

    xprev[ivar] = x[ivar];
    xnext[ivar] = x[ivar];

    grad[ivar] = grad_chosen;
  }

}

void SubGradientDescent::NextPoint(){

  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    xOld[ivar] = x[ivar];
    x[ivar] = x[ivar] - alpha*steps[ivar]*grad[ivar];
  }

  if ((*fObjFunc)(x) < (*fObjFunc)(xBest)){
    for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++) {
      xBest[ivar] = x[ivar];
      //gradBest[ivar] = grad[ivar];
    }
    ntimesSamePoint = 0;
  }

}

bool SubGradientDescent::CheckLimits(){

  //bool revertPos = false;

  bool moveToLimit = false;
  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    if (x[ivar] < xL[ivar]+steps[ivar]) {
      moveToLimit=true; 
      //x[ivar] = xL[ivar]+steps[ivar]; }
    }
    if (x[ivar] > xU[ivar]-steps[ivar]) {
      moveToLimit=true; 
      //x[ivar] = xU[ivar]-steps[ivar]; }
    }
  }
  if (moveToLimit) {
    if (verbose>0) cout << "Position above limit" << endl;
  }
/*
  if (revertPos){
    for (unsigned int ivar=0; ivar<function.NDim(); ivar++){
      xTmp[ivar] = x[ivar];
      x[ivar] = xOld[ivar];
      xOld[ivar] = xTmp[ivar];
    }
  }
*/
 return moveToLimit; 
}

bool SubGradientDescent::CheckNullGrad(){

  bool newPos = true;
  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    if (grad[ivar]!=0) {
      newPos = false;
    }
  }

  if (foundNonNullGrad==false && newPos==false) foundNonNullGrad=true;

  return newPos;

}

void SubGradientDescent::RandomInit(){

  if (verbose>0) cout << "Initialize randomly"<<endl;
  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    x[ivar] = xL[ivar] + rnd.Uniform()*(xU[ivar]-xL[ivar]);
  }
}

void SubGradientDescent::BackToBestPoint(){

 for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++) x[ivar] = xBest[ivar];

 ntimesSamePoint++;
 
 alpha /= 2;
 if (verbose>0) cout << "Back to best point, alpha="<<alpha<<endl;

}

bool SubGradientDescent::Minimize(){

  if (verbose>0) cout <<"kStopCriterion="<< kStopCriterion << " tolerance="<<tolerance<<" alpha="<<alpha<<" foundNonNullGrad="<<foundNonNullGrad<<" maxFunctionCalls="<<maxFunctionCalls<<" verbose="<<verbose<<endl;

  nFunctionCalls = 0;
  alpha = 1;

  bool c = false;
  bool isNull = false;
  bool isAboveLimit = false;
  ntimesSamePoint = 0;

  if ((*fObjFunc)(x) < 10000) for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++) xBest[ivar] = x[ivar]; 

  do {

    c = false;
    isNull = false;
    isAboveLimit = false;

    GetGradient();
    Printout();

    NextPoint();

    isAboveLimit = CheckLimits();
    isNull = CheckNullGrad();

    if (isAboveLimit && !foundNonNullGrad) RandomInit();
    if (isAboveLimit && foundNonNullGrad) BackToBestPoint();

    if (isNull && !foundNonNullGrad) RandomInit(); 
    else if (isNull && foundNonNullGrad) BackToBestPoint();

    if (!isNull) c = CheckConvergence();

    nFunctionCalls++;
    //if (nFunctionCalls==1 || nFunctionCalls % 100) Printout();

  } while (!c && nFunctionCalls<maxFunctionCalls);

  Printout();

  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++) x[ivar] = xBest[ivar];
  if (verbose>0) cout << "Minmization result:" <<endl;
  Printout();

  return c;
}

bool SubGradientDescent::CheckConvergence(){

  //bool isSamePoint = false;
  bool c = false;

  //if (function(x)==function(xBest)) ntimesSamePoint++;
  //else ntimesSamePoint = 0;

  if (kStopCriterion==kFunctionValue){
    if (fabs((*fObjFunc)(x)-(*fObjFunc)(xOld))<tolerance) c = true; 
  }

  if (ntimesSamePoint<5 && c) {if (verbose>0) cout << "Tolerance="<<tolerance<<", stopping"<<endl;}
  if (ntimesSamePoint==5) {if (verbose>0) cout <<"ntimesSamePoint="<<ntimesSamePoint<<", stopping "<<endl;}

   return ((ntimesSamePoint<5 && c) || ntimesSamePoint==5);
}

const double* SubGradientDescent::X() const{

  return xBest;

}

double SubGradientDescent::MinValue() const {

  return (*fObjFunc)(xBest);

}

void SubGradientDescent::Printout(){

  if (verbose>0) cout << "-> f(x)="<<(*fObjFunc)(x)<<" alpha="<<alpha<<endl;
  for (unsigned int ivar=0; ivar<fObjFunc->NDim(); ivar++){
    if (verbose>0) cout << "---> "<< xName[ivar]<<"="<< x[ivar] << " grad="<<grad[ivar]<< " step="<<steps[ivar]<<endl;
  }
}


