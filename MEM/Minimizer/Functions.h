

class Functions
{
  public:
  Functions();
  ~Functions();

  double Eval(const double* ) const;

  double GetMultiDimGradient(const double* ) const;
 
  double par[3];
 
  private:
};

Functions::Functions(){
}

Functions::~Functions(){

}

double Functions::Eval(const double* x) const {

  double xx = x[0];
  double yy = x[1];
  double f = par[0]*xx*xx + par[1]*yy*yy;// + par[1]*xx + par[2];
  
   if (f<-1) return 100;
   if (f>1) return 100;
   else return f;

}


