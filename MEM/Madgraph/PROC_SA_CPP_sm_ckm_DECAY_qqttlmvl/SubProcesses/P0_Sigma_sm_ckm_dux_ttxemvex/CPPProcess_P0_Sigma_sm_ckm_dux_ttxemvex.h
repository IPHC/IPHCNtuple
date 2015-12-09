//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_ckm_qqttlmvl_dux_ttxemvex_H
#define MG5_Sigma_sm_ckm_qqttlmvl_dux_ttxemvex_H

#include <complex> 
#include <vector> 

#include "Parameters_sm_ckm_qqttlmvl.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: d u~ > t t~ e- ve~ WEIGHTED=6
// Process: d u~ > t t~ mu- vm~ WEIGHTED=6
// Process: s c~ > t t~ e- ve~ WEIGHTED=6
// Process: s c~ > t t~ mu- vm~ WEIGHTED=6
//--------------------------------------------------------------------------

class CPPProcess_P0_Sigma_sm_ckm_dux_ttxemvex
{
  public:

    // Constructor.
    CPPProcess_P0_Sigma_sm_ckm_dux_ttxemvex() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "d u~ > t t~ e- ve~ (sm_ckm)";}

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 6; 
    static const int nprocesses = 2; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 10; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 2; 
    std::complex<double> amp[namplitudes]; 
    double matrix_dux_ttxemvex(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm_ckm_qqttlmvl * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_ckm_qqttlmvl_dux_ttxemvex_H
