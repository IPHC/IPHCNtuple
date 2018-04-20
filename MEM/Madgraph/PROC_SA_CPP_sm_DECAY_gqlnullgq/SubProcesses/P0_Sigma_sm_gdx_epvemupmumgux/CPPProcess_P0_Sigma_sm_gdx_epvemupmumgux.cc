//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.0, 2015-07-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux.h"
#include "HelAmps_sm_gqlnullgq.h"

using namespace MG5_sm_gqlnullgq; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g d~ > e+ ve mu+ mu- g u~ WEIGHTED=10
// Process: g d~ > mu+ vm e+ e- g u~ WEIGHTED=10
// Process: g s~ > e+ ve mu+ mu- g c~ WEIGHTED=10
// Process: g s~ > mu+ vm e+ e- g c~ WEIGHTED=10

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_gqlnullgq::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[2]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 2; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 256; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, -1, 1, -1},
      {-1, -1, -1, -1, -1, -1, 1, 1}, {-1, -1, -1, -1, -1, 1, -1, -1}, {-1, -1,
      -1, -1, -1, 1, -1, 1}, {-1, -1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1,
      -1, 1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      -1, 1}, {-1, -1, -1, -1, 1, -1, 1, -1}, {-1, -1, -1, -1, 1, -1, 1, 1},
      {-1, -1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, -1, 1, 1, -1, 1}, {-1, -1,
      -1, -1, 1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1, 1}, {-1, -1, -1, 1, -1,
      -1, -1, -1}, {-1, -1, -1, 1, -1, -1, -1, 1}, {-1, -1, -1, 1, -1, -1, 1,
      -1}, {-1, -1, -1, 1, -1, -1, 1, 1}, {-1, -1, -1, 1, -1, 1, -1, -1}, {-1,
      -1, -1, 1, -1, 1, -1, 1}, {-1, -1, -1, 1, -1, 1, 1, -1}, {-1, -1, -1, 1,
      -1, 1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1, -1}, {-1, -1, -1, 1, 1, -1, -1,
      1}, {-1, -1, -1, 1, 1, -1, 1, -1}, {-1, -1, -1, 1, 1, -1, 1, 1}, {-1, -1,
      -1, 1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1, 1, -1, 1}, {-1, -1, -1, 1, 1, 1,
      1, -1}, {-1, -1, -1, 1, 1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1, -1, -1},
      {-1, -1, 1, -1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, -1, 1, -1}, {-1, -1,
      1, -1, -1, -1, 1, 1}, {-1, -1, 1, -1, -1, 1, -1, -1}, {-1, -1, 1, -1, -1,
      1, -1, 1}, {-1, -1, 1, -1, -1, 1, 1, -1}, {-1, -1, 1, -1, -1, 1, 1, 1},
      {-1, -1, 1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, -1, 1, 1}, {-1, -1, 1, -1, 1, 1,
      -1, -1}, {-1, -1, 1, -1, 1, 1, -1, 1}, {-1, -1, 1, -1, 1, 1, 1, -1}, {-1,
      -1, 1, -1, 1, 1, 1, 1}, {-1, -1, 1, 1, -1, -1, -1, -1}, {-1, -1, 1, 1,
      -1, -1, -1, 1}, {-1, -1, 1, 1, -1, -1, 1, -1}, {-1, -1, 1, 1, -1, -1, 1,
      1}, {-1, -1, 1, 1, -1, 1, -1, -1}, {-1, -1, 1, 1, -1, 1, -1, 1}, {-1, -1,
      1, 1, -1, 1, 1, -1}, {-1, -1, 1, 1, -1, 1, 1, 1}, {-1, -1, 1, 1, 1, -1,
      -1, -1}, {-1, -1, 1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, 1, -1, 1, -1}, {-1,
      -1, 1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, -1, 1}, {-1, 1, -1,
      -1, -1, -1, 1, -1}, {-1, 1, -1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, -1, 1,
      -1, -1}, {-1, 1, -1, -1, -1, 1, -1, 1}, {-1, 1, -1, -1, -1, 1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1, -1}, {-1, 1, -1,
      -1, 1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1, 1, -1}, {-1, 1, -1, -1, 1, -1,
      1, 1}, {-1, 1, -1, -1, 1, 1, -1, -1}, {-1, 1, -1, -1, 1, 1, -1, 1}, {-1,
      1, -1, -1, 1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1, 1}, {-1, 1, -1, 1, -1,
      -1, -1, -1}, {-1, 1, -1, 1, -1, -1, -1, 1}, {-1, 1, -1, 1, -1, -1, 1,
      -1}, {-1, 1, -1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, 1, -1, -1}, {-1, 1,
      -1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, -1, 1,
      1, 1}, {-1, 1, -1, 1, 1, -1, -1, -1}, {-1, 1, -1, 1, 1, -1, -1, 1}, {-1,
      1, -1, 1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1,
      1, -1, -1}, {-1, 1, -1, 1, 1, 1, -1, 1}, {-1, 1, -1, 1, 1, 1, 1, -1},
      {-1, 1, -1, 1, 1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, -1, 1, -1}, {-1, 1, 1, -1, -1, -1,
      1, 1}, {-1, 1, 1, -1, -1, 1, -1, -1}, {-1, 1, 1, -1, -1, 1, -1, 1}, {-1,
      1, 1, -1, -1, 1, 1, -1}, {-1, 1, 1, -1, -1, 1, 1, 1}, {-1, 1, 1, -1, 1,
      -1, -1, -1}, {-1, 1, 1, -1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1, 1, -1},
      {-1, 1, 1, -1, 1, -1, 1, 1}, {-1, 1, 1, -1, 1, 1, -1, -1}, {-1, 1, 1, -1,
      1, 1, -1, 1}, {-1, 1, 1, -1, 1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1, 1},
      {-1, 1, 1, 1, -1, -1, -1, -1}, {-1, 1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1,
      1, -1, -1, 1, -1}, {-1, 1, 1, 1, -1, -1, 1, 1}, {-1, 1, 1, 1, -1, 1, -1,
      -1}, {-1, 1, 1, 1, -1, 1, -1, 1}, {-1, 1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1,
      1, -1, 1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1, -1}, {-1, 1, 1, 1, 1, -1, -1,
      1}, {-1, 1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1, 1, -1, 1, 1}, {-1, 1, 1,
      1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, 1, 1, -1},
      {-1, 1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1, -1}, {1, -1, -1,
      -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      -1, 1, 1}, {1, -1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, -1, 1, -1, 1},
      {1, -1, -1, -1, -1, 1, 1, -1}, {1, -1, -1, -1, -1, 1, 1, 1}, {1, -1, -1,
      -1, 1, -1, -1, -1}, {1, -1, -1, -1, 1, -1, -1, 1}, {1, -1, -1, -1, 1, -1,
      1, -1}, {1, -1, -1, -1, 1, -1, 1, 1}, {1, -1, -1, -1, 1, 1, -1, -1}, {1,
      -1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, -1, 1, 1, 1, -1}, {1, -1, -1, -1,
      1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1, -1}, {1, -1, -1, 1, -1, -1, -1,
      1}, {1, -1, -1, 1, -1, -1, 1, -1}, {1, -1, -1, 1, -1, -1, 1, 1}, {1, -1,
      -1, 1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1, -1, 1}, {1, -1, -1, 1, -1,
      1, 1, -1}, {1, -1, -1, 1, -1, 1, 1, 1}, {1, -1, -1, 1, 1, -1, -1, -1},
      {1, -1, -1, 1, 1, -1, -1, 1}, {1, -1, -1, 1, 1, -1, 1, -1}, {1, -1, -1,
      1, 1, -1, 1, 1}, {1, -1, -1, 1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, 1, -1,
      1}, {1, -1, -1, 1, 1, 1, 1, -1}, {1, -1, -1, 1, 1, 1, 1, 1}, {1, -1, 1,
      -1, -1, -1, -1, -1}, {1, -1, 1, -1, -1, -1, -1, 1}, {1, -1, 1, -1, -1,
      -1, 1, -1}, {1, -1, 1, -1, -1, -1, 1, 1}, {1, -1, 1, -1, -1, 1, -1, -1},
      {1, -1, 1, -1, -1, 1, -1, 1}, {1, -1, 1, -1, -1, 1, 1, -1}, {1, -1, 1,
      -1, -1, 1, 1, 1}, {1, -1, 1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, 1, -1,
      -1, 1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, -1, 1, 1}, {1,
      -1, 1, -1, 1, 1, -1, -1}, {1, -1, 1, -1, 1, 1, -1, 1}, {1, -1, 1, -1, 1,
      1, 1, -1}, {1, -1, 1, -1, 1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1, -1}, {1,
      -1, 1, 1, -1, -1, -1, 1}, {1, -1, 1, 1, -1, -1, 1, -1}, {1, -1, 1, 1, -1,
      -1, 1, 1}, {1, -1, 1, 1, -1, 1, -1, -1}, {1, -1, 1, 1, -1, 1, -1, 1}, {1,
      -1, 1, 1, -1, 1, 1, -1}, {1, -1, 1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, 1, -1,
      -1, -1}, {1, -1, 1, 1, 1, -1, -1, 1}, {1, -1, 1, 1, 1, -1, 1, -1}, {1,
      -1, 1, 1, 1, -1, 1, 1}, {1, -1, 1, 1, 1, 1, -1, -1}, {1, -1, 1, 1, 1, 1,
      -1, 1}, {1, -1, 1, 1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1, 1, 1}, {1, 1, -1,
      -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, -1, 1}, {1, 1, -1, -1, -1,
      -1, 1, -1}, {1, 1, -1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, -1, 1, -1, -1},
      {1, 1, -1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, -1, 1, 1, 1}, {1, 1, -1, -1, 1, -1, -1, -1}, {1, 1, -1, -1, 1, -1,
      -1, 1}, {1, 1, -1, -1, 1, -1, 1, -1}, {1, 1, -1, -1, 1, -1, 1, 1}, {1, 1,
      -1, -1, 1, 1, -1, -1}, {1, 1, -1, -1, 1, 1, -1, 1}, {1, 1, -1, -1, 1, 1,
      1, -1}, {1, 1, -1, -1, 1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1, -1}, {1, 1,
      -1, 1, -1, -1, -1, 1}, {1, 1, -1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1, -1, 1}, {1,
      1, -1, 1, -1, 1, 1, -1}, {1, 1, -1, 1, -1, 1, 1, 1}, {1, 1, -1, 1, 1, -1,
      -1, -1}, {1, 1, -1, 1, 1, -1, -1, 1}, {1, 1, -1, 1, 1, -1, 1, -1}, {1, 1,
      -1, 1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, 1, -1, -1}, {1, 1, -1, 1, 1, 1, -1,
      1}, {1, 1, -1, 1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1, 1}, {1, 1, 1, -1,
      -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, -1, 1, 1}, {1, 1, 1, -1, -1, 1, -1, -1}, {1, 1, 1,
      -1, -1, 1, -1, 1}, {1, 1, 1, -1, -1, 1, 1, -1}, {1, 1, 1, -1, -1, 1, 1,
      1}, {1, 1, 1, -1, 1, -1, -1, -1}, {1, 1, 1, -1, 1, -1, -1, 1}, {1, 1, 1,
      -1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, -1, 1, 1}, {1, 1, 1, -1, 1, 1, -1,
      -1}, {1, 1, 1, -1, 1, 1, -1, 1}, {1, 1, 1, -1, 1, 1, 1, -1}, {1, 1, 1,
      -1, 1, 1, 1, 1}, {1, 1, 1, 1, -1, -1, -1, -1}, {1, 1, 1, 1, -1, -1, -1,
      1}, {1, 1, 1, 1, -1, -1, 1, -1}, {1, 1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, 1,
      -1, 1, -1, -1}, {1, 1, 1, 1, -1, 1, -1, 1}, {1, 1, 1, 1, -1, 1, 1, -1},
      {1, 1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, 1, -1, -1, -1}, {1, 1, 1, 1, 1,
      -1, -1, 1}, {1, 1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1, 1, -1, 1, 1}, {1, 1,
      1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, 1, 1,
      -1}, {1, 1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {96}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_gdx_epvemupmumgux(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_gdx_epvemupmumgux(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == -1)
  {
    // Add matrix elements for processes with beams (21, -1)
    return matrix_element[0] * 2; 
  }
  else if(id1 == 21 && id2 == -3)
  {
    // Add matrix elements for processes with beams (21, -3)
    return matrix_element[0] * 2; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]); 
  oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  ixxxxx(p[perm[7]], mME[7], hel[7], -1, w[7]); 
  FFV1_2(w[7], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[8]); 
  FFV2_3(w[2], w[3], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[9]); 
  FFV1_2(w[8], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[10]); 
  FFV2_1(w[5], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[11]); 
  FFV2_3(w[10], w[1], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[12]); 
  FFV1P0_3(w[4], w[5], pars->GC_3, pars->ZERO, pars->ZERO, w[13]); 
  FFV2_2(w[8], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_1(w[1], w[13], pars->GC_1, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[8], w[13], pars->GC_2, pars->ZERO, pars->ZERO, w[16]); 
  FFV2_1(w[1], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[17]); 
  VVV1_2(w[13], w[9], pars->GC_4, pars->mdl_MW, pars->mdl_WW, w[18]); 
  FFV2_4_3(w[4], w[5], pars->GC_50, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[19]);
  FFV2_3_1(w[1], w[19], pars->GC_50, pars->GC_58, pars->ZERO, pars->ZERO,
      w[20]);
  FFV2_5_2(w[8], w[19], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[21]);
  VVV1_1(w[9], w[19], pars->GC_53, pars->mdl_MW, pars->mdl_WW, w[22]); 
  FFV1_1(w[1], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[23]); 
  FFV2_3(w[8], w[23], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[24]); 
  FFV1_2(w[2], w[13], pars->GC_3, pars->ZERO, pars->ZERO, w[25]); 
  FFV2_4_2(w[2], w[19], pars->GC_50, pars->GC_59, pars->ZERO, pars->ZERO,
      w[26]);
  FFV2_1(w[3], w[19], pars->GC_62, pars->ZERO, pars->ZERO, w[27]); 
  VVV1P0_1(w[0], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[28]); 
  FFV1_2(w[7], w[28], pars->GC_11, pars->ZERO, pars->ZERO, w[29]); 
  FFV2_3(w[29], w[1], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[30]); 
  FFV1_1(w[1], w[28], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  FFV2_3(w[7], w[31], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[32]); 
  FFV2_2(w[7], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_2(w[7], w[13], pars->GC_2, pars->ZERO, pars->ZERO, w[34]); 
  FFV2_5_2(w[7], w[19], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[35]);
  FFV1_1(w[1], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_2(w[7], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV2_3(w[37], w[36], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[38]); 
  FFV2_1(w[36], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[39]); 
  FFV1_1(w[36], w[13], pars->GC_1, pars->ZERO, pars->ZERO, w[40]); 
  FFV2_3_1(w[36], w[19], pars->GC_50, pars->GC_58, pars->ZERO, pars->ZERO,
      w[41]);
  FFV1_1(w[36], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV2_3(w[7], w[42], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[43]); 
  FFV1_2(w[37], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV2_3(w[44], w[1], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[45]); 
  FFV2_2(w[37], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[46]); 
  FFV1_2(w[37], w[13], pars->GC_2, pars->ZERO, pars->ZERO, w[47]); 
  FFV2_5_2(w[37], w[19], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[48]);
  FFV1_2(w[33], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[49]); 
  FFV1_1(w[15], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  FFV1_2(w[34], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV1_1(w[17], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[52]); 
  FFV1_1(w[20], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_2(w[35], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  FFV1_1(w[23], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[55]); 
  FFV1_1(w[23], w[13], pars->GC_1, pars->ZERO, pars->ZERO, w[56]); 
  FFV2_1(w[23], w[9], pars->GC_100, pars->ZERO, pars->ZERO, w[57]); 
  FFV2_3_1(w[23], w[19], pars->GC_50, pars->GC_58, pars->ZERO, pars->ZERO,
      w[58]);
  FFV2_3(w[7], w[55], pars->GC_100, pars->mdl_MW, pars->mdl_WW, w[59]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[4], w[11], w[12], pars->GC_100, amp[0]); 
  FFV1_0(w[14], w[15], w[6], pars->GC_11, amp[1]); 
  FFV1_0(w[16], w[17], w[6], pars->GC_11, amp[2]); 
  FFV2_0(w[10], w[1], w[18], pars->GC_100, amp[3]); 
  FFV1_0(w[10], w[17], w[13], pars->GC_2, amp[4]); 
  FFV2_0(w[10], w[15], w[9], pars->GC_100, amp[5]); 
  FFV1_0(w[14], w[20], w[6], pars->GC_11, amp[6]); 
  FFV1_0(w[21], w[17], w[6], pars->GC_11, amp[7]); 
  FFV2_0(w[10], w[1], w[22], pars->GC_100, amp[8]); 
  FFV2_5_0(w[10], w[17], w[19], pars->GC_51, pars->GC_58, amp[9]); 
  FFV2_0(w[10], w[20], w[9], pars->GC_100, amp[10]); 
  FFV1_0(w[14], w[23], w[13], pars->GC_1, amp[11]); 
  FFV2_0(w[16], w[23], w[9], pars->GC_100, amp[12]); 
  VVV1_0(w[13], w[24], w[9], pars->GC_4, amp[13]); 
  FFV2_3_0(w[14], w[23], w[19], pars->GC_50, pars->GC_58, amp[14]); 
  FFV2_0(w[21], w[23], w[9], pars->GC_100, amp[15]); 
  VVV1_0(w[24], w[9], w[19], pars->GC_53, amp[16]); 
  FFV2_0(w[4], w[11], w[24], pars->GC_100, amp[17]); 
  FFV2_0(w[25], w[3], w[12], pars->GC_100, amp[18]); 
  FFV2_0(w[26], w[3], w[12], pars->GC_100, amp[19]); 
  FFV2_0(w[2], w[27], w[12], pars->GC_100, amp[20]); 
  FFV2_0(w[25], w[3], w[24], pars->GC_100, amp[21]); 
  FFV2_0(w[26], w[3], w[24], pars->GC_100, amp[22]); 
  FFV2_0(w[2], w[27], w[24], pars->GC_100, amp[23]); 
  FFV2_0(w[4], w[11], w[30], pars->GC_100, amp[24]); 
  FFV2_0(w[4], w[11], w[32], pars->GC_100, amp[25]); 
  FFV2_0(w[29], w[1], w[18], pars->GC_100, amp[26]); 
  FFV1_0(w[29], w[17], w[13], pars->GC_2, amp[27]); 
  FFV2_0(w[29], w[15], w[9], pars->GC_100, amp[28]); 
  FFV1_0(w[33], w[31], w[13], pars->GC_1, amp[29]); 
  FFV2_0(w[34], w[31], w[9], pars->GC_100, amp[30]); 
  FFV2_0(w[7], w[31], w[18], pars->GC_100, amp[31]); 
  FFV1_0(w[33], w[15], w[28], pars->GC_11, amp[32]); 
  FFV1_0(w[34], w[17], w[28], pars->GC_11, amp[33]); 
  FFV2_0(w[29], w[1], w[22], pars->GC_100, amp[34]); 
  FFV2_5_0(w[29], w[17], w[19], pars->GC_51, pars->GC_58, amp[35]); 
  FFV2_0(w[29], w[20], w[9], pars->GC_100, amp[36]); 
  FFV2_3_0(w[33], w[31], w[19], pars->GC_50, pars->GC_58, amp[37]); 
  FFV2_0(w[35], w[31], w[9], pars->GC_100, amp[38]); 
  FFV2_0(w[7], w[31], w[22], pars->GC_100, amp[39]); 
  FFV1_0(w[33], w[20], w[28], pars->GC_11, amp[40]); 
  FFV1_0(w[35], w[17], w[28], pars->GC_11, amp[41]); 
  FFV2_0(w[25], w[3], w[30], pars->GC_100, amp[42]); 
  FFV2_0(w[25], w[3], w[32], pars->GC_100, amp[43]); 
  FFV2_0(w[26], w[3], w[30], pars->GC_100, amp[44]); 
  FFV2_0(w[2], w[27], w[30], pars->GC_100, amp[45]); 
  FFV2_0(w[26], w[3], w[32], pars->GC_100, amp[46]); 
  FFV2_0(w[2], w[27], w[32], pars->GC_100, amp[47]); 
  FFV2_0(w[4], w[11], w[38], pars->GC_100, amp[48]); 
  VVV1_0(w[13], w[38], w[9], pars->GC_4, amp[49]); 
  FFV1_0(w[37], w[39], w[13], pars->GC_2, amp[50]); 
  FFV2_0(w[37], w[40], w[9], pars->GC_100, amp[51]); 
  VVV1_0(w[38], w[9], w[19], pars->GC_53, amp[52]); 
  FFV2_5_0(w[37], w[39], w[19], pars->GC_51, pars->GC_58, amp[53]); 
  FFV2_0(w[37], w[41], w[9], pars->GC_100, amp[54]); 
  FFV2_0(w[25], w[3], w[38], pars->GC_100, amp[55]); 
  FFV2_0(w[26], w[3], w[38], pars->GC_100, amp[56]); 
  FFV2_0(w[2], w[27], w[38], pars->GC_100, amp[57]); 
  FFV2_0(w[4], w[11], w[43], pars->GC_100, amp[58]); 
  FFV1_0(w[34], w[39], w[6], pars->GC_11, amp[59]); 
  FFV1_0(w[33], w[40], w[6], pars->GC_11, amp[60]); 
  FFV1_0(w[33], w[42], w[13], pars->GC_1, amp[61]); 
  FFV2_0(w[34], w[42], w[9], pars->GC_100, amp[62]); 
  FFV2_0(w[7], w[42], w[18], pars->GC_100, amp[63]); 
  FFV1_0(w[35], w[39], w[6], pars->GC_11, amp[64]); 
  FFV1_0(w[33], w[41], w[6], pars->GC_11, amp[65]); 
  FFV2_3_0(w[33], w[42], w[19], pars->GC_50, pars->GC_58, amp[66]); 
  FFV2_0(w[35], w[42], w[9], pars->GC_100, amp[67]); 
  FFV2_0(w[7], w[42], w[22], pars->GC_100, amp[68]); 
  FFV2_0(w[25], w[3], w[43], pars->GC_100, amp[69]); 
  FFV2_0(w[26], w[3], w[43], pars->GC_100, amp[70]); 
  FFV2_0(w[2], w[27], w[43], pars->GC_100, amp[71]); 
  FFV2_0(w[4], w[11], w[45], pars->GC_100, amp[72]); 
  FFV2_0(w[44], w[1], w[18], pars->GC_100, amp[73]); 
  FFV1_0(w[44], w[17], w[13], pars->GC_2, amp[74]); 
  FFV2_0(w[44], w[15], w[9], pars->GC_100, amp[75]); 
  FFV1_0(w[46], w[15], w[0], pars->GC_11, amp[76]); 
  FFV1_0(w[47], w[17], w[0], pars->GC_11, amp[77]); 
  FFV2_0(w[44], w[1], w[22], pars->GC_100, amp[78]); 
  FFV2_5_0(w[44], w[17], w[19], pars->GC_51, pars->GC_58, amp[79]); 
  FFV2_0(w[44], w[20], w[9], pars->GC_100, amp[80]); 
  FFV1_0(w[46], w[20], w[0], pars->GC_11, amp[81]); 
  FFV1_0(w[48], w[17], w[0], pars->GC_11, amp[82]); 
  FFV2_0(w[25], w[3], w[45], pars->GC_100, amp[83]); 
  FFV2_0(w[26], w[3], w[45], pars->GC_100, amp[84]); 
  FFV2_0(w[2], w[27], w[45], pars->GC_100, amp[85]); 
  FFV1_0(w[49], w[15], w[6], pars->GC_11, amp[86]); 
  FFV1_0(w[33], w[50], w[6], pars->GC_11, amp[87]); 
  FFV1_0(w[51], w[17], w[6], pars->GC_11, amp[88]); 
  FFV1_0(w[34], w[52], w[6], pars->GC_11, amp[89]); 
  FFV1_0(w[49], w[20], w[6], pars->GC_11, amp[90]); 
  FFV1_0(w[33], w[53], w[6], pars->GC_11, amp[91]); 
  FFV1_0(w[54], w[17], w[6], pars->GC_11, amp[92]); 
  FFV1_0(w[35], w[52], w[6], pars->GC_11, amp[93]); 
  FFV1_0(w[33], w[55], w[13], pars->GC_1, amp[94]); 
  FFV2_0(w[34], w[55], w[9], pars->GC_100, amp[95]); 
  FFV2_0(w[7], w[55], w[18], pars->GC_100, amp[96]); 
  FFV1_0(w[33], w[56], w[0], pars->GC_11, amp[97]); 
  FFV1_0(w[34], w[57], w[0], pars->GC_11, amp[98]); 
  FFV2_3_0(w[33], w[55], w[19], pars->GC_50, pars->GC_58, amp[99]); 
  FFV2_0(w[35], w[55], w[9], pars->GC_100, amp[100]); 
  FFV2_0(w[7], w[55], w[22], pars->GC_100, amp[101]); 
  FFV1_0(w[33], w[58], w[0], pars->GC_11, amp[102]); 
  FFV1_0(w[35], w[57], w[0], pars->GC_11, amp[103]); 
  FFV2_0(w[4], w[11], w[59], pars->GC_100, amp[104]); 
  FFV2_0(w[25], w[3], w[59], pars->GC_100, amp[105]); 
  FFV2_0(w[26], w[3], w[59], pars->GC_100, amp[106]); 
  FFV2_0(w[2], w[27], w[59], pars->GC_100, amp[107]); 

}
double CPPProcess_P0_Sigma_sm_gdx_epvemupmumgux::matrix_gdx_epvemupmumgux() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 108; 
  const int ncolor = 2; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3}; 
  static const double cf[ncolor][ncolor] = {{16, -2}, {-2, 16}}; 

  // Calculate color flows
  jamp[0] = +std::complex<double> (0, 1) * amp[24] + std::complex<double> (0,
      1) * amp[25] + std::complex<double> (0, 1) * amp[26] +
      std::complex<double> (0, 1) * amp[27] + std::complex<double> (0, 1) *
      amp[28] + std::complex<double> (0, 1) * amp[29] + std::complex<double>
      (0, 1) * amp[30] + std::complex<double> (0, 1) * amp[31] +
      std::complex<double> (0, 1) * amp[32] + std::complex<double> (0, 1) *
      amp[33] + std::complex<double> (0, 1) * amp[34] + std::complex<double>
      (0, 1) * amp[35] + std::complex<double> (0, 1) * amp[36] +
      std::complex<double> (0, 1) * amp[37] + std::complex<double> (0, 1) *
      amp[38] + std::complex<double> (0, 1) * amp[39] + std::complex<double>
      (0, 1) * amp[40] + std::complex<double> (0, 1) * amp[41] +
      std::complex<double> (0, 1) * amp[42] + std::complex<double> (0, 1) *
      amp[43] + std::complex<double> (0, 1) * amp[44] + std::complex<double>
      (0, 1) * amp[45] + std::complex<double> (0, 1) * amp[46] +
      std::complex<double> (0, 1) * amp[47] - amp[48] - amp[49] - amp[50] -
      amp[51] - amp[52] - amp[53] - amp[54] - amp[55] - amp[56] - amp[57] -
      amp[58] - amp[59] - amp[60] - amp[61] - amp[62] - amp[63] - amp[64] -
      amp[65] - amp[66] - amp[67] - amp[68] - amp[69] - amp[70] - amp[71] -
      amp[72] - amp[73] - amp[74] - amp[75] - amp[76] - amp[77] - amp[78] -
      amp[79] - amp[80] - amp[81] - amp[82] - amp[83] - amp[84] - amp[85] -
      amp[87] - amp[89] - amp[91] - amp[93];
  jamp[1] = -amp[0] - amp[1] - amp[2] - amp[3] - amp[4] - amp[5] - amp[6] -
      amp[7] - amp[8] - amp[9] - amp[10] - amp[11] - amp[12] - amp[13] -
      amp[14] - amp[15] - amp[16] - amp[17] - amp[18] - amp[19] - amp[20] -
      amp[21] - amp[22] - amp[23] - std::complex<double> (0, 1) * amp[24] -
      std::complex<double> (0, 1) * amp[25] - std::complex<double> (0, 1) *
      amp[26] - std::complex<double> (0, 1) * amp[27] - std::complex<double>
      (0, 1) * amp[28] - std::complex<double> (0, 1) * amp[29] -
      std::complex<double> (0, 1) * amp[30] - std::complex<double> (0, 1) *
      amp[31] - std::complex<double> (0, 1) * amp[32] - std::complex<double>
      (0, 1) * amp[33] - std::complex<double> (0, 1) * amp[34] -
      std::complex<double> (0, 1) * amp[35] - std::complex<double> (0, 1) *
      amp[36] - std::complex<double> (0, 1) * amp[37] - std::complex<double>
      (0, 1) * amp[38] - std::complex<double> (0, 1) * amp[39] -
      std::complex<double> (0, 1) * amp[40] - std::complex<double> (0, 1) *
      amp[41] - std::complex<double> (0, 1) * amp[42] - std::complex<double>
      (0, 1) * amp[43] - std::complex<double> (0, 1) * amp[44] -
      std::complex<double> (0, 1) * amp[45] - std::complex<double> (0, 1) *
      amp[46] - std::complex<double> (0, 1) * amp[47] - amp[86] - amp[88] -
      amp[90] - amp[92] - amp[94] - amp[95] - amp[96] - amp[97] - amp[98] -
      amp[99] - amp[100] - amp[101] - amp[102] - amp[103] - amp[104] - amp[105]
      - amp[106] - amp[107];

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



