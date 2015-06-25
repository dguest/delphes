/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class ImpactParameterSmearing
 *
 *  Performs transverse impact parameter smearing.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */


#include "modules/ImpactParameterSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

// ROOFIT
#include "RooMultiVarGaussian.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

using namespace std;
using namespace RooFit;
using namespace TrackParam;
//------------------------------------------------------------------------------

ImpactParameterSmearing::ImpactParameterSmearing() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ImpactParameterSmearing::~ImpactParameterSmearing()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Init()
{
  // read resolution formula

  fFormula->Compile(GetString("ResolutionFormula", "0.0"));

  TString filename_IDPara = GetString("SmearParamFile", "Parametrisation/IDParametrisierung.root");

  file_para = new TFile(filename_IDPara.Data(),"READ");
  if (!file_para->IsOpen() || file_para->IsZombie()) {
    throw std::runtime_error("bad file: " + string(filename_IDPara));
  }


  ptbins.push_back(10000);
  ptbins.push_back(20000);
  ptbins.push_back(50000);   
  ptbins.push_back(100000);   
  ptbins.push_back(200000);   
  ptbins.push_back(250000);   
  ptbins.push_back(500000);   
  ptbins.push_back(750000);   

  etabins.push_back(0.0);
  etabins.push_back(0.4);
  etabins.push_back(0.8);
  etabins.push_back(1.05);
  etabins.push_back(1.5);
  etabins.push_back(1.7);
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.7);

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ImpactParameterSmearing::Process()
{
  Candidate *candidate, *particle, *mother;
  double xd, yd, zd, dxy, sx, sy, sz, ddxy;
  double pt, eta, px, py, phi, e;
  int charge;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    // take momentum, position before smearing (otherwise apply double smearing)
    particle = static_cast<Candidate*>(candidate->GetCandidates()->At(0));

    const TLorentzVector &candidateMomentum = particle->Momentum;

    charge = particle->Charge;

    eta = candidateMomentum.Eta();
    pt = candidateMomentum.Pt();
    phi = candidateMomentum.Phi();
    e = candidateMomentum.E();
    
    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd =  particle->Xd;
    yd =  particle->Yd;
    zd =  particle->Zd;


    // Compute qoverp and theta: Because matrix parametrisation is for (d0,z0,phi,theta,qoverp)
    double qoverp = 1./(pt*cosh(eta));
    double theta = 2.*TMath::ATan(TMath::Exp(-eta));

    // calculate impact parameter (after-smearing)
    double d0 = (xd*py - yd*px)/pt;
    double z0 = zd; 


// create 5 dim vector to store corrections
  RooArgList xVec;// this is filled with the output of the MVG; e.g. this are the corrections to our 5 track params
  RooRealVar* x;
  for (int i = 0; i < 5; i++) {
    string xname = "x";
    if(i == 0) xname = "d0_corr";
    else if(i == 1) xname = "z0_corr";
    else if(i == 2) xname = "phi_corr";
    else if(i == 3) xname = "theta_corr";
    else if(i == 4) xname = "qoverp_corr";
    else cout << "Dim is only 5 ... not " << i << endl;
    // WARNING: this line also leaks memory (although not as bad as below)
    x = new RooRealVar(xname.c_str(), xname.c_str(), 0., -5.,5.);

    xVec.addOwned(*x);
  } 

  // get pt and eta bins
  int ptbin = 0;
  for(unsigned int i=0;i< ptbins.size();i++){
     if(pt > ptbins.at(i)) ptbin=i;
  }
  int etabin = 0;
  for(unsigned int i=0;i< ptbins.size();i++){
     if(fabs(eta) > etabins.at(i)) etabin=i;
  }

  // place holders for corrections; filled in loop
  TMatrixDSym* cov = new TMatrixDSym(5);
  TVectorD* muVec = new TVectorD(5);
  TString name;
  name.Form("covmat_ptbin%.2i_etabin%.2i",ptbin,etabin);
  file_para->GetObject(name,cov);
  if(!cov){
    cout << "No covariance matrix available for pt bin : " << ptbin << " and eta bin : " << etabin << endl;
    return;
  } 
  name.Form("meanvec_ptbin%.2i_etabin%.2i",ptbin,etabin);
  file_para->GetObject(name,muVec);
  if(!muVec){
    cout << "No mean vector available for pt bin : " << ptbin << " and eta bin : " << etabin << endl;
    return;
  }
  // flip sign of mean qoverp, because mean was computed using absolute values of |qoverp|: |qoverp_reco| - |qoverp|
  // this is not needed for qoverp elements of covariance matrix, because the matrix was computed using no absolute values
  //(*muVec)[4] *= charge;

  // now make the multivariate Gaussian
  RooMultiVarGaussian mvg ("mvg", "mvg", xVec, *muVec, *cov);
  // WARNING: this line leaks a _lot_ of memory! Should try to fix!
  RooDataSet* data = mvg.generate(xVec,1);

  //momentum correction is in MeV. Convert to GeV

  float d0corr     = data->get(0)->getRealValue("d0_corr");
  float z0corr     = data->get(0)->getRealValue("z0_corr");
  float phicorr    = data->get(0)->getRealValue("phi_corr");
  float thetacorr  = data->get(0)->getRealValue("theta_corr");
  float qoverpcorr = data->get(0)->getRealValue("qoverp_corr"); 

  //clean memory (cov is needed further down)
  delete muVec;
  delete data;

//cout <<"SC: MatrixSmeear d: " << d0 <<" "<< d0corr <<" ;z0: "<< z0 <<" " << z0corr <<" phi "<< phi <<" "<< phicorr 
//     <<" theta " << theta <<" "<< thetacorr <<" qoverp "<< qoverp <<" "<<qoverpcorr << endl; 

  float d0_reco = d0 + d0corr;
  float z0_reco = z0 + z0corr;
  float phi_reco = phi + phicorr;
  float theta_reco = theta + thetacorr;
  float qoverp_reco = qoverp + qoverpcorr*1000; //convert from MeV to GeV

// cout << "qoverp " << qoverp_reco <<" mu "<< (*muVec)[4] <<" corr "<< qoverpcorr << endl;


    // calculate impact parameter (after-smearing)
    dxy = d0_reco; //(xd*py - yd*px)/pt;

    ddxy =  gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    // reference for truth particle that smeared from
    mother = candidate;

    candidate = static_cast<Candidate*>(candidate->Clone());

    candidate->Xd = xd;
    candidate->Yd = yd;
    candidate->Zd = zd;

    candidate->Dxy = dxy;
    candidate->SDxy = ddxy;

  float* trkPar = candidate->trkPar;
  trkPar[D0]=d0_reco;
  trkPar[Z0]=z0_reco;
  trkPar[PHI]=phi_reco;
  trkPar[THETA]=theta_reco;
  trkPar[QOVERP]=qoverp_reco;


  float* trkCov = candidate->trkCov;
  trkCov[D0D0]=(*cov)(0,0);
  trkCov[Z0D0]=(*cov)(1,0);    trkCov[Z0Z0]=(*cov)(1,1);
  trkCov[PHID0]=(*cov)(2,0);  trkCov[PHIZ0]=(*cov)(2,1); trkCov[PHID0]=(*cov)(2,2);
  trkCov[THETAD0]=(*cov)(3,0);  trkCov[THETAZ0]=(*cov)(3,1); trkCov[THETAPHI]=(*cov)(3,2); trkCov[THETATHETA]=(*cov)(3,3);
  trkCov[QOVERPD0]=(*cov)(4,0);trkCov[QOVERPZ0]=(*cov)(4,1);trkCov[QOVERPPHI]=(*cov)(4,2);trkCov[QOVERPTHETA]=(*cov)(4,3);trkCov[QOVERPQOVERP]=(*cov)(4,4);
  delete cov;
  cov = 0;
    //candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*c_light*1.0E3);


    candidate->Momentum.SetPtEtaPhiM(charge/(qoverp_reco*cosh(eta)), -TMath::Log(TMath::Tan(theta_reco/2)), phi_reco, candidateMomentum.M()); 
    candidate->Dxy = d0_reco;
    candidate->SDxy = TMath::Sqrt(fabs(trkCov[D0D0]));

/*
cout <<"SC: SmearFinish "<< mother->Momentum.Pt() <<" -> "<<  candidate->Momentum.Pt() <<" "
                     << mother->Momentum.Eta() <<" -> "<<  candidate->Momentum.Eta()<<" "
                     << mother->Momentum.Phi() <<" -> "<<  candidate->Momentum.Phi()<<" "
                     << mother->Momentum.M() <<" -> "<<  candidate->Momentum.M()<<" "
                     << mother->trkPar[D0] <<" -> "<< trkPar[D0] <<" "
                     << mother->trkPar[Z0] <<" -> "<< trkPar[Z0] <<" "
                     << mother->trkPar[PHI] <<" -> "<< trkPar[PHI] <<" "
                     << mother->trkPar[THETA] <<" -> "<< trkPar[THETA] <<" "
                     << mother->trkPar[QOVERP] <<" -> "<< trkPar[QOVERP] <<" "
                     << endl;
*/

    TObjArray* array = (TObjArray*) candidate->GetCandidates();
    array->Clear() ;
    array->Add(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
