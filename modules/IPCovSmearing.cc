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

/** \class IPCovSmearing
 *
 *  Performs track smearing
 *
 *  \author Shih-Chieh Hsu
 *  \author Dan Guest (bug fixes)
 *
 */


#include "modules/IPCovSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

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

const double pi = std::atan2(0, -1);

//------------------------------------------------------------------------------

IPCovSmearing::IPCovSmearing() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

IPCovSmearing::~IPCovSmearing()
{
}

//------------------------------------------------------------------------------

void IPCovSmearing::Init()
{
  // read resolution formula

  TString filename_IDPara = GetString("SmearParamFile", "Parametrisation/IDParametrisierung.root");

  file_para = new TFile(filename_IDPara.Data(),"READ");
  if (!file_para->IsOpen() || file_para->IsZombie()) {
    throw std::runtime_error("bad file: " + string(filename_IDPara));
  }


  ptbins.push_back(10);
  ptbins.push_back(20);
  ptbins.push_back(50);
  ptbins.push_back(100);
  ptbins.push_back(200);
  ptbins.push_back(250);
  ptbins.push_back(500);
  ptbins.push_back(750);

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

void IPCovSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

namespace {
  void do_low_pt_hack(TMatrixDSym& matrix);
}

void IPCovSmearing::Process()
{
  Candidate *candidate, *particle, *mother;
  double xd, yd, zd;
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
    xd =  candidate->Xd;
    yd =  candidate->Yd;
    zd =  candidate->Zd;

    // NOTE: the phi used here isn't _strictly_ correct, since it doesn't
    //       extrapolate all the way to the interaction point.
    //       We're using it for consistency with the rest of Delphes, though.
    double phid0 = phi - pi/2;
    // NOTE: this definition is more correct, but not consistent with
    //       previous use.
    // double phid0 = std::atan2(yd, xd);

    // Compute qoverp and theta: Because matrix parametrisation is for (d0,z0,phi,theta,qoverp)
    double qoverp = charge/(pt*cosh(eta));
    double theta = 2.*TMath::ATan(TMath::Exp(-eta));

    // calculate impact parameter (_before_ smearing)
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
    x = new RooRealVar(xname.c_str(), xname.c_str(), 0., -5.,5.);

    xVec.addOwned(*x);
  }

  // get pt and eta bins
  int ptbin = -1;
  for(unsigned int i=0;i< ptbins.size();i++){
     if(pt > ptbins.at(i)) ptbin=i;
  }
  bool low_pt_hack = false;
  if (ptbin == -1) {
    low_pt_hack = true;
    ptbin = 0;
  }
  int etabin = -1;
  for(unsigned int i=0;i< etabins.size();i++){
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
  if (low_pt_hack) do_low_pt_hack(*cov);
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
  float phid0_reco = phid0 + phicorr;
  float theta_reco = theta + thetacorr;
  float qoverp_reco = qoverp + qoverpcorr*1000; //convert from MeV to GeV

// cout << "qoverp " << qoverp_reco <<" mu "<< (*muVec)[4] <<" corr "<< qoverpcorr << endl;

    // reference for truth particle that smeared from
    mother = candidate;

    candidate = static_cast<Candidate*>(candidate->Clone());

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
  
    //end of SC pasted code
    double smeared_pt = charge/(qoverp_reco*cosh(eta));
    assert(smeared_pt >= 0);
    candidate->Momentum.SetPtEtaPhiM(smeared_pt, -TMath::Log(TMath::Tan(theta_reco/2)), phi_reco, candidateMomentum.M());
    candidate->Dxy = d0_reco;
    candidate->SDxy = TMath::Sqrt(fabs(trkCov[D0D0]));

    // smear the Xd and Yd consistent with d0 smearing
    candidate->Xd = d0_reco * std::cos(phid0_reco);
    candidate->Yd = d0_reco * std::sin(phid0_reco);
    candidate->Zd = z0_reco;
    // if (d0_reco > 10) {
    //   printf("x %f -> %f ", xd, candidate->Xd);
    //   printf("y %f -> %f ", yd, candidate->Yd);
    //   printf("change: %f\n", d0corr);
    // }

    // cout <<"SC: SmearFinish "
    // << mother->Momentum.Pt() <<" -> "<<  candidate->Momentum.Pt() <<" "
    //   << mother->Momentum.Eta() <<" -> "<<  candidate->Momentum.Eta()<<" "
     // << mother->Momentum.Phi() <<" -> "<<  candidate->Momentum.Phi()<<" "
     // << mother->Momentum.M() <<" -> "<<  candidate->Momentum.M()<<" "
	 // "D0 "<< mother->trkPar[D0] <<" -> "<< trkPar[D0] <<" "
	 // << "(" << candidate->Dxy << ") "
    // 	 << "phi_x: " << phid0_reco << std::endl;
    // cout
    //  << mother->trkPar[Z0] <<" -> "<< trkPar[Z0] <<" "
    //  << mother->trkPar[PHI] <<" -> "<< trkPar[PHI] <<" "
    //  << mother->trkPar[THETA] <<" -> "<< trkPar[THETA] <<" "
    //  << mother->trkPar[QOVERP] <<" -> "<< trkPar[QOVERP] <<" "
    //  << endl;

    TObjArray* array = (TObjArray*) candidate->GetCandidates();
    array->Clear() ;
    array->Add(mother);

    fOutputArray->Add(candidate);
  }
}

namespace {
  void do_low_pt_hack(TMatrixDSym& cov_matrix){
    // hack to give larger uncertainty to low pt bins
    const int rank = 5;
    const double unct_mul = 2.0; // uncertainty increase for low pt

    // build a diagonal matrix
    TMatrixDSym hack_matrix(rank);
    hack_matrix.Zero();
    for (int iii = 0; iii < rank; iii++) hack_matrix[iii][iii] = 1;

    // add non-1 entries
    for (auto comp: {D0, Z0} ) hack_matrix[comp][comp] = unct_mul;

    TMatrix out_matrix = hack_matrix * cov_matrix * hack_matrix;

    // ugh, root sucks... How is there no assignment operator?
    for (int iii = 0; iii < rank; iii++) {
      for (int jjj = iii; jjj < rank; jjj++) {
    	cov_matrix[iii][jjj] = out_matrix[iii][jjj];
      }
    }
  }
}

//------------------------------------------------------------------------------
