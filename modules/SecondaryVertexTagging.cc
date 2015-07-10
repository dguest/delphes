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


/** \class SecondaryVertexTagging
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/SecondaryVertexTagging.h"

#include "classes/DelphesClasses.h"

#include "TObjArray.h"
// #include "TLorentzVector.h"

#include "rave/Version.h"
#include "rave/Track.h"
#include "rave/Exception.h"
#include "RaveBase/Converters/interface/RaveStreamers.h"
#include "RaveBase/Converters/interface/PerigeeToRaveObjects.h"
#include "rave/Vector6D.h"
#include "rave/VertexFactory.h"
#include "rave/ConstantMagneticField.h"

#include <iostream>

// class to convert Delphes stuff to Rave stuff
class RaveConverter
{
public:
  RaveConverter(double Bz);
  std::vector<rave::Track> getRaveTracks(const std::vector<Candidate*>& in);
private:
  rave::Vector6D getState(const Candidate*);
  rave::Covariance6D getCov(const Candidate*);
  double rho(double pt_in_gev, int charge); // return in cm^-1
  double _bz;
  rave::PerigeeToRaveObjects _converter;
};


//------------------------------------------------------------------------------

SecondaryVertexTagging::SecondaryVertexTagging() :
  fItTrackInputArray(0), fItJetInputArray(0), fMagneticField(0),
  fVertexFactory(0)
{
}

//------------------------------------------------------------------------------

SecondaryVertexTagging::~SecondaryVertexTagging()
{
  delete fMagneticField;
  delete fVertexFactory;
  delete fRaveConverter;
}

//------------------------------------------------------------------------------

void SecondaryVertexTagging::Init()
{
  // make sure Rave is working

  // read parameters
  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);
  // magnetic field
  fBz = GetDouble("Bz", 2.0);

  // import input array(s)

  fTrackInputArray = ImportArray(
    GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fJetInputArray = ImportArray(
    GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "secondaryVertices"));

  // initalize Rave
  std::cout << "** INFO:     This is Rave Version " << rave::Version()
	    << std::endl;
  fMagneticField = new rave::ConstantMagneticField(0, 0, fBz);
  fVertexFactory = new rave::VertexFactory(*fMagneticField);
  fRaveConverter = new RaveConverter(fBz);
}

//------------------------------------------------------------------------------

void SecondaryVertexTagging::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
}

//------------------------------------------------------------------------------


std::vector<Candidate*> SecondaryVertexTagging::GetTracks(Candidate* jet) {
  // loop over all input jets
  std::vector<Candidate*> jet_tracks;

  const TLorentzVector &jetMomentum = jet->Momentum;
  double jpx = jetMomentum.Px();
  double jpy = jetMomentum.Py();

  // loop over all input tracks
  fItTrackInputArray->Reset();
  Candidate* track;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trkMomentum = track->Momentum;

    double dr = jetMomentum.DeltaR(trkMomentum);

    double tpt = trkMomentum.Pt();
    double dxy = std::abs(track->Dxy);
    double ddxy = track->SDxy;

    if(tpt < fPtMin) continue;
    if(dr > fDeltaR) continue;
    if(dxy > fIPmax) continue;

    jet_tracks.push_back(track);
  }
  return jet_tracks;
}

void SecondaryVertexTagging::Process()
{
  Candidate *candidate;

  // loop over all input candidates
  fItJetInputArray->Reset();
  Candidate* jet;
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    TLorentzVector jet_momentum = jet->Momentum;
    printf("pt: %f\n", jet_momentum.Pt());

    auto tracks = GetTracks(jet);
    printf("n_tracks: %i\n", tracks.size());
    fRaveConverter->getRaveTracks(tracks);
    // if (false)
    // {
    //   fOutputArray->Add(candidate);
    // }
  }
}

//------------------------------------------------------------------------------
// definition of RaveConverter

// constants copied from ParticlePropagator
namespace {
  const double c_light = 2.99792458E8; // in [m/sec]
}


RaveConverter::RaveConverter(double Bz): _bz(Bz)
{
}

rave::Vector6D RaveConverter::getState(const Candidate* cand) {
  using namespace TrackParam;
  double a_d0 = cand->trkPar[D0];
  double a_z0 = cand->trkPar[Z0];
  // this is the _momentum_ phi (check?)
  double a_phi = cand->trkPar[PHI];
  double a_theta = cand->trkPar[THETA];
  // double a_qoverp = cand->trkPar[QOVERP];
  double a_q = cand->Charge;
  double a_pt = cand->Momentum.Pt();

  // -- translate these to Rave coordinates
  // rave base units are cm and GeV, Delphes takes mm and GeV

  double r_rho = rho(a_pt, a_q);
  double r_theta = a_theta; // should check this sign too...
  double r_phip = a_phi; 	// may be off by 90 degrees...
  // d0 is x cross p, where x is at perigee, and p is the initial
  // particle momentum. This isn't strictly accurate but it _should_
  // be a small effect for high pT tracks.
  double r_epsilon = a_d0 * 0.1; // have to check sign
  double r_zp = a_z0 * 0.1;

  // build the parameters
  rave::PerigeeParameters5D pars(r_rho, r_theta, r_phip, r_epsilon, r_zp);
  rave::Point3D referencePoint(0,0,0); // what is this?
  return _converter.convert(pars, a_q, referencePoint);
}

double RaveConverter::rho(double pt_in_gev, int charge) {
  // rave wants rho = 1/r, where r is the radius curvature
  // this is copied from ParticlePropagator
  // compute radius in [m]
  double radius = pt_in_gev / (charge * _bz) * 1.0E9/c_light;
  // convert back to rho, in cm
  return 1 / (1e2 * radius);
}

std::vector<rave::Track> RaveConverter::getRaveTracks(
  const std::vector<Candidate*>& in) {
  std::vector<rave::Track> tracks;

  for (const auto& cand: in) {
    rave::Vector6D state = getState(cand);
    std::cout << "particle state: " << state << std::endl;
  }
  return tracks;
}
