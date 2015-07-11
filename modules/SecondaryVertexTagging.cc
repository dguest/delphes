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
 *  Builds secondary vertices from tracks. Uses Rave.
 *
 *  \author Dan Guest
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
  rave::PerigeeCovariance5D getPerigeeCov(const Candidate*);
  double getRho(double pt_in_gev, int charge); // return in cm^-1
  double getRhoAlt(double qoverp, double theta);
  double getDrhoDqoverp(double theta);
  double getDrhoDtheta(double qoverp, double theta);
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
  // to do list
  std::cout << "** TODO: - make a b-tagger that works in Delphes\n"
	    // << "         - check / fix the covariance terms involving rho\n"
	    << std::flush;
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
  // double jpx = jetMomentum.Px();
  // double jpy = jetMomentum.Py();

  // loop over all input tracks
  fItTrackInputArray->Reset();
  Candidate* track;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trkMomentum = track->Momentum;

    double dr = jetMomentum.DeltaR(trkMomentum);

    double tpt = trkMomentum.Pt();
    double dxy = std::abs(track->Dxy);
    // double ddxy = track->SDxy;

    if(tpt < fPtMin) continue;
    if(dr > fDeltaR) continue;
    if(dxy > fIPmax) continue;

    jet_tracks.push_back(track);
  }
  return jet_tracks;
}

void SecondaryVertexTagging::Process()
{
  // loop over all input candidates
  fItJetInputArray->Reset();
  Candidate* jet;
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    TLorentzVector jet_momentum = jet->Momentum;
    printf("pt: %f\n", jet_momentum.Pt());

    auto tracks = GetTracks(jet);
    printf("n_tracks: %lu\n", tracks.size());
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
  double a_qoverp = cand->trkPar[QOVERP];
  double a_theta = cand->trkPar[THETA];
  double a_q = cand->Charge;

  // -- translate these to Rave coordinates
  // rave base units are cm and GeV, Delphes takes mm and GeV
  double r_rho = getRhoAlt(a_qoverp, a_theta);
  double r_theta = a_theta;
  double r_phip = a_phi;

  // TODO: check the Delphes Dxy / D0 definition.
  // d0 is x cross p, where x is at perigee, and p is the initial
  // particle momentum. This isn't strictly accurate but it _should_
  // be a small effect for high pT tracks.

  // convert d0 and z0 to cm
  double r_epsilon = a_d0 * 0.1;
  double r_zp = a_z0 * 0.1;

  // build the parameters
  rave::PerigeeParameters5D pars(r_rho, r_theta, r_phip, r_epsilon, r_zp);
  rave::Point3D referencePoint(0,0,0); // what is this?
  return _converter.convert(pars, a_q, referencePoint);
}

rave::PerigeeCovariance5D RaveConverter::getPerigeeCov(const Candidate* cand) {
  using namespace TrackParam;
  // -- translate to Rave coordinates
  // rave base units are cm and GeV, Delphes takes mm and GeV
  // need to calculate some things for the jacobian
  const float* par = cand->trkPar;
  double drdq = getDrhoDqoverp(par[THETA]);
  double drdt = getDrhoDtheta(par[QOVERP], par[THETA]);

  const float* cov = cand->trkCov;
  float qq = cov[QOVERPQOVERP];
  float qt = cov[QOVERPTHETA];
  float tt = cov[THETATHETA];

  // now multiply out the J*cov*J^T
  float drr = drdq*drdq*qq + drdt*drdt*tt + 2*drdt*drdq*qt;
  float drt = drdt*tt + drdq*qt;

  float drp = cov[QOVERPPHI]*drdq + cov[THETAPHI]*drdt;
  float dtt = cov[THETATHETA];
  float dtp = cov[THETAPHI];
  float dpp = cov[PHIPHI];
  rave::PerigeeCovariance3D cov3d(drr, drt, drp, dtt, dtp, dpp);

  // now the remaining terms. lengths need to be converted to cm
  float dre = (cov[QOVERPD0]*drdq + cov[THETAD0]*drdt) * 0.1;
  float drz = (cov[QOVERPZ0]*drdq + cov[THETAZ0]*drdt) * 0.1;
  float dte = cov[THETAD0]  * 0.1;
  float dtz = cov[THETAZ0]  * 0.1;
  float dpe = cov[PHID0]    * 0.1;
  float dpz = cov[PHIZ0]    * 0.1;
  float dee = cov[D0D0]     * 0.01;
  float dez = cov[Z0D0]     * 0.01;
  float dzz = cov[Z0Z0]     * 0.01;
  rave::PerigeeCovariance5D cov5d(cov3d,
				  dre, drz,
				  dte, dtz,
				  dpe, dpz,
				  dee, dez, dzz);
  return cov5d;
}


double RaveConverter::getRho(double pt_in_gev, int charge) {
  // rave wants rho = 1/r, where r is the radius curvature
  // this is copied from ParticlePropagator
  // compute radius in [m]
  double radius = pt_in_gev / (charge * _bz) * 1.0E9/c_light;
  printf("radius: %f [m]\n", radius);
  // convert back to rho, in cm
  return -1 / (1e2 * radius);
}

double RaveConverter::getRhoAlt(double qoverp, double theta) {
  // double eta = -std::log(std::tan(theta/2));
  double cosh_eta = 1/std::sin(theta);
  double rho = - qoverp * _bz * cosh_eta * c_light * 1e-9 * 1e-2;
  return rho;
}
double RaveConverter::getDrhoDqoverp(double theta) {
  double cosh_eta = 1/std::sin(theta);
  return - _bz * cosh_eta * c_light * 1e-11;
}
double RaveConverter::getDrhoDtheta(double qoverp, double theta) {
  return c_light*1e-11 * _bz*qoverp / (std::tan(theta) * std::sin(theta));
}

namespace {
  const Candidate* get_part(const Candidate* const_cand) {
    Candidate* cand = const_cast<Candidate*>(const_cand);
    if (cand->GetCandidates()->GetEntriesFast() == 0) {
      return cand;
    }
    Candidate* mother = static_cast<Candidate*>(cand->GetCandidates()->At(0));
    return get_part(mother);
  }
  void print_candidate_info(const Candidate* cand) {
    const TLorentzVector& mom = cand->Momentum;
    const TLorentzVector& pos = cand->Position;
    const Candidate* mother = get_part(cand);
    const TLorentzVector& mpos = mother->Position;

    float d0 = cand->Dxy;
    assert(d0 == cand->trkPar[TrackParam::D0]);
    float phi = cand->trkPar[TrackParam::PHI];
    float phi0 = phi - std::copysign(3.14159/2, 1);
    float x = d0 * std::cos(phi0);
    float y = d0 * std::sin(phi0);
    std::cout << "dphi: " << (phi0 - std::atan2(cand->Yd, cand->Xd))/3.1415
	      << "pi " << std::endl;
    std::cout << "initial x, y, z, r: " << mpos.X() << " " << mpos.Y() << " "
	      << mpos.Z() << " " << mpos.Perp() << std::endl;
    std::cout << "det edge x, y, z, r: " << pos.X() << " " << pos.Y() << " "
	      << pos.Z() << " " << pos.Perp() << std::endl;
    std::cout << "d0: " << d0 << " mm, charge: " << cand->Charge << ", PID: "
	      << cand->PID << ", pt " << mom.Pt() << " GeV"<< std::endl;
    std::cout << "used  x, y, z " << x << " " << y << " "
	      << cand->Zd << std::endl;

    std::cout << "other x, y, z " << cand->Xd << " " << cand->Yd << " "
	      << cand->Zd << std::endl;
    std::cout << "momentum x, y, z " << mom.Px() << " " << mom.Py() << " "
	      << mom.Pz() << std::endl;
  }
}

std::vector<rave::Track> RaveConverter::getRaveTracks(
  const std::vector<Candidate*>& in) {
  std::vector<rave::Track> tracks;
  for (const auto& cand: in) {
    const Candidate* particle = get_part(cand);
    if (particle->Position.Perp() > 0.0) {
      print_candidate_info(cand);
      rave::Vector6D state = getState(cand);
      std::cout << "particle state: " << state << std::endl;
      rave::PerigeeCovariance5D cov5d = getPerigeeCov(cand);
      int charge = cand->Charge;
      rave::Covariance6D cov6d = _converter.convert(cov5d, state, charge);
      std::cout << "particle cov: " << cov6d << std::endl;
    }
  }
  return tracks;
}
