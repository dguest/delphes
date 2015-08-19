#include "hl_vars.hh"

// #include "classes/DelphesClasses.h"
#include "constants_jetprob.hh"
#include "enums_track.hh"

#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <algorithm>

// #include <iostream>

class TrackParameters;

namespace {
  const double pi = std::atan2(0, -1);
  template<typename T>
  bool by_descending_first(std::pair<double, T> v1, std::pair<double, T> v2) {
    return v1.first > v2.first;
  }
  double get_jet_prob(const std::vector<TrackParameters>&);

  // see hardcoded parameters in constants_jetprob.hh
  double get_track_prob(double d0sig);
}

// __________________________________________________________________________
// SVX

SecondaryVertex::SecondaryVertex()
{
  clear();
}
SecondaryVertex::SecondaryVertex(double x, double y, double z):
  TVector3(x, y, z)
{
  clear();
}
void SecondaryVertex::clear() {
  Lxy = -1;
  Lsig = -1;
  nTracks = -1;
  eFrac = -1;
  mass = -1;
  config = "null";
}

void HighLevelSvx::fill(const TVector3& jvec,
			const std::vector<SecondaryVertex>& vertices,
			size_t skip) {
  double sum_dr_tracks = 0;
  int sum_tracks = 0;
  int sum_vertices = 0;
  double sum_mass = 0;

  // copied these variables from jetfitter
  double sum_sig = 0;
  double sum_inverr = 0;
  const size_t n_vtx = vertices.size();
  for (size_t vxn = skip; vxn < n_vtx; vxn++) {
    const auto& vx = vertices.at(vxn);
    sum_vertices++;
    double delta_r = jvec.DeltaR(vx);
    sum_dr_tracks += vx.nTracks * delta_r;
    sum_tracks += vx.nTracks;

    sum_mass += vx.mass;

    sum_sig += vx.Mag() / vx.decayLengthVariance;
    sum_inverr += 1/vx.decayLengthVariance;
  }
  bool has_vx = (sum_vertices) > 0 && (sum_tracks > 0);

  // save summary info
  Lsig = has_vx ? sum_sig / sqrt(sum_inverr) : 0;
  NVertex = has_vx ? sum_vertices               : 0;
  NTracks = has_vx ? sum_tracks                 : 0;
  DrJet = has_vx ? sum_dr_tracks / sum_tracks : 3.0;
  Mass = has_vx ? sum_mass                   : 0;
}

std::ostream& operator<<(std::ostream& os, const HighLevelSvx& hl) {
#define DUMP(VAR) os << #VAR ": " << hl.VAR << ", "
  DUMP(Lsig);
  DUMP(NVertex);
  DUMP(NTracks);
  DUMP(DrJet);
  DUMP(Mass);
#undef DUMP
  return os;
}

// ________________________________________________________________________
// Tracking

TrackParameters::TrackParameters(const float trkPar[5],
				 const float trkCov[15]):
  d0(trkPar[trk::D0]),
  z0(trkPar[trk::Z0]),
  phi(trkPar[trk::PHI]),
  d0err(std::sqrt(trkCov[trk::D0D0])),
  z0err(std::sqrt(trkCov[trk::Z0Z0]))
{
}

void HighLevelTracking::fill(const TVector3& jet,
			     const std::vector<TrackParameters>& pars,
			     double ip_threshold) {
  // size_t n_tracks = pars.size();
  double jet_phi = jet.Phi();
  assert(std::abs(jet_phi) <= pi);
  std::vector<std::pair<double, TrackParameters> > tracks_by_ip;

  // zero some things
  track2d0sig = 0;
  track2z0sig = 0;
  track3d0sig = 0;
  track3z0sig = 0;
  tracksOverIpThreshold = 0;
  jetProb = 0;

  if (pars.size() == 0) return;
  jetProb = get_jet_prob(pars);

  // what follows uses numbered tracks
  if (pars.size() < 2) return;

  for (const auto& par: pars) {
    double diff = std::abs(jet_phi - par.phi);
    int sign = (diff > 3*pi/4 || diff < pi/2) ? 1 : -1;
    double ip = std::copysign(par.d0, sign);
    tracks_by_ip.emplace_back(ip, par);
    if (ip > ip_threshold) tracksOverIpThreshold++;
  }

  std::sort(tracks_by_ip.begin(), tracks_by_ip.end(),
	    by_descending_first<TrackParameters>);
  track2d0sig = tracks_by_ip.at(1).second.d0;
  track2z0sig = tracks_by_ip.at(1).second.z0;
  if (tracks_by_ip.size() < 3) return;
  track3d0sig = tracks_by_ip.at(2).second.d0;
  track3z0sig = tracks_by_ip.at(2).second.z0;
}

namespace {
  double gauss_prob(double sig, const double norm, const double width) {
    const double sqp = std::sqrt(pi);
    const double sq2 = std::sqrt(2);
    return sqp / 2 * norm * width * (1 - std::erf( sig / (sq2 * width) ));
  }
  double exp_prob(double sig, const double off, const double mult) {
    return 1 / mult * std::exp(-off - mult*std::abs(sig));
  }
  double get_track_prob(double sig) {
    using namespace jetprob;
    double prob = gauss_prob(sig, P0, P1) + gauss_prob(sig, P2, P3) +
      exp_prob(sig, P4, P5) + exp_prob(sig, P6, P7);
    // std::cout << "prob for track with sig: "
    // 	      << sig << ": " << prob << std::endl;
    return prob;
  }
  double get_jet_prob(const std::vector<TrackParameters>& pars) {
    double p0 = 1.0;
    for (const auto& par: pars) {
      p0 *= get_track_prob(par.d0 / par.d0err);
    }
    int n_trk = pars.size();
    double corrections = 0;
    for (int k = 0; k < n_trk; k++) {
      corrections += std::pow( -std::log(p0), k) / std::tgamma(k + 1);
    }
    return p0 * corrections;
  }
}
