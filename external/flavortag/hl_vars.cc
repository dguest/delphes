#include "hl_vars.hh"

#include "classes/DelphesClasses.h"

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
