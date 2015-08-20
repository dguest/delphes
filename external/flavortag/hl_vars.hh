#ifndef HL_VARS_HH
#define HL_VARS_HH

// Compute some high-level tagging variables.
//
// These structures are initially filled with garbage values: -1 for
// counters and NaN for floating point types. When the `fill()` method
// is called, these are replaced with either zero or +-inf when the
// tagger doesn't apply. In other words, NaN signals that something
// went wrong.

#include <vector>
#include <ostream>

#include "TVector3.h"

class SecondaryVertex;
// class TVector3;

// __________________________________________________________________________
// SVX

class SecondaryVertex: public TVector3
{
public:
  SecondaryVertex();
  SecondaryVertex(double, double, double);
  float Lxy;
  float Lsig;
  float decayLengthVariance;
  int nTracks;
  float eFrac;
  float mass;
  std::string config;
  // void Copy(SecondaryVertex& object) const;
  void clear();
  // TODO: get rid of ClassDef here, store as a flat vertex object in output
  ClassDef(SecondaryVertex, 1)
};

struct HighLevelSvx
{
  HighLevelSvx();
  void fill(const TVector3& jet, const std::vector<SecondaryVertex>&,
	    size_t skip_vx = 1);
  double Lsig;
  int NVertex;
  int NTracks;
  double DrJet;
  double Mass;
};

std::ostream& operator<<(std::ostream& os, const HighLevelSvx&);

// function to copy some info into a jet object
template<typename T>
void copy(const HighLevelSvx& from, T& to) {
#define CP(NAME) to.sv ## NAME = from.NAME
  CP(Lsig);
  CP(NVertex);
  CP(NTracks);
  CP(DrJet);
  CP(Mass);
#undef CP
}

// _________________________________________________________________________
// IP based

struct TrackParameters
{
  TrackParameters(const float[5], const float[15]);
  double d0;
  double z0;
  double phi;
  double d0err;
  double z0err;
};
std::ostream& operator<<(std::ostream& os, const TrackParameters&);

struct HighLevelTracking
{
  HighLevelTracking();
  void fill(const TVector3& jet, const std::vector<TrackParameters>&,
	    double ip_threshold = 1.8);
  double track2d0sig;
  double track3d0sig;
  double track2z0sig;
  double track3z0sig;
  int tracksOverIpThreshold;
  double jetProb;
};

std::ostream& operator<<(std::ostream& os, const HighLevelTracking&);

// function to copy some info into a jet object
template<typename T>
void copy(const HighLevelTracking& from, T& to) {
#define CP(NAME) to.NAME = from.NAME
  CP(track2d0sig);
  CP(track3d0sig);
  CP(track2z0sig);
  CP(track3z0sig);
  CP(tracksOverIpThreshold);
  CP(jetProb);
#undef CP
}


#endif
