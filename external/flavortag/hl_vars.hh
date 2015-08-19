#ifndef HL_VARS_HH
#define HL_VARS_HH

#include <vector>
#include <ostream>

class SecondaryVertex;
class TVector3;
class Candidate;

// __________________________________________________________________________
// SVX

struct HighLevelSvx
{
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
  TrackParameters(const Candidate&);
  // void fillFromParArray(const float[5], const float[15])
  double d0;
  double z0;
  double phi;
  double d0err;
  double z0err;
};

struct HighLevelTracking
{
  void fill(const TVector3& jet, const std::vector<TrackParameters>&,
	    double ip_threshold = 1.8);
  double track2d0sig;
  double track3d0sig;
  double track2z0sig;
  double track3z0sig;
  double tracksOverIpThreshold;
  double jetProb;
};

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
