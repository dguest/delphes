#ifndef SECONDARY_VERTEX_HH
#define SECONDARY_VERTEX_HH

#include "TVector3.h"

class Candidate;

class SecondaryVertex: public TVector3
{
public:
  SecondaryVertex();
  SecondaryVertex(double, double, double);
  double Lxy;
  double Lsig;
  double decayLengthVariance;
  int nTracks;
  double eFrac;
  double mass;
  std::string config;
  std::vector<std::pair<double, Candidate*> > tracks;
  void clear();
};

#endif
