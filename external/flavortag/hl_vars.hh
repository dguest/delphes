#ifndef HL_VARS_HH
#define HL_VARS_HH

#include <vector>
#include <ostream>

class SecondaryVertex;
class TVector3;

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

#endif
