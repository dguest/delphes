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
 *  MERCHANTABILITY || FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class JetFlavorAssociation
 *
 *  Find origin of jet && evaluate jet flavor
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/SecondaryVertexAssociator.h"

#include "classes/DelphesClasses.h"

#include "TLorentzVector.h"

#include <iostream>
#include <cassert>


//------------------------------------------------------------------------------

SecondaryVertexAssociator::SecondaryVertexAssociator() :
  fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

SecondaryVertexAssociator::~SecondaryVertexAssociator()
{
}

//------------------------------------------------------------------------------

void SecondaryVertexAssociator::Init()
{
  ExRootConfParam param;

  // import input array(s)
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void SecondaryVertexAssociator::Finish()
{
  delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void SecondaryVertexAssociator::Process(){

  Candidate *jet;
  // TObjArray *partonArray = 0;
  std::cout << "new event" << std::endl;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    // loop over tracks
    Candidate* track;
    TIter itTracks(jet->GetTracks());
    // std::cout << "new jet " << jet->Flavor << std::endl;
    while ((track = static_cast<Candidate*>(itTracks.Next()))) {
      auto vertices = getHeavyFlavorVertices(track);
      for (const auto& vx: vertices) {
	std::cout << vx << std::endl;
      }
    }
  }
}

//------------------------------------------------------------------------------

namespace {
  int get_heaviest_particle(int pid);
  bool is_metastable(int pid);
}

SecondaryVertexAssociator::HFVs
SecondaryVertexAssociator::getHeavyFlavorVertices(Candidate* track) {
  if (track->GetCandidates()->GetEntriesFast() == 0) {
    return walkTruthRecord(track, {4, 5});
  }
  Candidate* mother = static_cast<Candidate*>(track->GetCandidates()->At(0));
  return getHeavyFlavorVertices(mother);
}

SecondaryVertexAssociator::HFVs
SecondaryVertexAssociator::walkTruthRecord(Candidate* genPart,
					   const std::set<int>& targets) {
  if (targets.size() == 0) {
    return {};
  }
  HFVs found;
  // check mother particles
  for (int mid: {genPart->M1, genPart->M2}) {
    if (mid == -1) continue;
    Candidate* mother = static_cast<Candidate*>(fParticleInputArray->At(mid));
    int heaviest = get_heaviest_particle(mother->PID);
    // if this is a (stable) target particle, record vertex...
    if (targets.count(heaviest) && is_metastable(mother->PID)) {
      auto newtarg = targets;
      const auto& pos = genPart->Position;
      HeavyFlavorVertex vx;
      vx.x = pos.X();
      vx.y = pos.Y();
      vx.z = pos.Z();
      vx.pdgid = mother->PID;
      vx.idx = mid;
      found.push_back(vx);
      // then remove from targets and call this function on the mother
      newtarg.erase(heaviest);
      for (auto& new_vx: walkTruthRecord(mother, newtarg)) {
	found.push_back(new_vx);
      }
    } else {
      // if this isn't a target, just call function on the mother
      for (auto& new_vx: walkTruthRecord(mother, targets)) {
	found.push_back(new_vx);
      }
    }
  }
  return found;
}

namespace {
  int get_heaviest_particle(int pid) {
    int absid = std::abs(pid);
    if (absid < 10) return absid;
    int tens = (absid % 100) / 10;
    int hundreds = (absid % 1000) / 100;
    int thousands = (absid % 10000) / 1000;
    // things with no quarks (hundreds and thousands == 0) just get abs pid
    if (hundreds == 0 && thousands == 0) {
      return absid;
    }
    // the leading digit should be larger than or equal to the lower ones
    // special exception for K_short, number 130
    assert(tens <= hundreds || tens <= thousands || absid == 130);
    assert(thousands == 0 || hundreds <= thousands);
    return std::max({tens, hundreds, thousands});
  }
  bool is_metastable(int pid) {
    int absid = std::abs(pid);
    // bare quarks should hadronize
    if (absid < 10) return false;

    // the only metastable particles we're interested in are c, s, b, tau
    if (absid == 15) return true;
    int quark = get_heaviest_particle(absid);
    bool metastable_quark = (
      (quark == 3) || (quark == 4) || (quark == 5) );
    if (!metastable_quark) return false;

    // high-spin hadrons decay fast
    int spin = absid % 10;
    if (spin > 2) return false;
    return true;
  }
}

std::ostream& operator<<(std::ostream& out, const HeavyFlavorVertex& vx) {
  out << "#" << vx.idx;
  out << " (" << vx.x << " " << vx.y << " " << vx.z << "): " << vx.pdgid;
  return out;
}
