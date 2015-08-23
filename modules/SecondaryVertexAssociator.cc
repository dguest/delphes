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
namespace {
  // generic pid-based information
  int get_heaviest_particle(int pid);
  int get_charge(int pid);
  bool is_metastable(int pid);

  // delphes specific
  int getNCharged(const std::vector<Candidate*>,
		  double pt_threshold = 1, double eta_threshold = 2.5);
  Candidate* getMother(Candidate* cand);

  // just for checks
  double getSmearingAngle(Candidate* track);
}

void SecondaryVertexAssociator::Process(){

  Candidate *jet;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    // loop over tracks
    Candidate* track;
    TIter itTracks(jet->GetTracks());
    std::map<TruthVertex, int> track_count;
    while ((track = static_cast<Candidate*>(itTracks.Next()))) {
      auto vertices = getHeavyFlavorVertices(track);
      for (const auto& vx: vertices) {
	track_count[vx]++;
      }
    }
    for (const auto& vx: track_count) {
      jet->truthVertices.push_back(vx.first);
    }
  }
}

//------------------------------------------------------------------------------


SecondaryVertexAssociator::HFVs
SecondaryVertexAssociator::getHeavyFlavorVertices(Candidate* track) {
  Candidate* mother = getMother(track);
  if (mother == 0) {
    return walkTruthRecord(track, {4, 5});
  }
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
      auto children = getStableChildren(getGenPart(mid));
      TruthVertex vx;
      vx.x = pos.X();
      vx.y = pos.Y();
      vx.z = pos.Z();
      vx.pdgid = mother->PID;
      vx.idx = mid;
      vx.n_charged_tracks = getNCharged(children);
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

std::vector<Candidate*> SecondaryVertexAssociator::getStableChildren(Candidate* mother)
{
  if (mother->Status == 1) {
    return {mother};
  }
  std::vector<Candidate*> children;
  for (int idx = mother->D1; idx <= mother->D2; idx++) {
    Candidate* child = getGenPart(idx);
    for (auto stable: getStableChildren(child)) {
      children.push_back(stable);
    }
  }
  return children;
}
Candidate* SecondaryVertexAssociator::getGenPart(int idx) {
  return static_cast<Candidate*>(fParticleInputArray->At(idx));
}


namespace {
  // _________________________________________________________________
  // generic pid-based information
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

  int lept_charge(int pid) {
    return std::copysign(std::abs(pid) % 2, -pid);
  }
  int quark_charge(int num) {
    return std::abs(num) % 2 == 1 ? -1: 2;
  }
  int had_charge(int pid) {
    int absid = std::abs(pid);
    int tens = (absid % 100) / 10;
    int hundreds = (absid % 1000) / 100;
    int thousands = (absid % 10000) / 1000;
    int frac_abs = 0;
    if (thousands == 0) {
      // mesons, one is an anti-particle
      frac_abs = std::abs(quark_charge(tens) - quark_charge(hundreds));
    } else {
      // otherwise, all have same charge
      frac_abs = quark_charge(tens) + quark_charge(hundreds) +
	quark_charge(thousands);
    }
    assert(frac_abs % 3 == 0);
    return std::copysign(frac_abs / 3, pid);
  }
  int get_charge(int pid) {
    // this is only supposed to work with (meta)stable particles,
    // mesons, baryons, and leptons
    int aid = std::abs(pid);
    if (aid > 10 && aid < 20) return lept_charge(pid);
    if (aid > 100 && aid < 1000000) return had_charge(pid);
    if (aid == 22) return 0;
    throw std::logic_error(__FILE__ ": no charge defined for " +
			   std::to_string(pid));
  }

  // _______________________________________________________________________
  // delphes specific
  int getNCharged(const std::vector<Candidate*> parts,
		  double pt_threshold,
		  double eta_threshold)
  {
    int n_charged = 0;
    for (const auto& cand: parts) {
      auto& mom = cand->Momentum;
      if (get_charge(cand->PID) && mom.Pt() > pt_threshold &&
	  std::abs(mom.Eta()) < eta_threshold) {
	n_charged++;
      }
    }
    return n_charged;
  }
  Candidate* getMother(Candidate* cand) {
    if (cand == 0) return 0;
    auto* subcand = cand->GetCandidates();
    if (subcand->GetEntriesFast() == 0) return 0;
    return static_cast<Candidate*>(subcand->At(0));
  }
  double getSmearingAngle(Candidate* track) {
    Candidate* unsmeared = getMother(track);
    Candidate* genpart = getMother(unsmeared);
    assert(genpart);
    Candidate* great_grandma = getMother(genpart);
    assert(!great_grandma);
    return unsmeared->Momentum.DeltaR(track->Momentum);
  }
}

