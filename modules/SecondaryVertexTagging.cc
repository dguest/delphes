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
#include "rave/Version.h"

#include "TObjArray.h"
// #include "TLorentzVector.h"

#include <iostream>


//------------------------------------------------------------------------------

SecondaryVertexTagging::SecondaryVertexTagging() :
  fItTrackInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

SecondaryVertexTagging::~SecondaryVertexTagging()
{
}

//------------------------------------------------------------------------------

void SecondaryVertexTagging::Init()
{
  // make sure Rave is working
  std::cout << "This is Rave Version " << rave::Version() << std::endl;

  // read parameters
  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);

  // import input array(s)

  fTrackInputArray = ImportArray(
    GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fJetInputArray = ImportArray(
    GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "secondaryVertices"));

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
    // if (false)
    // {
    //   fOutputArray->Add(candidate);
    // }
  }
}

//------------------------------------------------------------------------------
