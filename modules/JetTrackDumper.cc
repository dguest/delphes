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


/** \class JetTrackDumper
 *
 *  b-tagging algorithm based on counting tracks with large impact parameter
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetTrackDumper.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

//------------------------------------------------------------------------------

JetTrackDumper::JetTrackDumper() :
  fItTrackInputArray(0), fItJetInputArray(0)
{
}

//------------------------------------------------------------------------------

JetTrackDumper::~JetTrackDumper()
{
}

//------------------------------------------------------------------------------

void JetTrackDumper::Init()
{

  fPtMin = GetDouble("TrackPtMin", 1.0);
  fDeltaR = GetDouble("DeltaR", 0.3);
  fIPmax = GetDouble("TrackIPMax", 2.0);

  // import input array(s)

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/eflowTracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void JetTrackDumper::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void JetTrackDumper::Process()
{

  // loop over all input jets
  fItJetInputArray->Reset();
  Candidate* jet;
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
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

      // add tracks as jet candidates
      // TODO: make sure this doesn't mess with the downstream variables
      jet->AddCandidate(track);
    }
  }
}

//------------------------------------------------------------------------------
