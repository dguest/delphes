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

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>


//------------------------------------------------------------------------------

SecondaryVertexAssociator::SecondaryVertexAssociator() :
  fItParticleInputArray(0), fItJetInputArray(0)
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
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void SecondaryVertexAssociator::Finish()
{
  delete fItJetInputArray;
  delete fItParticleInputArray;
}

//------------------------------------------------------------------------------

void SecondaryVertexAssociator::Process(){

  Candidate *jet;
  TObjArray *partonArray = 0;

  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate *>(fItJetInputArray->Next())))
  {
    // DO THE STUFF
  }
}

//------------------------------------------------------------------------------
