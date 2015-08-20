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


/** \class HDF5Writer
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/HDF5Writer.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootConfReader.h"
// as a hack we get the output file name from ExRootTreeWriter
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TObjArray.h"
#include "TFolder.h"

#include <iostream>
#include <string>

namespace {
  std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
  }
}

//------------------------------------------------------------------------------

HDF5Writer::HDF5Writer() :
  fItInputArray(0), m_out_file(0), m_jet_buffer(0)
{
}

//------------------------------------------------------------------------------

HDF5Writer::~HDF5Writer()
{
  delete m_out_file;
  delete m_jet_buffer;
  delete fItInputArray;
}

//------------------------------------------------------------------------------

void HDF5Writer::Init()
{
  fInputArray = ImportArray(
    GetString("JetInputArray", "UniqueObjectFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  auto* treeWriter = static_cast<ExRootTreeWriter*>(
    GetFolder()->FindObject("TreeWriter"));
  std::string output_file = remove_extension(treeWriter->GetOutputFileName());
  std::string output_ext = GetString("OutputExtension", ".ntuple.h5");
  output_file.append(output_ext);
  m_out_file = new H5::H5File(output_file, H5F_ACC_TRUNC);

  auto jtype = out::getJetType();

  m_jet_buffer = new OneDimBuffer<out::Jet>(*m_out_file, "jets", jtype, 1000);
}

namespace out {
  H5::CompType getJetType() {
    auto dtype = H5::PredType::NATIVE_DOUBLE;
    H5::CompType vertexType(sizeof(Vertex));
    vertexType.insertMember("mass", offsetof(Vertex, mass), dtype);
    vertexType.insertMember("dr_jet", offsetof(Vertex, dr_jet), dtype);
    auto verticesType = H5::VarLenType(&vertexType);
    H5::CompType jetType(sizeof(Jet));
    jetType.insertMember("pt", offsetof(Jet, pt), dtype);
    jetType.insertMember("eta", offsetof(Jet, eta), dtype);
    jetType.insertMember("vertices", offsetof(Jet, vertices.h5), verticesType);
    return jetType;
  }
}

//------------------------------------------------------------------------------

void HDF5Writer::Finish()
{
  m_jet_buffer->flush();
}

//------------------------------------------------------------------------------

void HDF5Writer::Process()
{
  fItInputArray->Reset();
  Candidate* jet;
  while ((jet = static_cast<Candidate*>(fItInputArray->Next()))) {
    const TLorentzVector& jvec = jet->Momentum;
    out::Jet outjet {jvec.Pt(), jvec.Eta()};
    const TVector3 j3vec = jvec.Vect();
    for (auto vx: jet->secondaryVertices) {
      outjet.vertices.push_back(out::Vertex{vx.mass, jvec.Vect().DeltaR(vx)});
    }
    m_jet_buffer->push_back(outjet);
  }
}

//------------------------------------------------------------------------------
