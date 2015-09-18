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
 * writes jet information to HDF5
 *
 *  \author Dan Guest
 *
 */

#include "modules/HDF5Writer.h"
#include "external/h5/h5types.hh"
#include "external/flavortag/SecondaryVertex.hh"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootConfReader.h"
// as a hack we get the output file name from ExRootTreeWriter
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TObjArray.h"
#include "TFolder.h"

#include <iostream>
#include <string>
#include <limits>

namespace {
  // double inf = std::numeric_limits<double>::infinity();
  std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
  }
}

//------------------------------------------------------------------------------

HDF5Writer::HDF5Writer() :
  fItInputArray(0), m_out_file(0), m_hl_jet_buffer(0),
  m_ml_jet_buffer(0)
{
}

//------------------------------------------------------------------------------

HDF5Writer::~HDF5Writer()
{
  delete m_out_file;
  delete m_hl_jet_buffer;
  delete m_ml_jet_buffer;
  delete fItInputArray;
}

//------------------------------------------------------------------------------
namespace out {
  // insering a compound type requires that `type(Class)` is defined
  H5::CompType type(Vertex);
  H5::CompType type(Jet);
}

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

  auto hl_jtype = out::type(out::HighLevelJet());
  auto ml_jtype = out::type(out::MediumLevelJet());

  m_hl_jet_buffer = new OneDimBuffer<out::HighLevelJet>(
    *m_out_file, "high_level_jets", hl_jtype, 1000);
  m_ml_jet_buffer = new OneDimBuffer<out::MediumLevelJet>(
    *m_out_file, "medium_level_jets", ml_jtype, 1000);
}

namespace out {
  int simple_flavor(int flav) {
    switch (flav) {
    case 4: return 4;
    case 5: return 5;
    case 15: return 15;
    default: return 0;
    }
  }

  JetParameters::JetParameters(Candidate& jet):
    pt(jet.Momentum.Pt()),
    eta(jet.Momentum.Eta()),
    flavor(simple_flavor(jet.Flavor))
  {
  }
  HighLevelTracking::HighLevelTracking(const ::HighLevelTracking& hlTrk):
    track_2_d0_significance(hlTrk.track2d0sig),
    track_3_d0_significance(hlTrk.track3d0sig),
    track_2_z0_significance(hlTrk.track2z0sig),
    track_3_z0_significance(hlTrk.track3z0sig),
    n_tracks_over_d0_threshold(hlTrk.tracksOverIpThreshold),
    jet_prob(hlTrk.jetProb),
    jet_width_eta(hlTrk.jetWidthEta),
    jet_width_phi(hlTrk.jetWidthPhi)
  {
  }
  HighLevelSecondaryVertex::HighLevelSecondaryVertex(
    const ::HighLevelSvx& hlSvx):
    vertex_significance(hlSvx.Lsig),
    n_secondary_vertices(hlSvx.NVertex),
    n_secondary_vertex_tracks(hlSvx.NTracks),
    delta_r_vertex(hlSvx.DrJet),
    vertex_mass(hlSvx.Mass),
    vertex_energy_fraction(hlSvx.EnergyFraction)
  {
  }

  HighLevelJet::HighLevelJet(Candidate& jet):
    jet_parameters(jet),
    tracking(jet.hlTrk),
    vertex(jet.hlSvx)
  {
  }

  VertexTrack::VertexTrack(const SecondaryVertexTrack& tk):
    d0(tk.d0), z0(tk.z0),
    d0_uncertainty(tk.d0err), z0_uncertainty(tk.z0err),
    pt(tk.pt),
    delta_phi_jet(tk.dphi), delta_eta_jet(tk.deta)
  {
  }
  SecondaryVertex::SecondaryVertex(const ::SecondaryVertex& vx):
    mass(vx.mass),
    displacement(vx.Mag()),
    delta_eta_jet(vx.deta),
    delta_phi_jet(vx.dphi),
    displacement_significance(vx.Lsig)
  {
    for (const auto& track: vx.tracks_along_jet) {
      associated_tracks.push_back(track);
    }
  }
  MediumLevelJet::MediumLevelJet(Candidate& jet):
    jet_parameters(jet)
  {
    for (const auto& trk: jet.primaryVertexTracks) {
      primary_vertex_tracks.push_back(trk);
    }
    for (const auto& vx: jet.secondaryVertices) {
      secondary_vertices.push_back(SecondaryVertex(vx));
    }
  }

  // ____________________________________________________________________
  // HDF5 types

  // insering a compound type requires that `type(Class)` is defined
  // high level variables
  H5::CompType type(JetParameters) {
    H5::CompType out(sizeof(JetParameters));
    H5_INSERT(out, JetParameters, pt);
    H5_INSERT(out, JetParameters, eta);
    H5_INSERT(out, JetParameters, flavor);
    return out;
  }

  H5::CompType type(HighLevelTracking) {
    H5::CompType out(sizeof(HighLevelTracking));
    H5_INSERT(out, HighLevelTracking, track_2_d0_significance);
    H5_INSERT(out, HighLevelTracking, track_3_d0_significance);
    H5_INSERT(out, HighLevelTracking, track_2_z0_significance);
    H5_INSERT(out, HighLevelTracking, track_3_z0_significance);
    H5_INSERT(out, HighLevelTracking, n_tracks_over_d0_threshold);
    H5_INSERT(out, HighLevelTracking, jet_prob);
    H5_INSERT(out, HighLevelTracking, jet_width_eta);
    H5_INSERT(out, HighLevelTracking, jet_width_phi);
    return out;
  }

  H5::CompType type(HighLevelSecondaryVertex) {
    H5::CompType out(sizeof(HighLevelSecondaryVertex));
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_significance);
    H5_INSERT(out, HighLevelSecondaryVertex, n_secondary_vertices);
    H5_INSERT(out, HighLevelSecondaryVertex, n_secondary_vertex_tracks);
    H5_INSERT(out, HighLevelSecondaryVertex, delta_r_vertex);
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_mass);
    H5_INSERT(out, HighLevelSecondaryVertex, vertex_energy_fraction);
    return out;
  }
  H5::CompType type(HighLevelJet) {
    H5::CompType out(sizeof(HighLevelJet));
    H5_INSERT(out, HighLevelJet, jet_parameters);
    H5_INSERT(out, HighLevelJet, tracking);
    H5_INSERT(out, HighLevelJet, vertex);
    return out;
  }

  // medium-level variables
  H5::CompType type(VertexTrack) {
    H5::CompType out(sizeof(VertexTrack));
    H5_INSERT(out, VertexTrack, d0);
    H5_INSERT(out, VertexTrack, z0);
    H5_INSERT(out, VertexTrack, d0_uncertainty);
    H5_INSERT(out, VertexTrack, z0_uncertainty);
    H5_INSERT(out, VertexTrack, pt);
    H5_INSERT(out, VertexTrack, delta_eta_jet);
    H5_INSERT(out, VertexTrack, delta_phi_jet);
    return out;
  }
  H5::CompType type(SecondaryVertex) {
    H5::CompType out(sizeof(SecondaryVertex));
    H5_INSERT(out, SecondaryVertex, mass);
    H5_INSERT(out, SecondaryVertex, displacement);
    H5_INSERT(out, SecondaryVertex, delta_eta_jet);
    H5_INSERT(out, SecondaryVertex, delta_phi_jet);
    H5_INSERT(out, SecondaryVertex, displacement_significance);
    H5_INSERT(out, SecondaryVertex, associated_tracks);
    return out;
  }
  H5::CompType type(MediumLevelJet) {
    H5::CompType out(sizeof(MediumLevelJet));
    H5_INSERT(out, MediumLevelJet, jet_parameters);
    H5_INSERT(out, MediumLevelJet, primary_vertex_tracks);
    H5_INSERT(out, MediumLevelJet, secondary_vertices);
    return out;
  }
}

//------------------------------------------------------------------------------

void HDF5Writer::Finish()
{
  m_hl_jet_buffer->flush();
  m_hl_jet_buffer->close();

  m_ml_jet_buffer->flush();
  m_ml_jet_buffer->close();
}

//------------------------------------------------------------------------------

void HDF5Writer::Process()
{
  fItInputArray->Reset();
  Candidate* jet;
  while ((jet = static_cast<Candidate*>(fItInputArray->Next()))) {
    m_hl_jet_buffer->push_back(*jet);
    m_ml_jet_buffer->push_back(*jet);
  }
}

//------------------------------------------------------------------------------
