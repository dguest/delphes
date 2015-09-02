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

#ifndef HDF5Writer_h
#define HDF5Writer_h

/** \class HDF5Writer
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */



class TObjArray;
class DelphesFormula;

#ifndef __CINT__

#include "external/h5/OneDimBuffer.hh"
#include "external/h5/h5container.hh"

#include "classes/DelphesModule.h"
#include "h5/bork.hh"

#include "H5Cpp.h"

namespace out {
  struct Vertex {
    double mass;
    double dr_jet;
  };

  struct Jet {
    double pt;
    double eta;
    h5::vector<Vertex> vertices;
  };
  H5::CompType type(Jet);
}

#else  // CINT include dummy

namespace H5 {
  class H5File;
}

#endif

class HDF5Writer: public DelphesModule
{
public:

  HDF5Writer();
  ~HDF5Writer();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  H5::H5File* m_out_file;

#ifndef __CINT__
  OneDimBuffer<out::Jet>* m_jet_buffer;
#endif

  ClassDef(HDF5Writer, 1)
};

#endif
