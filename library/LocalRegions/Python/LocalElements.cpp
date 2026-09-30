///////////////////////////////////////////////////////////////////////////////
//
// File: LocalElements.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Python wrapper for Expansion elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/TriExp.h>

using namespace Nektar;
using namespace Nektar::LocalRegions;

void export_LocalElements(py::module &m)
{
    py::class_<SegExp, Expansion, StdRegions::StdSegExp,
               std::shared_ptr<SegExp>>(m, "SegExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      SpatialDomains::SegGeom *>());
    py::class_<TriExp, Expansion, StdRegions::StdTriExp,
               std::shared_ptr<TriExp>>(m, "TriExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      SpatialDomains::TriGeom *>());
    py::class_<QuadExp, Expansion, StdRegions::StdQuadExp,
               std::shared_ptr<QuadExp>>(m, "QuadExp",
                                         py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      SpatialDomains::QuadGeom *>());
    py::class_<TetExp, Expansion, StdRegions::StdTetExp,
               std::shared_ptr<TetExp>>(m, "TetExp", py::multiple_inheritance())
        .def(py::init<
             const LibUtilities::BasisKey &, const LibUtilities::BasisKey &,
             const LibUtilities::BasisKey &, SpatialDomains::TetGeom *>());
    py::class_<PrismExp, Expansion, StdRegions::StdPrismExp,
               std::shared_ptr<PrismExp>>(m, "PrismExp",
                                          py::multiple_inheritance())
        .def(py::init<
             const LibUtilities::BasisKey &, const LibUtilities::BasisKey &,
             const LibUtilities::BasisKey &, SpatialDomains::PrismGeom *>());
    py::class_<PyrExp, Expansion, StdRegions::StdPyrExp,
               std::shared_ptr<PyrExp>>(m, "PyrExp", py::multiple_inheritance())
        .def(py::init<
             const LibUtilities::BasisKey &, const LibUtilities::BasisKey &,
             const LibUtilities::BasisKey &, SpatialDomains::PyrGeom *>());
    py::class_<HexExp, Expansion, StdRegions::StdHexExp,
               std::shared_ptr<HexExp>>(m, "HexExp", py::multiple_inheritance())
        .def(py::init<
             const LibUtilities::BasisKey &, const LibUtilities::BasisKey &,
             const LibUtilities::BasisKey &, SpatialDomains::HexGeom *>());
}
