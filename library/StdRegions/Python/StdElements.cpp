///////////////////////////////////////////////////////////////////////////////
//
// File: StdElements.cpp
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
// Description: Python wrapper for all Nektar++ elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdTriExp.h>

#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>

using namespace Nektar;
using namespace Nektar::StdRegions;

void export_StdElements(py::module &m)
{
    py::class_<StdPointExp, StdExpansion, std::shared_ptr<StdPointExp>>(
        m, "StdPointExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &>());
    py::class_<StdSegExp, StdExpansion, std::shared_ptr<StdSegExp>>(
        m, "StdSegExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &>());
    py::class_<StdQuadExp, StdExpansion, std::shared_ptr<StdQuadExp>>(
        m, "StdQuadExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
    py::class_<StdTriExp, StdExpansion, std::shared_ptr<StdTriExp>>(
        m, "StdTriExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
    py::class_<StdTetExp, StdExpansion, std::shared_ptr<StdTetExp>>(
        m, "StdTetExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
    py::class_<StdPrismExp, StdExpansion, std::shared_ptr<StdPrismExp>>(
        m, "StdPrismExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
    py::class_<StdPyrExp, StdExpansion, std::shared_ptr<StdPyrExp>>(
        m, "StdPyrExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
    py::class_<StdHexExp, StdExpansion, std::shared_ptr<StdHexExp>>(
        m, "StdHexExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
}
