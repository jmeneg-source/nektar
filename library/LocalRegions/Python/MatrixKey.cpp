////////////////////////////////////////////////////////////////////////////////
//
//  File: MatrixKey.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Python wrapper for MatrixKey.
//
////////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/MatrixKey.h>
#include <vector>

#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::LocalRegions;

PYBIND11_MAKE_OPAQUE(StdRegions::ConstFactorMap)
PYBIND11_MAKE_OPAQUE(StdRegions::VarCoeffMap)

std::unique_ptr<MatrixKey> MatrixKey_Init(
    const StdRegions::MatrixType matType,
    const LibUtilities::ShapeType shapeType,
    const StdRegions::StdExpansionSharedPtr exp,
    const StdRegions::ConstFactorMap &constFactorMap,
    const StdRegions::VarCoeffMap &varCoeffMap)
{
    return std::make_unique<MatrixKey>(matType, shapeType, *exp, constFactorMap,
                                       varCoeffMap);
}

/**
 * @brief Export for MatrixKey enumeration.
 */
void export_MatrixKey(py::module &m)
{
    py::class_<MatrixKey, StdRegions::StdMatrixKey>(m, "MatrixKey")
        .def(py::init<>(&MatrixKey_Init), py::arg("matType"),
             py::arg("shapeType"), py::arg("exp"),
             py::arg("constFactorMap") = StdRegions::NullConstFactorMap,
             py::arg("varCoeffMap")    = StdRegions::NullVarCoeffMap);
}
