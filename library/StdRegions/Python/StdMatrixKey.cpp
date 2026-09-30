////////////////////////////////////////////////////////////////////////////////
//
//  File: StdMatrixKey.cpp
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
//  Description: Python wrapper for StdMatrixKey.
//
////////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdMatrixKey.h>

#include <vector>

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <StdRegions/Python/VarCoeffEntry.hpp>

using namespace Nektar;
using namespace Nektar::StdRegions;

PYBIND11_MAKE_OPAQUE(ConstFactorMap)
PYBIND11_MAKE_OPAQUE(VarCoeffMap)

std::unique_ptr<StdMatrixKey> StdMatrixKey_Init(
    const MatrixType matType, const LibUtilities::ShapeType shapeType,
    const StdExpansionSharedPtr exp, const ConstFactorMap &constFactorMap,
    const VarCoeffMap &varCoeffMap)
{
    return std::make_unique<StdMatrixKey>(matType, shapeType, *exp,
                                          constFactorMap, varCoeffMap);
}

/**
 * @brief Export for StdMatrixKey enumeration.
 */
void export_StdMatrixKey(py::module &m)
{
    NEKPY_WRAP_ENUM(m, MatrixType, MatrixTypeMap);
    NEKPY_WRAP_ENUM(m, ConstFactorType, ConstFactorTypeMap);
    NEKPY_WRAP_ENUM(m, VarCoeffType, VarCoeffTypeMap);

    // Wrapper for constant factor map.
    py::bind_map<ConstFactorMap>(m, "ConstFactorMap");

    // Wrapper for variable coefficients map.
    py::bind_map<VarCoeffMap>(m, "VarCoeffMap");

    py::class_<StdMatrixKey>(m, "StdMatrixKey")
        .def(py::init(&StdMatrixKey_Init), py::arg("matType"),
             py::arg("shapeType"), py::arg("exp"),
             py::arg("constFactorMap") = NullConstFactorMap,
             py::arg("varCoeffMap")    = NullVarCoeffMap)

        .def("GetMatrixType", &StdMatrixKey::GetMatrixType)
        .def("GetShapeType", &StdMatrixKey::GetShapeType)
        .def("GetNcoeffs", &StdMatrixKey::GetNcoeffs)
        .def("GetBasis", &StdMatrixKey::GetBasis);
}
