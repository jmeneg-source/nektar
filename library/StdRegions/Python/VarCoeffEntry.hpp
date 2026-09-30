///////////////////////////////////////////////////////////////////////////////
//
// File: VarCoeffEntry.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Python wrapper for VarCoeffEntry.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PYTHON_BASICUITLS_VARCOEFFENTRY_HPP
#define NEKTAR_LIBUTILITIES_PYTHON_BASICUITLS_VARCOEFFENTRY_HPP

#include <type_traits>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <StdRegions/StdRegions.hpp>

#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

namespace pybind11::detail
{

template <> struct type_caster<VarCoeffEntry>
{
public:
    PYBIND11_TYPE_CASTER(VarCoeffEntry, const_name("VarCoeffEntry"));

    bool load(handle src, bool)
    {
        if (!array_t<NekDouble>::check_(src))
        {
            return false;
        }

        auto arr =
            array_t<NekDouble, array::c_style | array::forcecast>::ensure(src);
        if (!arr)
        {
            return false;
        }

        if (arr.ndim() != 1)
        {
            return false;
        }

        Array<OneD, NekDouble> nekArray = py::cast<Array<OneD, NekDouble>>(src);
        value                           = VarCoeffEntry(nekArray);

        return true;
    }

    static handle cast(VarCoeffEntry const &src, return_value_policy, handle)
    {
        // Create a Python capsule to hold a pointer that contains a lightweight
        // copy of src. That way we guarantee Python will still have access to
        // the memory allocated inside arr even if arr is deallocated in C++.
        py::capsule c(new VarCoeffEntry(src), [](void *ptr) {
            VarCoeffEntry *tmp = (VarCoeffEntry *)ptr;
            delete tmp;
        });

        array_t<NekDouble> array({src.GetValue().size()}, {},
                                 src.GetValue().data(), c);

        return array.release();
    }
};

} // namespace pybind11::detail

#endif
