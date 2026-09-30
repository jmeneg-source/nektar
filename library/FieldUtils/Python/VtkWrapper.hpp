///////////////////////////////////////////////////////////////////////////////
//
// File: VtkWrapper.hpp
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
// Description: Helper routines from extracting VTK object pointers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_FIELDUTILS_PYTHON_VTKWRAPPER_HPP
#define NEKTAR_FIELDUTILS_PYTHON_VTKWRAPPER_HPP

#include <vtkObjectBase.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

/*
 * Adapted from the code at:
 * https://www.paraview.org/Wiki/Example_from_and_to_python_converters
 *
 * and: https://github.com/EricCousineau-TRI/repro/tree/master/python/vtk_pybind
 */

namespace pybind11::detail
{

template <typename T>
struct type_caster<T, enable_if_t<std::is_base_of_v<vtkObjectBase, T>>>
{
private:
    // Store value as a pointer
    T *value;

public:
    static constexpr auto name = _<T>();

    bool load(handle src, bool)
    {
        std::string thisStr = "__this__";

        PyObject *obj = src.ptr();

        // first we need to get the __this__ attribute from the Python Object
        if (!PyObject_HasAttrString(obj, thisStr.c_str()))
        {
            return false;
        }

        PyObject *thisAttr = PyObject_GetAttrString(obj, thisStr.c_str());
        if (thisAttr == nullptr)
        {
            return false;
        }

        std::string str = PyUnicode_AsUTF8(thisAttr);

        if (str.size() < 21)
        {
            return false;
        }

        auto iter = str.find("_p_vtk");
        if (iter == std::string::npos)
        {
            return false;
        }

        auto className = str.substr(iter).find("vtk");
        if (className == std::string::npos)
        {
            return false;
        }

        long address             = stol(str.substr(1, 17));
        vtkObjectBase *vtkObject = (vtkObjectBase *)((void *)address);

        if (vtkObject->IsA(str.substr(className).c_str()))
        {
            value = vtkObject;
            return true;
        }

        return false;
    }

    static handle cast(const T *src, return_value_policy policy, handle parent)
    {
        if (src == nullptr)
        {
            return none().release();
        }

        std::ostringstream oss;
        oss << (vtkObjectBase *)src; // here don't get address
        std::string address_str = oss.str();

        try
        {
            py::object obj = py::module_::import("vtkmodules.vtkCommonCore")
                                 .attr("vtkObjectBase")(address_str);
            return obj.release();
        }
        catch (py::error_already_set &)
        {
            // Clear error to avoid potential failure to catch error in next
            // try-catch block.
            PyErr_Clear();
        }

        try
        {
            py::object obj =
                py::module_::import("vtk").attr("vtkObjectBase")(address_str);
            return obj.release();
        }
        catch (py::error_already_set &)
        {
            using NekError = Nektar::ErrorUtil::NekError;
            throw NekError("Unable to import VTK.");
        }

        // Should never get here!
        return none().release();
    }
};

/// VTK Pointer-like object - may be non-copyable.
template <typename Ptr> struct vtk_ptr_cast_only
{
protected:
    using Class             = intrinsic_t<decltype(*std::declval<Ptr>())>;
    using value_caster_type = type_caster<Class>;

public:
    static constexpr auto name = _<Class>();
    static handle cast(const Ptr &ptr, return_value_policy policy,
                       handle parent)
    {
        return value_caster_type::cast(*ptr, policy, parent);
        ;
    }
};

/// VTK Pointer-like object - copyable / movable.
template <typename Ptr>
struct vtk_ptr_cast_and_load : public vtk_ptr_cast_only<Ptr>
{
private:
    Ptr value;
    // N.B. Can't easily access base versions...
    using Class             = intrinsic_t<decltype(*std::declval<Ptr>())>;
    using value_caster_type = type_caster<Class>;

public:
    operator Ptr &()
    {
        return value;
    }
    // Does this even make sense in VTK?
    operator Ptr &&() &&
    {
        return std::move(value);
    }

    template <typename T_>
    using cast_op_type = pybind11::detail::movable_cast_op_type<T_>;

    bool load(handle src, bool convert)
    {
        value_caster_type value_caster;
        if (!value_caster.load(src, convert))
        {
            return false;
        }
        value = Ptr(value_caster.operator Class *());
        return true;
    }
};

template <typename Class>
struct type_caster<vtkSmartPointer<Class>>
    : public vtk_ptr_cast_and_load<vtkSmartPointer<Class>>
{
};

template <typename Class>
struct type_caster<vtkNew<Class>> : public vtk_ptr_cast_only<vtkNew<Class>>
{
};

} // namespace pybind11::detail

#endif
