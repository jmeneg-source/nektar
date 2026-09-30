///////////////////////////////////////////////////////////////////////////////
//
// File: NekFactory.hpp
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
// Description: Python convenience wrappers for NekFactory.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PYTHON_BASICUTILS_NEKFACTORY_HPP
#define NEKTAR_LIBUTILITIES_PYTHON_BASICUTILS_NEKFACTORY_HPP

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <functional>
#include <type_traits>
#include <utility>

//
// The following code is adapted from:
// https://stackoverflow.com/questions/21192659/variadic-templates-and-stdbind
//
// The objective is to define an alternative to std::placeholders that can be
// used in combination with variadic templates to wrap arbitrary functions with
// arguments defined through a parameter pack.
//

template <std::size_t> struct placeholder_template
{
};

namespace std
{
template <std::size_t N>
struct is_placeholder<placeholder_template<N>>
    : integral_constant<std::size_t, N + 1>
{
};
} // namespace std

/**
 * @brief Helper class to be used in NekFactory_Register. Stores the Python
 * function used to call the creator function.
 */
#pragma GCC visibility push(hidden)
template <class T> class NekFactoryRegisterHelper
{
public:
    NekFactoryRegisterHelper(py::object obj) : m_obj(obj)
    {
    }

    ~NekFactoryRegisterHelper() = default;

    /**
     * @brief Callback to Python instantiator function (typically a class
     * constructor). Assumes that all arguments can be converted to boost.python
     * through converter infrastructure.
     */
    template <class... Args> std::shared_ptr<T> create(Args... args)
    {
        // Create an object using the Python callable.
        py::object inst = m_obj(args...);

        // Assume that it returns an object of the appropriate type.
        return py::cast<std::shared_ptr<T>>(inst);
    }

protected:
    /// Python function that is used to construct objects.
    py::object m_obj;
};
#pragma GCC visibility pop

template <typename tFac> class NekFactory_Register;

/**
 * @brief Lightweight wrapper for the factory RegisterCreatorFunction, to
 * support the ability for Python subclasses of @tparam tBase to register
 * themselves with the NekFactory.
 *
 * This function wraps NekFactory::RegisterCreatorFunction. This function
 * expects a function pointer to a C++ object that will construct an instance of
 * @tparam tBase. In this case we therefore need to construct a function call
 * that will construct our Python object (which is a subclass of @tparam tBase),
 * and then pass this back to Boost.Python to give the Python object back.
 *
 * We have to do some indirection here to get this to work, but we can
 * achieve this with the following strategy:
 *
 * - Create a @c NekFactoryCapsuleDestructor object, which as an argument will
 *   store the Python class instance that will be instantiated from the Python
 *   side.
 * - Using std::bind, construct a function pointer to the helper's creation
 *   function, NekFactoryCapsuleDestructor<tBase>::create.
 * - Create a Python capsule that will contain the @c
 *   NekFactoryCapsuleDestructor<tBase> instance, and register this in the
 *   global namespace of the current module. This then ties the capsule to the
 *   lifetime of the module.
 */
template <template <typename, typename, typename...> typename tFac,
          typename tBase, typename tKey, typename... tParam>
class NekFactory_Register<tFac<tKey, tBase, tParam...>>
{
public:
    NekFactory_Register(tFac<tKey, tBase, tParam...> &fac) : m_fac(fac)
    {
    }

    void operator()(tKey &name, py::object &obj, std::string const &nameKey)
    {
        // Create a factory register helper, which will call the C++ function to
        // create the factory.
        NekFactoryRegisterHelper<tBase> *helper =
            new NekFactoryRegisterHelper<tBase>(obj);

        // Register this with the factory using std::bind to grab a function
        // pointer to that particular object's function. We use
        // placeholder_template to bind to some number of arguments.
        DoRegister(name, helper, std::index_sequence_for<tParam...>{});

        // Create a capsule that will be embedded in the __main__ namespace. So
        // deallocation will occur, but only once Python ends or the Python
        // module is deallocated.
        std::string key = "_" + nameKey;

        // Allocate a capsule to ensure memory is cleared up upon deallocation
        // (depends on Python 2 or 3).
        py::capsule capsule(helper, [](void *ptr) {
            NekFactoryRegisterHelper<tBase> *tmp =
                (NekFactoryRegisterHelper<tBase> *)ptr;
            delete tmp;
        });

        // Embed the capsule in __main__.
        py::globals()[key.c_str()] = capsule;
    }

private:
    /**
     * @brief Helper function to unpack parameter arguments into the factory's
     * register creation function using the #placeholder_template helper struct.
     */
    template <std::size_t... Is>
    void DoRegister(tKey &name, NekFactoryRegisterHelper<tBase> *helper,
                    std::integer_sequence<std::size_t, Is...>)
    {
        m_fac.RegisterCreatorFunction(
            name,
            std::bind(
                &NekFactoryRegisterHelper<tBase>::template create<tParam...>,
                helper, placeholder_template<Is>{}...));
    }

    /// Reference to the NekFactory.
    tFac<tKey, tBase, tParam...> &m_fac;
};

#endif
