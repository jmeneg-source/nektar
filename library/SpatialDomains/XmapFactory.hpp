///////////////////////////////////////////////////////////////////////////////
//
// File: XmapFactory.hpp
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
// Description: Simple factory pattern for creation of geometry xmap objects.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_XMAPFACTORY_HPP
#define NEKTAR_SPATIALDOMAINS_XMAPFACTORY_HPP

#include <functional>
#include <memory>
#include <type_traits>
#include <unordered_map>

namespace Nektar::SpatialDomains
{

/**
 * @brief A simple factory for Xmap objects that is based on the element type,
 * the basis and quadrature selection.
 */
template <class StdExp, int dim> class XmapFactory
{
    /// Key for the Xmap: an array of dimension dim which contains the basis
    /// keys for each coordinate direction.
    using key_t = std::array<LibUtilities::BasisKey, dim>;

    /// Hash for key_t
    struct KeyHash
    {
        /**
         * We base key_t on the idea that number of points and modes are almost
         * certainly < 2^8. Therefore total hash size is at most given by dim *
         * 4 * 2 <= 48 bits, which comfortably fits in a 64-bit size_t.
         */
        size_t operator()(key_t const &key) const
        {

            size_t hash = 0;
            for (int d = 0; d < dim; ++d)
            {
                // This loop should be unrolled by the compiler.
                hash = (hash << 8) | (((key[d].GetNumPoints() & 0xFF) << 8) |
                                      (key[d].GetNumPoints() & 0xFF));
            }

            return hash;
        }
    };

    struct KeyEqual
    {
        bool operator()(key_t const &a, key_t const &b) const
        {
            bool equal = true;
            for (int d = 0; d < dim; ++d)
            {
                equal = equal && (a[d].GetNumPoints() == b[d].GetNumPoints() &&
                                  a[d].GetNumModes() == b[d].GetNumModes());
            }
            return equal;
        }
    };

    template <int d, typename std::enable_if<d == 1, int>::type = 0>
    std::shared_ptr<StdExp> CreateStdExp(const key_t &args)
    {
        return std::make_shared<StdExp>(args[0]);
    }

    template <int d, typename std::enable_if<d == 2, int>::type = 0>
    std::shared_ptr<StdExp> CreateStdExp(const key_t &args)
    {
        return std::make_shared<StdExp>(args[0], args[1]);
    }

    template <int d, typename std::enable_if<d == 3, int>::type = 0>
    std::shared_ptr<StdExp> CreateStdExp(const key_t &args)
    {
        return std::make_shared<StdExp>(args[0], args[1], args[2]);
    }

public:
    XmapFactory() = default;

    /**
     * @brief Returns (and creates, if necessary) a standard expansion
     * corresponding to @tparam StdExp and the supplied @p basis.
     */
    StdRegions::StdExpansionSharedPtr CreateInstance(const key_t &basis)
    {
        auto it = m_xmaps.find(basis);

        if (it != m_xmaps.end())
        {
            return it->second;
        }

        std::shared_ptr<StdExp> xmap = CreateStdExp<dim>(basis);
        m_xmaps[basis]               = xmap;

        return xmap;
    }

private:
    /// Storage for created xmap objects.
    std::unordered_map<key_t, StdRegions::StdExpansionSharedPtr, KeyHash,
                       KeyEqual>
        m_xmaps;
};

} // namespace Nektar::SpatialDomains

#endif
