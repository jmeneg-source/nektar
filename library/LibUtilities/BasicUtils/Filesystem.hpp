///////////////////////////////////////////////////////////////////////////////
//
// File: Filesystem.hpp
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
// Description: Header file to simplify use of <filesystem>
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_FILESYSTEM_HPP
#define NEKTAR_LIBUTILITIES_FILESYSTEM_HPP

#ifdef NEKTAR_USE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
namespace fs     = boost::filesystem;
namespace fserrc = boost::system::errc;
#else
#ifdef NEKTAR_USE_GNU_FS_EXPERIMENTAL
#include <experimental/filesystem>
#include <system_error>
namespace fs = std::experimental::filesystem;
using fserrc = std::errc;
#else
#include <filesystem>
#include <system_error>
namespace fs = std::filesystem;
using fserrc = std::errc;
#endif
#endif

#include <algorithm>
#include <chrono>
#include <random>

namespace Nektar::LibUtilities
{

/**
 * \brief create portable path on different platforms for std::filesystem path.
 */
static inline std::string PortablePath(const fs::path &path)
{
    return fs::path(path).make_preferred().string();
}

/**
 * @brief Create a unique (random) path, based on an input stem string. The
 * returned string is a filename or directory that does not exist.
 *
 * Given @p specFileStem as a parameter, this returns a string in the form
 * `tmp_<stem>_abcdef` where the final 6 characters are random characters
 * and digits.
 */
static inline fs::path UniquePath(std::string specFileStem)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    fs::path tmp;

    do
    {
        std::string chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                            "abcdefghijklmnopqrstuvwxyz";

        std::shuffle(chars.begin(), chars.end(), generator);
        tmp = fs::path("tmp_" + specFileStem + "_" + chars.substr(0, 6));
    } while (fs::exists(tmp));

    return tmp;
}

} // namespace Nektar::LibUtilities

#endif
