///////////////////////////////////////////////////////////////////////////////
//
//  File: CppCommandLine.hpp
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
//  Description: Utility class to define C++ command line arguments.
//
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#ifndef NEKTAR_LIBUTILITIES_BASICUTILS_CPPCOMMANDLINE_HPP
#define NEKTAR_LIBUTILITIES_BASICUTILS_CPPCOMMANDLINE_HPP

namespace Nektar::LibUtilities
{

/**
 * @brief Helper structure to construct C++ command line `argc` and `argv`
 * variables from a C++ vector.
 *
 * This is a useful class when setting up a call to
 * SessionReader::CreateInstance, since arguments must be constructed to the
 * right memory addresses, otherwise MPI (specifically OpenMPI) often tends to
 * segfault.
 */
struct CppCommandLine
{
    CppCommandLine() = default;

    /**
     * @brief Constructor.
     *
     * @param argv  List of command line arguments.
     */
    CppCommandLine(std::vector<std::string> argv)
    {
        Setup(argv);
    }

    /**
     * @brief Destructor.
     */
    ~CppCommandLine()
    {
        if (m_argv == nullptr)
        {
            return;
        }

        // Only single pointer delete is required since storage is in m_buf.
        delete[] m_argv;
    }

    /**
     * @brief Returns the constructed `argv`.
     */
    char **GetArgv()
    {
        return m_argv;
    }

    /**
     * @brief Returns the constructed `argc`.
     */
    int GetArgc()
    {
        return m_argc;
    }

protected:
    void Setup(std::vector<std::string> &argv)
    {
        int i          = 0;
        size_t bufSize = 0;
        char *p;

        m_argc = argv.size();
        m_argv = new char *[m_argc + 1];

        // Create argc, argv to give to the session reader. Note that this needs
        // to be a contiguous block in memory, otherwise MPI (specifically
        // OpenMPI) will likely segfault.
        for (i = 0; i < m_argc; ++i)
        {
            bufSize += argv[i].size() + 1;
        }

        m_buf.resize(bufSize);
        for (i = 0, p = &m_buf[0]; i < m_argc; ++i)
        {
            std::string tmp = argv[i];
            std::copy(tmp.begin(), tmp.end(), p);
            p[tmp.size()] = '\0';
            m_argv[i]     = p;
            p += tmp.size() + 1;
        }

        m_argv[m_argc] = nullptr;
    }

private:
    /// Pointers for strings `argv`.
    char **m_argv = nullptr;
    /// Number of arguments `argc`.
    int m_argc = 0;
    /// Buffer for storage of the argument strings.
    std::vector<char> m_buf;
};

} // namespace Nektar::LibUtilities

#endif
