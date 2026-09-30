///////////////////////////////////////////////////////////////////////////////
//
// File: TestFile.h
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
// Description: Encapsulation of test XML file.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTER_TESTFILE
#define NEKTAR_TESTER_TESTFILE

#include <boost/program_options.hpp>

#include <string>
#include <vector>

#include <TestData.h>
#include <tinyxml.h>

#include <LibUtilities/BasicUtils/Filesystem.hpp>

namespace po = boost::program_options;

namespace Nektar
{

/**
 * @brief The TestData class is responsible for parsing a test XML file and
 * storing the data.
 */

class TestFile
{
public:
    TestFile(const fs::path &pFilename, po::variables_map &pVm);
    TestFile(const TestFile &pSrc);

    void Parse(TiXmlDocument *pDoc);

    void SaveFile();

    std::vector<TestData *> GetTests()
    {
        return m_tests;
    }

private:
    po::variables_map m_cmdoptions;
    TiXmlDocument *m_doc;

    /// Vector of one or more tests to run.
    std::vector<TestData *> m_tests;
};
} // namespace Nektar

#endif
