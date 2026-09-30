///////////////////////////////////////////////////////////////////////////////
//
// File: TestFile.cpp
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

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>

#include <TestException.hpp>
#include <TestFile.h>

using namespace std;

namespace Nektar
{

/**
 * @brief TestData constructor.
 *
 * The class is constructed with the path to the test XML file and a
 * `po::variables_map` object containing the command-line options passed to the
 * program.
 *
 * @param pFilename
 * @param pVm
 */

TestFile::TestFile(const fs::path &pFilename, po::variables_map &pVm)
    : m_cmdoptions(pVm)
{
    // Process test file format.
    m_doc = new TiXmlDocument(pFilename.string().c_str());

    bool loadOkay = m_doc->LoadFile();

    ASSERTL0(loadOkay,
             "Failed to load test definition file: " + pFilename.string() +
                 "\n" + string(m_doc->ErrorDesc()));

    Parse(m_doc);
}

TestFile::TestFile(const TestFile &pSrc)
{
    boost::ignore_unused(pSrc);
}

/// Parse the test file and populate member variables for the test.
void TestFile::Parse(TiXmlDocument *pDoc)
{
    TiXmlHandle handle(pDoc);
    TiXmlElement *tmp, *testElement;

    // Check to see whether the .tst file includes multiple <test> elements
    // contained within the root document.
    tmp = handle.FirstChildElement("tests").Element();
    if (tmp)
    {
        testElement = tmp->FirstChildElement("test");
    }
    else
    {
        testElement = handle.FirstChildElement("test").Element();
    }

    ASSERTL0(testElement, "Cannot find 'test' or 'tests' root element.");

    while (testElement)
    {
        m_tests.push_back(new TestData(testElement, m_cmdoptions));
        testElement = testElement->NextSiblingElement("test");
    }
}

void TestFile::SaveFile()
{
    m_doc->SaveFile();
}

} // namespace Nektar
