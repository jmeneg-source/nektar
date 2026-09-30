///////////////////////////////////////////////////////////////////////////////
//
// File: TestTetrahedron.cpp
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
/// The above copyright notice and this permission notice shall be included
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <NekMesh/MeshElements/Element.h>
#include <NekMesh/NekMeshDeclspec.h>

#include <boost/test/unit_test.hpp>

namespace Nektar::NekMeshTetrahedronUnitTest
{
using namespace Nektar::NekMesh;

BOOST_AUTO_TEST_CASE(TestConstruction)
{
    // Setup element config
    ElmtConfig el_conf(LibUtilities::eTetrahedron, 1, false, false, false);

    // Setup element node list -- order matters!
    std::vector<NodeSharedPtr> node_set(4);
    node_set[0] = std::make_shared<Node>(0, -1., 0., 0.);
    node_set[1] = std::make_shared<Node>(1, 1., 0., 0.);
    node_set[2] = std::make_shared<Node>(2, 0., 1., 0.);
    node_set[3] = std::make_shared<Node>(3, 0., 0.5, 1.);

    std::vector<int> tags;
    tags.push_back(0);

    ElementSharedPtr TestTet = GetElementFactory().CreateInstance(
        LibUtilities::eTetrahedron, el_conf, node_set, tags);

    BOOST_CHECK_EQUAL(TestTet->GetVertex(3)->m_x, 0.);
    BOOST_CHECK_EQUAL(TestTet->GetVertex(3)->m_y, 0.5);
    BOOST_CHECK_EQUAL(TestTet->GetVertex(3)->m_z, 1.);
}

} // namespace Nektar::NekMeshTetrahedronUnitTest
