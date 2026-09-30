////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessScale.cpp
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
//  Description: scales the mesh based on user input.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/Interpreter/Interpreter.h>
#include <NekMesh/MeshElements/Element.h>

#include "ProcessScale.h"

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar::NekMesh
{

ModuleKey ProcessScale::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "scale"), ProcessScale::create,
    "Scales the full mesh based on an input coefficient 'scale' ");

ProcessScale::ProcessScale(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["scaleX"] = ConfigOption(
        false, "1.0", "Define the coefficient to scale the mesh with.");
    m_config["scaleY"] = ConfigOption(
        false, "1.0", "Define the coefficient to scale the mesh with.");
    m_config["scaleZ"] = ConfigOption(
        false, "1.0", "Define the coefficient to scale the mesh with.");
}

ProcessScale::~ProcessScale()
{
}

void ProcessScale::Process()
{
    NekDouble scaleCoeff_X = m_config["scaleX"].as<NekDouble>();
    NekDouble scaleCoeff_Y = m_config["scaleY"].as<NekDouble>();
    NekDouble scaleCoeff_Z = m_config["scaleZ"].as<NekDouble>();

    m_log(VERBOSE) << "Scaled in X = " << scaleCoeff_X << endl;
    m_log(VERBOSE) << "Scaled in Y = " << scaleCoeff_Y << endl;
    m_log(VERBOSE) << "Scaled in Z = " << scaleCoeff_Z << endl;

    // Vertices
    for (auto vertex : m_mesh->m_vertexSet)
    {
        vertex->m_x *= scaleCoeff_X;
        vertex->m_y *= scaleCoeff_Y;
        vertex->m_z *= scaleCoeff_Z;
    }

    // EdgesNodes
    for (auto edge : m_mesh->m_edgeSet)
    {
        for (auto node : edge->m_edgeNodes)
        {
            node->m_x *= scaleCoeff_X;
            node->m_y *= scaleCoeff_Y;
            node->m_z *= scaleCoeff_Z;
        }
    }

    // FacesNodes
    for (auto face : m_mesh->m_faceSet)
    {
        for (auto node : face->m_faceNodes)
        {
            node->m_x *= scaleCoeff_X;
            node->m_y *= scaleCoeff_Y;
            node->m_z *= scaleCoeff_Z;
        }
    }
}
} // namespace Nektar::NekMesh
