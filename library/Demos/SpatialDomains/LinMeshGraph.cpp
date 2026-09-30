///////////////////////////////////////////////////////////////////////////////
//
// File: LinMeshGraph.cpp
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
// Description: Small utility to test a linear mesh graph generated
// from an input mesh
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/CppCommandLine.hpp>
#include <LibUtilities/BasicUtils/Filesystem.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/LinearMeshGraph.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <SpatialDomains/MeshPartition.h>

#include <boost/format.hpp>
#include <iostream>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::SpatialDomains;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: LinMeshGraph nsplit useSimplex "
                  << " <xml file1> [xml file 2..n]" << std::endl;
        return 1;
    }

    int nSplit      = atoi(argv[1]);
    bool useSimplex = atoi(argv[2]);
    std::vector<std::string> filenames(argv + 3, argv + argc);

    filenames.insert(filenames.begin(), argv[0]);

    CppCommandLine cmd(filenames);

    LibUtilities::SessionReaderSharedPtr vSession =
        LibUtilities::SessionReader::CreateInstance(cmd.GetArgc(),
                                                    cmd.GetArgv());

    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraphIO::Read(vSession);

    std::map<int, std::pair<int, std::vector<int>>> coeffmap;

    SpatialDomains::LinearMeshGraph linear(graph);
    linear.CreateLinearGraph(nSplit, coeffmap, true, useSimplex);
    SpatialDomains::MeshGraphSharedPtr lingraph = linear.GetLinearGraph();

    std::string outname = vSession->GetSessionName() + "_LinMesh";

    outname += ".xml";
    auto lingraphIO =
        SpatialDomains::GetMeshGraphIOFactory().CreateInstance("Xml");
    lingraphIO->SetMeshGraph(lingraph);
    lingraphIO->WriteGeometry(outname, true);

    vSession->Finalise();

    return 0;
}
