///////////////////////////////////////////////////////////////////////////////
//
// File: PEMFCSolver.cpp
//
// PEM fuel cell solver driver. Part of the nektar++_fuel_cells project
// (EPUSP / RCGI). Equation systems provided: PEMFCSpecies (stage 1).
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Driver.h>
#include <SpatialDomains/MeshGraphIO.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
    string vDriverModule;
    DriverSharedPtr drv;

    try
    {
        // Create session reader and mesh graph.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        graph   = SpatialDomains::MeshGraphIO::Read(session);

        // Create and execute driver.
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);
        drv->Execute();

        session->Finalise();
    }
    catch (const std::runtime_error &e)
    {
        return 1;
    }
    catch (const std::string &eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
