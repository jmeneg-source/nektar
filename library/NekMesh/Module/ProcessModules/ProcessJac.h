////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSJAC
#define UTILITIES_NEKMESH_PROCESSJAC

#include <NekMesh/Module/Module.h>

namespace Nektar::NekMesh
{

/**
 * @brief This processing module calculates the Jacobian of elements
 * using %SpatialDomains::GeomFactors and the %Element::GetGeom
 * method. For now it simply prints a list of elements which have
 * negative Jacobian.
 */
class ProcessJac : public NekMesh::ProcessModule
{
public:
    static std::shared_ptr<Module> create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<ProcessJac>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className;

    ProcessJac(NekMesh::MeshSharedPtr m);
    ~ProcessJac() override;

    /// Write mesh to output file.
    void Process() override;

    std::string GetModuleName() override
    {
        return "ProcessJac";
    }
    /**
     * @brief Outputing space in the jac.txt using std::setw().
     */
    enum
    {
        ElementID        = 10,
        Jac              = 15,
        type             = 15,
        Boundary_edge_ID = 20,
        VertexID         = 15,
        Boundary_face_ID = 20,
        EdgeID           = 10,
        CoordX           = 15,
        CoordY           = 15,
        CoordZ           = 15,
        CompositeName    = 15,
        CompositeID      = 15,
        NeighborElID     = 30
    };
    typedef struct
    {
        NekDouble Jac;
        ElementSharedPtr El;
    } element_reorder;

private:
    /// The maximum value shown on the on-screen histogram. Defined by the user
    /// in the first value of "histo".
    NekDouble m_histo_max_value;
    /// The number of bins with positive Jacobians of the histogram shown on the
    /// screen. Defined by the user in the second value of "histo".
    int m_ScalePos;
    /// The number of bins with negative Jacobians of the histogram shown on the
    /// screen. Defined by the user in the third value of "histo"
    int m_ScaleNeg;
    /// The interval of the bin of the histogram shown on the screen, which is
    /// calculated by m_histo_max_value/m_ScalePos
    NekDouble m_interval;

    /**
     * @brief Get the coordinates of the boundary vertices of the boundary
     * face/edge and export them to the text file.
     *
     * @param el SharedPointer to the evaluated element.
     * @param output_file The file where information is writtern to.
     * @param detail True also gets and writes the composite name the element
     * boundary is on.
     */
    void GetBoundaryCoordinate(const ElementSharedPtr &el,
                               std::ofstream &output_file, bool detail);
    /**
     * @brief Get the composite name and ID of a boundary element, which can be
     * further exported into a text file.
     *
     * @param el SharedPointer to the evaluated element.
     * @param CompositeNamed The composite name of the boundary element.
     * @param CompositeID The ID of the composite where the boundary element
     * lies on.
     * @return False if the imported mesh does not set name for
     * composites.
     */

    bool GetCompositeName(const ElementSharedPtr &el,
                          std::string &CompositeNamed,
                          unsigned int &CompositeID);
    /**
     * @brief  Output histogram of scaled Jacobian on the screen.
     *
     * @param histo Array which contains the number of elements in each interval
     * between 0.0 and 1.0.
     */
    void GetHistogram(const Array<OneD, NekDouble> &histo);

    /**
     * @brief Check quality based on scaled Jacobians. The scaled Jacobian is
     * between 0.0 to 1.0, and is divided to 10 columns (by default, this number
     * is the same as that in the function GetHistogram), this function shows
     * the number of elements in each columns
     *
     * @param QuaHisto Array which contains the number of elements in each
     * interval between 0.0 and 1.0.
     */
    void Qualitycheck(Array<OneD, NekDouble> QuaHisto);
};
} // namespace Nektar::NekMesh

#endif
