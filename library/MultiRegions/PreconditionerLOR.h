///////////////////////////////////////////////////////////////////////////////
//
// File: PreconditionerLOR.h
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
// Description: PreconditionerLOR header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOR_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOR_H

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/ContField.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/Preconditioner.h>

namespace Nektar::MultiRegions
{

class PreconditionerLOR;

/**
 * @brief Low-order representation (LOR) preconditioner.
 *
 * This preconditioner works by taking a high-order element and splitting into
 * P1 elements that 'join the dots' of a quadrature distribution on the element.
 */
class PreconditionerLOR : public Preconditioner
{
public:
    /// Creates an instance of this class
    static PreconditionerSharedPtr create(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
    {
        PreconditionerSharedPtr p =
            MemoryManager<PreconditionerLOR>::AllocateSharedPtr(plinsys,
                                                                pLocToGloMap);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    /// Constructor.
    MULTI_REGIONS_EXPORT PreconditionerLOR(
        const std::shared_ptr<GlobalLinSys> &plinsys,
        const AssemblyMapSharedPtr &pLocToGloMap);

    MULTI_REGIONS_EXPORT ~PreconditionerLOR() override
    {
    }

protected:
    void v_InitObject() override;

    void v_DoPreconditioner(const Array<OneD, NekDouble> &pInput,
                            Array<OneD, NekDouble> &pOutput,
                            const bool &isLocal = false) override;

    void v_BuildPreconditioner() override;

private:
    static std::string lookupIds[];
    static std::string def;

    /// Linear solver for LOR system
    GlobalLinSysSharedPtr m_LORLinSys;
    /// Solution type for LOR space, DirectFull, IterativeFull, PETScFull
    std::string m_slvType;
    /// Flag that determines whether to use equispaced point distributions or
    /// GLL.
    bool m_equiSpaced = false;
    /// Use a simplex distribution in the LOR system so that e.g. quadrilaterals
    /// are split into triangles.
    bool m_useSimplex = false;
    //// Number of points to use for splitting high-order elements, usually
    //// nummodes-1.
    unsigned int m_nsplit;
    /// MeshGraph for the LOR system
    SpatialDomains::MeshGraphSharedPtr m_lor_graph;
    /// Field object for the LOR system
    MultiRegions::ContFieldSharedPtr m_lor_field;
    /// Number of Dirichlet values on original mesh
    unsigned int m_nDir;
    /// Map from high-order node numbers to LOR representation
    Array<OneD, int> m_ho2lor;
    /// Inverse of multiplicity of each global DOF in the high-order mesh
    Array<OneD, NekDouble> m_invMultiplicity;
    /// Inverse of multiplicity of each global DOF in the LOR mesh
    Array<OneD, NekDouble> m_invLinMeshMultiplicity;
    /// Whether to be verbose.
    bool m_verboseIter;

    void CreateInvMultiplicity(void);
};

} // namespace Nektar::MultiRegions
#endif
