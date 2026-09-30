///////////////////////////////////////////////////////////////////////////////
//
// File: Filter.cpp
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
// Description: Base class for filters.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/Filter.h>

namespace Nektar::SolverUtils
{
FilterFactory &GetFilterFactory()
{
    static FilterFactory instance;
    return instance;
}

Filter::Filter(const LibUtilities::SessionReaderSharedPtr &pSession,
               const std::shared_ptr<EquationSystem> &pEquation)
    : m_session(pSession), m_equ(pEquation)
{
}

std::string Filter::v_SetupOutput(const std::string ext,
                                  const ParamMap &pParams)
{
    // get the backup switch
    bool backup = m_session->GetBackups();

    // determinte the root
    bool root = m_session->GetComm()->TreatAsRankZero();

    std::string outname;

    // set up the base output file name
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        outname = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        outname = it->second;
    }

    // remove any aextension from the input name and add the correct one
    outname = fs::path(outname).replace_extension("").string() + ext;
    fs::path specPath(outname);
    // filters output are handled by root (rank 0) only. check if the path is
    // already exists
    if (backup && root && fs::exists(specPath))
    {
        // rename. foo/output-name.ext -> foo/output-name.bak0.ext and in case
        // foo/output-name.bak0.chk already exists,
        // foo/output-name.bak0.ext ->foo/output-name.bak1.chk
        fs::path bakPath = specPath;
        int cnt          = 0;
        while (fs::exists(bakPath))
        {
            bakPath = specPath.parent_path();
            // this is for parallel in time
            bakPath += specPath.stem();
            bakPath += fs::path(".bak" + std::to_string(cnt++));
            bakPath += specPath.extension();
        }
        std::cout << "renaming " << specPath << " -> " << bakPath << std::endl;
        try
        {
            fs::rename(specPath, bakPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(e.code() == std::errc::no_such_file_or_directory,
                     "Filesystem error: " + std::string(e.what()));
        }
    }
    // return the output file name
    return LibUtilities::PortablePath(specPath);
}

std::string Filter::v_SetupOutput(const std::string ext,
                                  const std::string inname)
{
    // get the backup switch
    bool backup = m_session->GetBackups();

    // determinte the root
    bool root = m_session->GetComm()->TreatAsRankZero();

    // for the history point filter we have already established the file name
    std::string outname = inname;

    // remove any aextension from the input name and add the correct one
    outname = fs::path(outname).replace_extension("").string() + ext;
    // Path to output
    fs::path specPath(outname);

    // filters output are handled by root (rank 0) only. check if the path is
    // already exists
    if (backup && root && fs::exists(specPath))
    {
        // rename. foo/output-name.ext -> foo/output-name.bak0.ext and in case
        // foo/output-name.bak0.chk already exists,
        // foo/output-name.bak0.ext ->foo/output-name.bak1.chk
        //
        fs::path bakPath = specPath;
        int cnt          = 0;
        while (fs::exists(bakPath))
        {
            bakPath = specPath.parent_path();
            // this is for parallel in time
            bakPath += specPath.stem();
            bakPath += fs::path(".bak" + std::to_string(cnt++));
            bakPath += specPath.extension();
        }
        std::cout << "renaming " << specPath << " -> " << bakPath << std::endl;
        try
        {
            fs::rename(specPath, bakPath);
        }
        catch (fs::filesystem_error &e)
        {
            ASSERTL0(e.code() == std::errc::no_such_file_or_directory,
                     "Filesystem error: " + std::string(e.what()));
        }
    }

    return LibUtilities::PortablePath(specPath);
}

} // namespace Nektar::SolverUtils
