///////////////////////////////////////////////////////////////////////////////
//
// File: Metric.h
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
// Description: Definition of the metric base class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRIC_H
#define NEKTAR_TESTS_METRIC_H

#include <map>
#include <memory>
#include <string>
#include <tinyxml.h>

#include <TestException.hpp>

namespace Nektar
{
/**
 * @brief Check to see whether the given string @p s is empty (or null).
 */
inline bool EmptyString(const char *s)
{
    if (!s)
    {
        return true;
    }
    return std::string(s) == "";
}

/**
 * @brief Base class for all metrics.
 * Metric represents a test metric that can be used to evaluate the
 * functionality or performance of a Nektar++ executable.
 */
class Metric
{
public:
    Metric(TiXmlElement *metric, bool generate);

    virtual ~Metric() = default;

    bool Test(std::istream &pStdout, std::istream &pStderr);
    void Generate(std::istream &pStdout, std::istream &pStderr);
    /// Return metric type
    std::string GetType()
    {
        return m_type;
    }
    /// Return metric ID
    int GetID()
    {
        return m_id;
    }
    /// Return whether this metric supports averaging results from multiple
    /// runs.
    bool SupportsAverage() const
    {
        return m_average;
    }

protected:
    /// Stores the ID of this metric.
    int m_id;
    /// Stores the type of this metric (uppercase).
    std::string m_type;
    /// Determines whether to generate this metric or not.
    bool m_generate;
    /// Indicates whether a metric supports averaging results from multiple
    /// runs.
    bool m_average = false;
    /// Pointer to XML structure containing metric definition.
    TiXmlElement *m_metric;

    /**
     * @brief Virtual function to test the metric. Should be redefined in
     * derived classes.
     *
     * @param pStdout Reference to test output stream.
     * @param pStderr Reference to test error stream.
     * @return \p true if the test passes, \p false otherwise.
     */
    virtual bool v_Test(std::istream &pStdout, std::istream &pStderr) = 0;

    /**
     * @brief Virtual function to generate the metric. Should be redefined in
     * derived classes.
     *
     * @param pStdout Reference to test output stream.
     * @param pSrderr Reference to test error stream.
     */
    virtual void v_Generate(std::istream &pStdout, std::istream &pSrderr) = 0;
};

/// A shared pointer to an EquationSystem object
typedef std::shared_ptr<Metric> MetricSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from the
/// Advection class.
class MetricFactory
{
public:
    typedef MetricSharedPtr (*CreatorFunction)(TiXmlElement *, bool);

    std::string RegisterCreatorFunction(std::string key, CreatorFunction func)
    {
        m_map[key] = func;
        return key;
    }

    MetricSharedPtr CreateInstance(std::string key, TiXmlElement *elmt,
                                   bool generate)
    {
        return m_map[key](elmt, generate);
    }

private:
    std::map<std::string, CreatorFunction> m_map;
};

MetricFactory &GetMetricFactory();
} // namespace Nektar

#endif
