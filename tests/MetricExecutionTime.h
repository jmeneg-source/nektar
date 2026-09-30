///////////////////////////////////////////////////////////////////////////////
//
// File: MetricExecutionTime.h
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
// Description: Implementation of the execution time metric. A test will fail
// if the execution time of the test falls outside of an assigned tolerance.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTS_METRICEXECUTIONTIME_H
#define NEKTAR_TESTS_METRICEXECUTIONTIME_H

#include <Metric.h>

#include <regex>
#include <vector>

namespace Nektar
{
/**
 * @brief Data structure for an execution time field value.
 */
struct MetricExecutionTimeFieldValue
{
    MetricExecutionTimeFieldValue()
    {
    }

    MetricExecutionTimeFieldValue(std::string value) : m_value(value)
    {
    }

    /// The value to match. Defaults to empty string.
    std::string m_value = "";
    /// Indicates whether the metric should be skipped. Defaults to false.
    bool m_skip = false;
    /// The tolerance to use for checking the execution time. Defaults to 5.0.
    double m_tolerance = 5.0;
};

/**
 * @brief Metric that finds the execution time in an output and tests it against
 * an accepted value and tolerance.
 */
class MetricExecutionTime : public Metric
{
public:
    ~MetricExecutionTime() override
    {
    }

    static MetricSharedPtr create(TiXmlElement *metric, bool generate)
    {
        return MetricSharedPtr(new MetricExecutionTime(metric, generate));
    }

    static std::string type;

protected:
    /// Regex used to match an execution time in a test output.
    std::regex m_regex;
    /// Stores each execution time found in the test output.
    MetricExecutionTimeFieldValue m_match;
    /// If true, use stderr for testing/generation instead of stdout.
    bool m_useStderr = false;

    MetricExecutionTime(TiXmlElement *metric, bool generate);

    bool v_Test(std::istream &pStdout, std::istream &pStderr) override;
    void v_Generate(std::istream &pStdout, std::istream &pStderr) override;
};
} // namespace Nektar

#endif
