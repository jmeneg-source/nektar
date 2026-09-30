///////////////////////////////////////////////////////////////////////////////
//
// File: MetricExecutionTime.cpp
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

#include <MetricExecutionTime.h>

#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace Nektar
{
std::string MetricExecutionTime::type =
    GetMetricFactory().RegisterCreatorFunction("EXECUTIONTIME",
                                               MetricExecutionTime::create);

/**
 * @brief Construct a new MetricExecutionTime object.
 *
 * If the metric is not a derived class, then this constructor will parse the
 * regular expression from the metric's element in the test file, or use a
 * default value if this doesn't exist. It will also see if an expected
 * execution time and tolerance have been provided for the current computer's
 * hostname. If this is the case, these will be used by \p v_Test. If not, then
 * a skip flag \p m_skip will be set to \p true and the test will automatically
 * pass.
 *
 * @param metric
 * @param generate
 */
MetricExecutionTime::MetricExecutionTime(TiXmlElement *metric, bool generate)
    : Metric(metric, generate)
{
    // If we are a derived class, do nothing
    if (m_type != "EXECUTIONTIME")
    {
        return;
    }

    // Parse Regex expression
    TiXmlElement *regex = metric->FirstChildElement("regex");

    if (regex)
    {
        ASSERTL0(regex->GetText(), "Failed to read regex element.");
        m_regex = regex->GetText();
    }
    else
    {
        // Set the regex to a default value.
        m_regex = R"(^.*Total Computation Time\s*=\s*(\d+\.?\d*).*)";
    }

    // Inform the Tester this metric supports averaging data from multiple runs.
    m_average = true;
    m_match   = MetricExecutionTimeFieldValue("nil");

    // Parse matching values if not generating.
    if (m_generate)
    {
        return;
    }

    TiXmlElement *value = metric->FirstChildElement("value");
    ASSERTL0(value || m_generate, "Missing value tag for metric.");

    bool hostnameMatch = false;
    while (value)
    {
        ASSERTL0(value->Attribute("tolerance"),
                 "Missing tolerance in execution time metric.");
        ASSERTL0(value->Attribute("hostname"),
                 "Missing hostname in execution time metric.");
        ASSERTL0(!EmptyString(value->GetText()),
                 "Missing value in execution time metric.");

        // Only use the value if it matches the runner's hostname.
        if (value->Attribute("hostname") == boost::asio::ip::host_name())
        {
            MetricExecutionTimeFieldValue val;
            val.m_value     = value->GetText();
            val.m_tolerance = atof(value->Attribute("tolerance"));
            m_match         = val;

            hostnameMatch = true;

            // No use finding any more values
            break;
        }

        value = value->NextSiblingElement("value");
    }

    // If no hostname match was found, pass a skip flag to the test.
    if (!hostnameMatch)
    {
        cerr << "WARNING: No execution time value provided for host "
             << boost::asio::ip::host_name() << ". Skipping metric." << endl;
        m_match.m_skip = true;
    }
}

/**
 * @brief Test output against a regular expression and its expected value.
 *
 * This function will find the execution time for each test run in the output,
 * using \b m_regex (provided in the test file). It will then average these
 * times and compare them to the expected execution time for the current
 * computer's hostname (this behaviour is defined in the constructor). If the
 * average execution time falls outside of the tolerance then the test will
 * fail.
 *
 * If no execution time value is provided for the current computer's hostname,
 * \b m_skip will be set to \b true and the test will automatically pass.
 *
 * @param pStdout Reference to test output stream.
 * @param pStderr Reference to test error stream.
 * @return \b true if the test passes, \b false otherwise.
 */
bool MetricExecutionTime::v_Test(istream &pStdout, istream &pStderr)
{
    boost::ignore_unused(pStdout, pStderr);

    bool success = true;

    // Select istream to use.
    istream &is = m_useStderr ? pStderr : pStdout;

    std::smatch matches;
    string line;
    // Vector of execution times found in the output.
    vector<double> times;
    bool matched = false;

    // Process output file line by line searching for regex matches
    while (getline(is, line))
    {
        // Test to see if we have a match on this line.
        if (std::regex_match(line, matches, m_regex))
        {
            // If no matches are found then throw an error.
            if (matches.size() == 1)
            {
                cerr << "No test sections in regex!" << endl;
                return false;
            }

            // Check each regex capture group in turn
            for (int i = 1; i < matches.size(); ++i)
            {
                // Extract the captured string
                string match(matches[i].first, matches[i].second);
                double val;
                try
                {
                    val = boost::lexical_cast<double>(match);

                    // Add value to the list of found execution times.
                    times.push_back(val);
                    matched = true;
                }
                catch (boost::bad_lexical_cast &e)
                {
                    cerr << "Could not convert match " << match << " to double"
                         << endl;
                    success = false;
                    continue;
                }
            }
        }
    }

    if (!matched)
    {
        cerr << "No execution times were found in the test output. Consider "
                "providing a custom regex for this test metric."
             << endl;
        success = false;
    }

    // Find the minimum of the execution times.
    double minTime = times[0];
    for (unsigned int i = 0; i < times.size(); ++i)
    {
        minTime = (times[i] < minTime) ? times[i] : minTime;
    }

    // If a skip flag is set (due to no matching hostname) then return
    // successful test and output the average time.
    if (m_match.m_skip)
    {
        cerr << "Minimum execution time for host "
             << boost::asio::ip::host_name() << ": " << minTime << endl;
        return success;
    }

    ASSERTL0(m_match.m_value != "nil",
             "No test conditions defined for execution time.");

    // Check that the average time is within the tolerance.
    if (fabs(minTime - boost::lexical_cast<double>(m_match.m_value)) >
            m_match.m_tolerance &&
        !m_match.m_skip)
    {
        cerr << endl;
        cerr << "Failed tolerance match." << endl;
        cerr << "  Expected min execution time: " << m_match.m_value << " +/- "
             << m_match.m_tolerance << endl;
        cerr << "  Actual min: " << minTime << endl;

        for (unsigned int i = 0; i < times.size(); ++i)
        {
            cerr << "       Run " << i << ": " << times[i] << endl;
        }

        success = false;
    }

    return success;
}

/**
 * @brief Generate an accepted execution time value using a test output or error
 * stream.
 *
 * This function generate an execution time match by finding the first execution
 * time in the output that matches \p m_regex (provided in the test file).
 *
 * @param pStdout Reference to test output stream.
 * @param pStderr Reference to test error stream.
 */
void MetricExecutionTime::v_Generate(istream &pStdout, istream &pStderr)
{
    boost::ignore_unused(pStderr);

    // Select istream to use
    istream &is = m_useStderr ? pStderr : pStdout;

    std::smatch matches;

    string line;
    // Vector of execution times found in the output
    vector<double> sampleTimes;
    bool matched = false;

    // Process output file line by line searching for regex matches
    while (getline(is, line))
    {
        // Test to see if we have a match on this line
        if (std::regex_match(line, matches, m_regex))
        {
            // If no fields in regex then throw an error
            ASSERTL0(matches.size() != 1, "No test sections in regex!");

            // Check each regex capture group in turn
            for (int i = 1; i < matches.size(); ++i)
            {
                // Extract the captured string
                string match(matches[i].first, matches[i].second);
                double val;
                try
                {
                    val = boost::lexical_cast<double>(match);

                    sampleTimes.push_back(val);
                    matched = true;
                }
                catch (boost::bad_lexical_cast &e)
                {
                    cerr << "Could not convert match " << match << " to double"
                         << endl;
                    continue;
                }
            }
        }
    }

    // Stores the accepted average execution time and a default tolerance.
    MetricExecutionTimeFieldValue okValue;

    if (matched)
    {
        // Average the found execution times and set it as the accepted value.
        double minTime = sampleTimes[0];
        for (unsigned int i = 0; i < sampleTimes.size(); ++i)
        {
            minTime = (sampleTimes[i] < minTime ? sampleTimes[i] : minTime);
        }

        okValue.m_value = to_string(minTime);
        m_match         = okValue;
    }
    else
    {
        cerr << "No execution times were found in the test output. Value set "
                "to 'nil'."
             << endl;
    }

    // If we are not a derived class then create a new structure.
    if (m_type == "EXECUTIONTIME")
    {
        // Remove values if they already exist.
        while (m_metric->FirstChildElement("value"))
        {
            ASSERTL0(
                m_metric->RemoveChild(m_metric->FirstChildElement("value")),
                "Couldn't remove value from metric!");
        }

        TiXmlElement *val = new TiXmlElement("value");

        val->SetAttribute("tolerance",
                          boost::lexical_cast<string>(m_match.m_tolerance));
        val->SetAttribute("hostname", boost::asio::ip::host_name());
        val->LinkEndChild(new TiXmlText(m_match.m_value));
        m_metric->LinkEndChild(val);
    }
}
} // namespace Nektar
