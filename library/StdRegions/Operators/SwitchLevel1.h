///////////////////////////////////////////////////////////////////////////////
//
// File: StdSwitchLevel1.h
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
// Description: One level switch statement with definable bounds for 3D elmts
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#include "Operators/SwitchLimits.h"

/* The following code provide BOOST_PP macrros to sets up a switch
   statement that goes from  SMIN to SMAX */

/** this macro tests the values of the tuple 'state' to see if the
    first element, given by BOOST_PP_TUPLE_ELEM(0, state), is not
    equal to the second element plus one, given by
    BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(1, state)) and returns 1 if not
    equal otherwise returns zero. */
#define STDLEV1TEST(r, state)                                                  \
    BOOST_PP_NOT_EQUAL(BOOST_PP_TUPLE_ELEM(0, state),                          \
                       BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(1, state)))

/** This macro returns an updated tuple where the first element, given
    by BOOST_PP_TUPLE_ELEM(0, state), is incremented by one and the
    second element, given by BOOST_PP_TUPLE_ELEM(1, state) remains the
    same. */
#define STDLEV1UPDATE(r, state)                                                \
    (BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(0, state)), BOOST_PP_TUPLE_ELEM(1, state))

#undef NQ1
#define NQ1(i) BOOST_PP_TUPLE_ELEM(0, i)
#undef NQ1_M1
#define NQ1_M1(i) BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(0, i))
#undef NQ1_P1
#define NQ1_P1(i) BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(0, i))
