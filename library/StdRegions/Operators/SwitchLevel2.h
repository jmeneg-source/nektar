///////////////////////////////////////////////////////////////////////////////
//
// File: SwitchLevel2.h
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
// Description: Boost preprocessing macro for two level switch
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#include "Operators/SwitchLimits.h"

/* Switch macros for two level switch over modes and then quadrature
   orders wrapped around a Deformed check */
#undef NM
#define NM(i) BOOST_PP_TUPLE_ELEM(0, i)
#undef NM_P1
#define NM_P1(i) BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(0, i))
#undef NQ
#define NQ(i) BOOST_PP_TUPLE_ELEM(1, i)
#undef NQ_M1
#define NQ_M1(i) BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(1, i))

/** Macro tests the values of the tuple 'state' to see if the first
   element, given by BOOST_PP_TUPLE_ELEM(0, state), is not equal to
   the second element plus one, given by
   BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(1, state)) and returns 1 if not
   equal otherwise zero */
#define STDLEV2TEST(r, state)                                                  \
    BOOST_PP_NOT_EQUAL(BOOST_PP_TUPLE_ELEM(0, state),                          \
                       BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, state)))

#define STDLEV2TEST1(r, state)                                                 \
    BOOST_PP_NOT_EQUAL(BOOST_PP_TUPLE_ELEM(1, state),                          \
                       BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, state)))

/** Macro returns an updated tuple where the first element, given by
   BOOST_PP_TUPLE_ELEM(0, state), is incremented by one and the
   second element, given by BOOST_PP_TUPLE_ELEM(1, state) remains
   the same*/
#define STDLEV2UPDATE(r, state)                                                \
    (BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(0, state)),                              \
     BOOST_PP_TUPLE_ELEM(1, state), BOOST_PP_TUPLE_ELEM(2, state))

#define STDLEV2UPDATE1(r, state)                                               \
    (BOOST_PP_TUPLE_ELEM(0, state),                                            \
     BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(1, state)),                              \
     BOOST_PP_TUPLE_ELEM(2, state))
