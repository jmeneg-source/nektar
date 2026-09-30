///////////////////////////////////////////////////////////////////////////////
//
// File: SwitchPoints.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

// This code is the mombo switch statement that is used in mutliple
// operators. It uses preprocessor directives based on the shape
// type and dimension to limit the code inclusion.

// This header is included in the OperatorType.h file. Ideally the
// operator() method would be a template in the Operator base class
// but becasue the operator{1,23}D is both a function and template that
// need to be in the inherited class it is not possible.

#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

#include "SwitchLimits.h"
/* Macro compares the values of the tuple 'state' to see if the first
   element, given by BOOST_PP_TUPLE_ELEM(2, 0, state), is not equal to
   the second element plus one, given by
   BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)) and returns 1 if not
   equal otherwise zero */
#define TEST(r, state)                                                         \
    BOOST_PP_NOT_EQUAL(BOOST_PP_TUPLE_ELEM(2, 0, state),                       \
                       BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state)))

/* Macro returns a tuple where the first element, given by
   BOOST_PP_TUPLE_ELEM(2, 0, state), is incremented by one and the
   second element, given by BOOST_PP_TUPLE_ELEM(2, 1, state) remains
   the same*/
#define UPDATE(r, state)                                                       \
    (BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),                           \
     BOOST_PP_TUPLE_ELEM(2, 1, state))

#if defined(SHAPE_DIMENSION_1D)

/* Define the switch statement for 1D shape operators and min, max range */
#define OPERATOR1D(r, state)                                                   \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator1D<BOOST_PP_TUPLE_ELEM(2, 0, state)>(input, output);           \
        break;

/* start of switch include for operator1D */
{
#if defined(SHAPE_TYPE_SEG)
    switch (nq0)
    {
        /*
          expand OPERATOR1D case from 2 to 10, i.e.
          ..
          case 4:
             operator1D<4>(input,output);
             break;
          ...
        */
        BOOST_PP_FOR((MIN1D, MAX1D), TEST, UPDATE, OPERATOR1D)
        default:
            operator1D(input, output);
            break;
    }
#endif
}

#elif defined(SHAPE_DIMENSION_2D)

/* Define the switch statement for 2D shape operators and min, max range */
#define OPERATOR2D_TRI(r, state)                                               \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator2D<BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(2, 0, state))>(input,      \
                                                                   output);    \
        break;

#define OPERATOR2D_QUAD(r, state)                                              \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator2D<BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_TUPLE_ELEM(2, 0, state)>(input, output);           \
        break;

/* start of switch include for operator2D */
{
#if defined(SHAPE_TYPE_TRI)

    if (nq0 == nq1 + 1)
    {
        switch (m_basis[0]->GetNumPoints())
        {
            /*
              expand OPERATOR2D case from 2 to 10, i.e.
              ..
              case 4:
                 operator2D<4,3>(input,output);
                 break;
              ...
            */
            BOOST_PP_FOR((MIN2D, MAX2D), TEST, UPDATE, OPERATOR2D_TRI);

            default:
                operator2D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_QUAD)

    if (nq0 == nq1)
    {
        switch (m_basis[0]->GetNumPoints())
        {
            /*
              expand OPERATOR2D case from 2 to 10, i.e.
              ..
              case 4:
                 operator2D<4,4>(input,output);
                 break;
              ...
            */
            BOOST_PP_FOR((MIN2D, MAX2D), TEST, UPDATE, OPERATOR2D_QUAD)
            default:
                operator2D(input, output);
                break;
        }
    }

#endif // SHAPE_TYPE
    else
    {
        operator2D(input, output);
    }
}

#elif defined(SHAPE_DIMENSION_3D)

/* Define the switch statement for 3D shape operators and min, max range */
#define OPERATOR3D_HEX(r, state)                                               \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator3D<BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_TUPLE_ELEM(2, 0, state)>(input, output);           \
        break;

#define OPERATOR3D_TET(r, state)                                               \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator3D<BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(2, 0, state)),             \
                   BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(2, 0, state))>(input,      \
                                                                   output);    \
        break;

#define OPERATOR3D_PYR_PRISM(r, state)                                         \
    case BOOST_PP_TUPLE_ELEM(2, 0, state):                                     \
        operator3D<BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_TUPLE_ELEM(2, 0, state),                           \
                   BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(2, 0, state))>(input,      \
                                                                   output);    \
        break;

/* start of include switch details for operator3D */
{
#if defined(SHAPE_TYPE_HEX)

    if (nq0 == nq1 && nq0 == nq2)
    {
        switch (m_basis[0]->GetNumPoints())
        {
            /*
              expand OPERATOR3D case from 2 to 16, i.e.
              ..
              case 4:
                 operator3D<4,4,4>(input,output);
                 break;
              ...
            */
            BOOST_PP_FOR((MIN3D, MAX3D), TEST, UPDATE, OPERATOR3D_HEX)

            default:
                operator3D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_TET)

    if (nq0 == nq1 + 1 && nq0 == nq2 + 1)
    {
        switch (nq0)
        {
            /*
              expand OPERATOR3D case from 2 to 16, i.e.
              ..
              case 4:
                 operator3D<4,3,3>(input,output);
                 break;
              ...
            */
            BOOST_PP_FOR((MIN3D, MAX3D), TEST, UPDATE, OPERATOR3D_TET)
            default:
                operator3D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_PYR) || defined(SHAPE_TYPE_PRISM)

    if (nq0 == nq1 && nq0 == nq2 + 1)
    {
        switch (m_basis[0]->GetNumPoints())
        {
            /*
              expand OPERATOR3D case from 3 to 16, i.e.
              ..
              case 4:
                 operator3D<4,3,3>(input,output);
                 break;
              ...
            */
            BOOST_PP_FOR((MIN3D, MAX3D), TEST, UPDATE, OPERATOR3D_PYR_PRISM)
            default:
                operator3D(input, output);
                break;
        }
    }

#endif // SHAPE_TYPE

    else
    {
        operator3D(input, output);
    }
}

#endif // SHAPE_DIMENSION
