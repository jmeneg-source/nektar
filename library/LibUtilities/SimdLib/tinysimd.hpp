///////////////////////////////////////////////////////////////////////////////
//
// File: tinysimd.hpp
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
// Description: Light wrapper for automatic selection of available SIMD type.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_TINYSIMD_H
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_TINYSIMD_H

#include "avx2.hpp"
#include "avx512.hpp"
#include "scalar.hpp"
#include "sse2.hpp"
#include "sve.hpp"
#include <array>
#include <cstdint>
#include <vector>

namespace tinysimd
{

template <typename...> struct first_not_void_of
{
    using type = void;
};

template <typename... Rest> struct first_not_void_of<void, Rest...>
{
    using type = typename first_not_void_of<Rest...>::type;
};

template <typename T, typename... Rest> struct first_not_void_of<T, Rest...>
{
    using type = T;
};

template <unsigned int Width> struct width_tag
{
    static constexpr unsigned int width = Width;
};

template <typename T, typename = void>
struct width_or_zero : std::integral_constant<unsigned int, 0>
{
};

template <typename T>
struct width_or_zero<T, std::void_t<decltype(T::width)>>
    : std::integral_constant<unsigned int, T::width>
{
};

namespace abi
{

// pick the most specialiazed available match out of available ABIs.
template <typename T, int width> struct default_abi
{
    using type = typename first_not_void_of<
        typename sve<T, width>::type, typename avx512<T, width>::type,
        typename avx2<T, width>::type, typename sse2<T, width>::type,
        typename simd64<T>::type, typename scalar<T>::type>::type;

    static_assert(!std::is_void_v<type>, "unsupported SIMD type");
};

} // namespace abi

// Global long-SIMD scaling knob. PACKSIZE is derived from this multiplier,
// the selected architecture, and the scalar type.
static constexpr unsigned int PACKMULTIPLIER = 1;

namespace abi
{

template <typename T> struct default_longsimd_width
{
private:
    using native_type = typename default_abi<T, 0>::type;

public:
    static constexpr unsigned int value = native_type::width * PACKMULTIPLIER;

    static_assert(value > 0, "unsupported SIMD type");
};

} // namespace abi

// Global long-SIMD width used by library code for NekDouble packs.
static constexpr unsigned int PACKSIZE =
    abi::default_longsimd_width<double>::value;

#if defined(__clang__)
#define TINYSIMD_PRAGMA_UNROLL _Pragma("clang loop unroll(full)")
#elif defined(__GNUC__)
#define TINYSIMD_PRAGMA_UNROLL _Pragma("GCC unroll 16")
#else
#define TINYSIMD_PRAGMA_UNROLL
#endif

// light wrapper for default types
template <typename ScalarType, int width = 0,
          template <typename, int> class abi = abi::default_abi>
using simd = typename abi<ScalarType, width>::type;

namespace details
{

template <typename T, typename = void> struct scalar_index_type
{
    using type = std::conditional_t<(sizeof(typename T::scalarType) <= 4),
                                    std::uint32_t, std::uint64_t>;
};

template <typename T>
struct scalar_index_type<T, std::enable_if_t<has_scalarIndexType<T>::value>>
{
    using type = typename T::scalarIndexType;
};

template <typename T>
using scalar_index_type_t = typename scalar_index_type<T>::type;

template <typename SimdType, unsigned int Width> struct long_simd
{
    static_assert(is_vector_floating_point_v<SimdType>,
                  "longsimd requires a floating-point SIMD chunk type");

    static constexpr unsigned int width      = Width;
    static constexpr unsigned int alignment  = SimdType::alignment;
    static constexpr unsigned int chunkWidth = SimdType::width;
    static constexpr unsigned int numChunks  = width / chunkWidth;

    static_assert(width >= chunkWidth,
                  "longsimd width must be at least one native SIMD chunk");
    static_assert(width % chunkWidth == 0,
                  "longsimd width must be a multiple of the chunk width");

    using scalarType      = typename SimdType::scalarType;
    using scalarIndexType = scalar_index_type_t<SimdType>;
    using chunkType       = SimdType;
    using chunkIndexType  = simd<scalarIndexType, chunkWidth>;
    using vectorType      = std::array<chunkType, numChunks>;
    using scalarArray     = scalarType[width];

    static_assert(sizeof(chunkType) == sizeof(scalarType) * chunkWidth,
                  "longsimd requires tightly packed SIMD chunk storage");

    alignas(alignment) chunkType _data[numChunks];

    inline long_simd()                     = default;
    inline long_simd(const long_simd &rhs) = default;

    inline void broadcast(const chunkType &rhs)
    {
        _data[0] = rhs;

        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 1; i < numChunks; ++i)
        {
            _data[i] = _data[0];
        }
    }

    inline long_simd(const vectorType &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] = rhs[i];
        }
    }

    inline long_simd(const chunkType *rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] = rhs[i];
        }
    }

    inline long_simd(const scalarType rhs)
    {
        broadcast(chunkType(rhs));
    }

    inline long_simd &operator=(const long_simd &) = default;

    inline void store(scalarType *p) const
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].store(p + i * chunkWidth);
        }
    }

    template <class flag,
              typename std::enable_if<is_load_tag_v<flag>, bool>::type = 0>
    inline void store(scalarType *p, flag f) const
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].store(p + i * chunkWidth, f);
        }
    }

    inline void load(const scalarType *p)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].load(p + i * chunkWidth);
        }
    }

    template <class flag,
              typename std::enable_if<is_load_tag_v<flag>, bool>::type = 0>
    inline void load(const scalarType *p, flag f)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].load(p + i * chunkWidth, f);
        }
    }

    inline void broadcast(const scalarType rhs)
    {
        broadcast(chunkType(rhs));
    }

    template <typename IndexType>
    inline void gather(const scalarType *p, const IndexType *indices)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].gather(p, indices[i]);
        }
    }

    template <typename IndexType>
    inline void scatter(scalarType *out, const IndexType *indices) const
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].scatter(out, indices[i]);
        }
    }

    inline void fma(const long_simd &a, const long_simd &b)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].fma(a._data[i], b._data[i]);
        }
    }

    inline void fma(const long_simd &a, const chunkType &b)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].fma(a._data[i], b);
        }
    }

    inline void fma(const chunkType &a, const long_simd &b)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i].fma(a, b._data[i]);
        }
    }

    inline void fma(const chunkType &a, const chunkType &b)
    {
        const chunkType ab = a * b;
        *this += ab;
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void fma(const long_simd &a, U b)
    {
        fma(a, chunkType(static_cast<scalarType>(b)));
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void fma(U a, const long_simd &b)
    {
        fma(chunkType(static_cast<scalarType>(a)), b);
    }

    template <
        typename U, typename V,
        typename std::enable_if<
            std::is_arithmetic_v<U> && std::is_arithmetic_v<V>, bool>::type = 0>
    inline void fma(U a, V b)
    {
        fma(chunkType(static_cast<scalarType>(a)),
            chunkType(static_cast<scalarType>(b)));
    }

    inline scalarType operator[](size_t i) const
    {
        alignas(alignment) scalarArray tmp;
        store(tmp, is_aligned);
        return tmp[i];
    }

    inline scalarType &operator[](size_t i)
    {
        scalarType *tmp = reinterpret_cast<scalarType *>(_data);
        return tmp[i];
    }

    inline void operator+=(const long_simd &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] += rhs._data[i];
        }
    }

    inline void operator+=(const chunkType &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] += rhs;
        }
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void operator+=(U rhs)
    {
        *this += chunkType(static_cast<scalarType>(rhs));
    }

    inline void operator-=(const long_simd &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] -= rhs._data[i];
        }
    }

    inline void operator-=(const chunkType &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] -= rhs;
        }
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void operator-=(U rhs)
    {
        *this -= chunkType(static_cast<scalarType>(rhs));
    }

    inline void operator*=(const long_simd &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] *= rhs._data[i];
        }
    }

    inline void operator*=(const chunkType &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] *= rhs;
        }
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void operator*=(U rhs)
    {
        *this *= chunkType(static_cast<scalarType>(rhs));
    }

    inline void operator/=(const long_simd &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] /= rhs._data[i];
        }
    }

    inline void operator/=(const chunkType &rhs)
    {
        TINYSIMD_PRAGMA_UNROLL
        for (size_t i = 0; i < numChunks; ++i)
        {
            _data[i] /= rhs;
        }
    }

    template <typename U,
              typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
    inline void operator/=(U rhs)
    {
        *this /= chunkType(static_cast<scalarType>(rhs));
    }
};

template <typename ScalarType, int width,
          template <typename, int> class abi_policy>
struct longsimd_selector
{
private:
    using native_type = simd<ScalarType, 0, abi_policy>;

    static constexpr int default_width =
        static_cast<int>(abi::default_longsimd_width<ScalarType>::value);

public:
    static constexpr int requested_width = width == 0 ? default_width : width;

    static_assert(requested_width > 0, "longsimd width must be positive");
    static_assert(
        requested_width >= static_cast<int>(native_type::width),
        "longsimd width must be at least the native SIMD width; use simd<> "
        "for narrower types");
    static_assert(requested_width % static_cast<int>(native_type::width) == 0,
                  "longsimd width must be a multiple of the native SIMD width");

    using type = std::conditional_t<
        requested_width == static_cast<int>(native_type::width), native_type,
        long_simd<native_type, static_cast<unsigned int>(requested_width)>>;
};

} // namespace details

template <typename ScalarType, int width = 0,
          template <typename, int> class abi = abi::default_abi>
using longsimd =
    typename details::longsimd_selector<ScalarType, width, abi>::type;

#if defined(__AVX2__) && defined(NEKTAR_ENABLE_SIMD_AVX2)
using avx2Double8 = details::long_simd<avx2Double4, 8>;
#endif

#if defined(__AVX512F__) && defined(NEKTAR_ENABLE_SIMD_AVX512)
using avx512Double16 = details::long_simd<avx512Double8, 16>;
#endif

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator+(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator+(
    details::long_simd<SimdType, Width> lhs,
    const typename details::long_simd<SimdType, Width>::chunkType &rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator+(
    const typename details::long_simd<SimdType, Width>::chunkType &lhs,
    details::long_simd<SimdType, Width> rhs)
{
    rhs += lhs;
    return rhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator+(
    details::long_simd<SimdType, Width> lhs, U rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator+(
    U lhs, details::long_simd<SimdType, Width> rhs)
{
    using vec_t = details::long_simd<SimdType, Width>;
    return typename vec_t::chunkType(
               static_cast<typename vec_t::scalarType>(lhs)) +
           rhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator-(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator-(
    details::long_simd<SimdType, Width> lhs,
    const typename details::long_simd<SimdType, Width>::chunkType &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator-(
    const typename details::long_simd<SimdType, Width>::chunkType &lhs,
    details::long_simd<SimdType, Width> rhs)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < rhs.numChunks; ++i)
    {
        rhs._data[i] = lhs - rhs._data[i];
    }
    return rhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator-(
    details::long_simd<SimdType, Width> lhs, U rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator-(
    U lhs, details::long_simd<SimdType, Width> rhs)
{
    using vec_t = details::long_simd<SimdType, Width>;
    return typename vec_t::chunkType(
               static_cast<typename vec_t::scalarType>(lhs)) -
           rhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator-(
    details::long_simd<SimdType, Width> in)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < in.numChunks; ++i)
    {
        in._data[i] = -in._data[i];
    }
    return in;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator*(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator*(
    details::long_simd<SimdType, Width> lhs,
    const typename details::long_simd<SimdType, Width>::chunkType &rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator*(
    const typename details::long_simd<SimdType, Width>::chunkType &lhs,
    details::long_simd<SimdType, Width> rhs)
{
    rhs *= lhs;
    return rhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator*(
    details::long_simd<SimdType, Width> lhs, U rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator*(
    U lhs, details::long_simd<SimdType, Width> rhs)
{
    using vec_t = details::long_simd<SimdType, Width>;
    return typename vec_t::chunkType(
               static_cast<typename vec_t::scalarType>(lhs)) *
           rhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator/(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    lhs /= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator/(
    details::long_simd<SimdType, Width> lhs,
    const typename details::long_simd<SimdType, Width>::chunkType &rhs)
{
    lhs /= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> operator/(
    const typename details::long_simd<SimdType, Width>::chunkType &lhs,
    details::long_simd<SimdType, Width> rhs)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < rhs.numChunks; ++i)
    {
        rhs._data[i] = lhs / rhs._data[i];
    }
    return rhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator/(
    details::long_simd<SimdType, Width> lhs, U rhs)
{
    lhs /= rhs;
    return lhs;
}

template <typename SimdType, unsigned int Width, typename U,
          typename std::enable_if<std::is_arithmetic_v<U>, bool>::type = 0>
inline details::long_simd<SimdType, Width> operator/(
    U lhs, details::long_simd<SimdType, Width> rhs)
{
    using vec_t = details::long_simd<SimdType, Width>;
    return typename vec_t::chunkType(
               static_cast<typename vec_t::scalarType>(lhs)) /
           rhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> sqrt(
    details::long_simd<SimdType, Width> in)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < in.numChunks; ++i)
    {
        in._data[i] = sqrt(in._data[i]);
    }
    return in;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> abs(
    details::long_simd<SimdType, Width> in)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < in.numChunks; ++i)
    {
        in._data[i] = abs(in._data[i]);
    }
    return in;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> min(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < lhs.numChunks; ++i)
    {
        lhs._data[i] = min(lhs._data[i], rhs._data[i]);
    }
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> max(
    details::long_simd<SimdType, Width> lhs,
    const details::long_simd<SimdType, Width> &rhs)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < lhs.numChunks; ++i)
    {
        lhs._data[i] = max(lhs._data[i], rhs._data[i]);
    }
    return lhs;
}

template <typename SimdType, unsigned int Width>
inline details::long_simd<SimdType, Width> log(
    details::long_simd<SimdType, Width> in)
{
    TINYSIMD_PRAGMA_UNROLL
    for (size_t i = 0; i < in.numChunks; ++i)
    {
        in._data[i] = log(in._data[i]);
    }
    return in;
}

template <typename SimdType, unsigned int Width>
inline void load_unalign_interleave(
    const typename details::long_simd<SimdType, Width>::scalarType *in,
    const std::uint32_t dataLen,
    std::vector<details::long_simd<SimdType, Width>,
                allocator<details::long_simd<SimdType, Width>>> &out)
{
    using vec_t = details::long_simd<SimdType, Width>;

    alignas(vec_t::alignment) typename vec_t::scalarArray tmp;
    for (size_t i = 0; i < dataLen; ++i)
    {
        for (size_t j = 0; j < vec_t::width; ++j)
        {
            tmp[j] = in[i + j * dataLen];
        }
        out[i].load(tmp);
    }
}

template <typename SimdType, unsigned int Width>
inline void load_unalign_interleave_skipPads(
    const typename details::long_simd<SimdType, Width>::scalarType *in,
    const std::uint32_t dataLen, const std::uint32_t skipPads,
    std::vector<details::long_simd<SimdType, Width>,
                allocator<details::long_simd<SimdType, Width>>> &out)
{
    using vec_t = details::long_simd<SimdType, Width>;

    alignas(vec_t::alignment) typename vec_t::scalarArray tmp;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp1;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp2;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp3;

    size_t nBlocks = dataLen / 4;
    const vec_t zero{static_cast<typename vec_t::scalarType>(0)};

    for (size_t i = 0; i < nBlocks; ++i)
    {
        zero.store(tmp);
        zero.store(tmp1);
        zero.store(tmp2);
        zero.store(tmp3);

        for (size_t j = 0; j < vec_t::width - skipPads; ++j)
        {
            tmp[j]  = in[j * dataLen + 4 * i];
            tmp1[j] = in[j * dataLen + 4 * i + 1];
            tmp2[j] = in[j * dataLen + 4 * i + 2];
            tmp3[j] = in[j * dataLen + 4 * i + 3];
        }

        out[4 * i].load(tmp);
        out[4 * i + 1].load(tmp1);
        out[4 * i + 2].load(tmp2);
        out[4 * i + 3].load(tmp3);
    }

    for (size_t i = nBlocks * 4; i < dataLen; ++i)
    {
        zero.store(tmp);
        for (size_t j = 0; j < vec_t::width - skipPads; ++j)
        {
            tmp[j] = in[i + j * dataLen];
        }
        out[i].load(tmp);
    }
}

template <typename SimdType, unsigned int Width>
inline void load_interleave(
    const typename details::long_simd<SimdType, Width>::scalarType *in,
    std::uint32_t dataLen,
    std::vector<details::long_simd<SimdType, Width>,
                allocator<details::long_simd<SimdType, Width>>> &out)
{
    using vec_t   = details::long_simd<SimdType, Width>;
    using index_t = typename vec_t::chunkIndexType;

    alignas(index_t::alignment) typename index_t::scalarArray tmp;
    index_t index[vec_t::numChunks];

    for (size_t chunk = 0; chunk < vec_t::numChunks; ++chunk)
    {
        for (size_t lane = 0; lane < index_t::width; ++lane)
        {
            tmp[lane] = static_cast<typename vec_t::scalarIndexType>(
                (chunk * index_t::width + lane) * dataLen);
        }
        index[chunk].load(tmp);
    }

    for (size_t i = 0; i < dataLen; ++i)
    {
        out[i].gather(in, index);
        for (auto &idx : index)
        {
            idx = idx + index_t(1);
        }
    }
}

template <typename SimdType, unsigned int Width>
inline void deinterleave_unalign_store(
    const std::vector<details::long_simd<SimdType, Width>,
                      allocator<details::long_simd<SimdType, Width>>> &in,
    const std::uint32_t dataLen,
    typename details::long_simd<SimdType, Width>::scalarType *out)
{
    using vec_t = details::long_simd<SimdType, Width>;

    alignas(vec_t::alignment) typename vec_t::scalarArray tmp;
    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].store(tmp);
        for (size_t j = 0; j < vec_t::width; ++j)
        {
            out[i + j * dataLen] = tmp[j];
        }
    }
}

template <typename SimdType, unsigned int Width>
inline void deinterleave_unalign_store_skipPads(
    const std::vector<details::long_simd<SimdType, Width>,
                      allocator<details::long_simd<SimdType, Width>>> &in,
    const std::uint32_t dataLen, const std::uint32_t skipPads,
    typename details::long_simd<SimdType, Width>::scalarType *out)
{
    using vec_t = details::long_simd<SimdType, Width>;

    alignas(vec_t::alignment) typename vec_t::scalarArray tmp;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp1;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp2;
    alignas(vec_t::alignment) typename vec_t::scalarArray tmp3;

    size_t nBlocks = dataLen / 4;

    for (size_t i = 0; i < nBlocks; ++i)
    {
        in[4 * i].store(tmp);
        in[4 * i + 1].store(tmp1);
        in[4 * i + 2].store(tmp2);
        in[4 * i + 3].store(tmp3);

        for (size_t j = 0; j < vec_t::width - skipPads; ++j)
        {
            out[j * dataLen + 4 * i]     = tmp[j];
            out[j * dataLen + 4 * i + 1] = tmp1[j];
            out[j * dataLen + 4 * i + 2] = tmp2[j];
            out[j * dataLen + 4 * i + 3] = tmp3[j];
        }
    }

    for (size_t i = nBlocks * 4; i < dataLen; ++i)
    {
        in[i].store(tmp);
        for (size_t j = 0; j < vec_t::width - skipPads; ++j)
        {
            out[j * dataLen + i] = tmp[j];
        }
    }
}

template <typename SimdType, unsigned int Width>
inline void deinterleave_store(
    const std::vector<details::long_simd<SimdType, Width>,
                      allocator<details::long_simd<SimdType, Width>>> &in,
    std::uint32_t dataLen,
    typename details::long_simd<SimdType, Width>::scalarType *out)
{
    using vec_t   = details::long_simd<SimdType, Width>;
    using index_t = typename vec_t::chunkIndexType;

    alignas(index_t::alignment) typename index_t::scalarArray tmp;
    index_t index[vec_t::numChunks];

    for (size_t chunk = 0; chunk < vec_t::numChunks; ++chunk)
    {
        for (size_t lane = 0; lane < index_t::width; ++lane)
        {
            tmp[lane] = static_cast<typename vec_t::scalarIndexType>(
                (chunk * index_t::width + lane) * dataLen);
        }
        index[chunk].load(tmp);
    }

    for (size_t i = 0; i < dataLen; ++i)
    {
        in[i].scatter(out, index);
        for (auto &idx : index)
        {
            idx = idx + index_t(1);
        }
    }
}

#undef TINYSIMD_PRAGMA_UNROLL

} // namespace tinysimd
#endif
