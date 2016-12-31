/*=====================================================================*
 *                   Copyright (C) 2012 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

#ifndef __CAST_H_

#ifdef __cplusplus
#define cast_uint32_t static_cast<uint32_t>
#else
#define cast_uint32_t (uint32_t)
#endif

#endif // __CAST_H_
/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

#ifndef __SSE_H_
#define __SSE_H_

#ifdef __SSE2__

#include <emmintrin.h>

#ifdef __cplusplus
namespace {
#endif // __cplusplus

typedef __m128 v4sf;
typedef __m128i v4si;

#define v4si_to_v4sf _mm_cvtepi32_ps
#define v4sf_to_v4si _mm_cvttps_epi32

#define v4sfl(x) ((const v4sf) { (x), (x), (x), (x) })
#define v2dil(x) ((const v4si) { (x), (x) })
#define v4sil(x) v2dil((((unsigned long long) (x)) << 32) | (x))

typedef union { v4sf f; float array[4]; } v4sfindexer;
#define v4sf_index(_findx, _findi)      \
  ({                                    \
     v4sfindexer _findvx = { _findx } ; \
     _findvx.array[_findi];             \
  })
typedef union { v4si i; int array[4]; } v4siindexer;
#define v4si_index(_iindx, _iindi)      \
  ({                                    \
     v4siindexer _iindvx = { _iindx } ; \
     _iindvx.array[_iindi];             \
  })

typedef union { v4sf f; v4si i; } v4sfv4sipun;
#define v4sf_fabs(x)                    \
  ({                                    \
     v4sfv4sipun vx;                    \
     vx.f = x;                          \
     vx.i &= v4sil (0x7FFFFFFF);        \
     vx.f;                              \
  })

#ifdef __cplusplus
} // end namespace
#endif // __cplusplus

#endif // __SSE2__

#endif // __SSE_H_
/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

#ifndef __FAST_EXP_H_
#define __FAST_EXP_H_

#include <stdint.h>

// Underflow of exponential is common practice in numerical routines,
// so handle it here.

static inline float
fastpow2 (float p)
{
  float offset = (p < 0) ? 1.0f : 0.0f;
  float clipp = (p < -126) ? -126.0f : p;
  int w = clipp;
  float z = clipp - w + offset;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) ) };

  return v.f;
}

static inline float
fastexp (float p)
{
  return fastpow2 (1.442695040f * p);
}

static inline float
fasterpow2 (float p)
{
  float clipp = (p < -126) ? -126.0f : p;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 126.94269504f) ) };
  return v.f;
}

static inline float
fasterexp (float p)
{
  return fasterpow2 (1.442695040f * p);
}

#ifdef __SSE2__

static inline v4sf
vfastpow2 (const v4sf p)
{
  v4sf ltzero = _mm_cmplt_ps (p, v4sfl (0.0f));
  v4sf offset = _mm_and_ps (ltzero, v4sfl (1.0f));
  v4sf lt126 = _mm_cmplt_ps (p, v4sfl (-126.0f));
  v4sf clipp = _mm_or_ps (_mm_andnot_ps (lt126, p), _mm_and_ps (lt126, v4sfl (-126.0f)));
  v4si w = v4sf_to_v4si (clipp);
  v4sf z = clipp - v4si_to_v4sf (w) + offset;

  const v4sf c_121_2740838 = v4sfl (121.2740575f);
  const v4sf c_27_7280233 = v4sfl (27.7280233f);
  const v4sf c_4_84252568 = v4sfl (4.84252568f);
  const v4sf c_1_49012907 = v4sfl (1.49012907f);
  union { v4si i; v4sf f; } v = {
    v4sf_to_v4si (
      v4sfl (1 << 23) * 
      (clipp + c_121_2740838 + c_27_7280233 / (c_4_84252568 - z) - c_1_49012907 * z)
    )
  };

  return v.f;
}

static inline v4sf
vfastexp (const v4sf p)
{
  const v4sf c_invlog_2 = v4sfl (1.442695040f);

  return vfastpow2 (c_invlog_2 * p);
}

static inline v4sf
vfasterpow2 (const v4sf p)
{
  const v4sf c_126_94269504 = v4sfl (126.94269504f);
  v4sf lt126 = _mm_cmplt_ps (p, v4sfl (-126.0f));
  v4sf clipp = _mm_or_ps (_mm_andnot_ps (lt126, p), _mm_and_ps (lt126, v4sfl (-126.0f)));
  union { v4si i; v4sf f; } v = { v4sf_to_v4si (v4sfl (1 << 23) * (clipp + c_126_94269504)) };
  return v.f;
}

static inline v4sf
vfasterexp (const v4sf p)
{
  const v4sf c_invlog_2 = v4sfl (1.442695040f);

  return vfasterpow2 (c_invlog_2 * p);
}

#endif //__SSE2__

#endif // __FAST_EXP_H_
/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

#ifndef __FAST_LOG_H_
#define __FAST_LOG_H_

#include <stdint.h>

static inline float 
fastlog2 (float x)
{
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
           - 1.498030302f * mx.f 
           - 1.72587999f / (0.3520887068f + mx.f);
}

static inline float
fastlog (float x)
{
  return 0.69314718f * fastlog2 (x);
}

static inline float 
fasterlog2 (float x)
{
  union { float f; uint32_t i; } vx = { x };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;
  return y - 126.94269504f;
}

static inline float
fasterlog (float x)
{
//  return 0.69314718f * fasterlog2 (x);

  union { float f; uint32_t i; } vx = { x };
  float y = vx.i;
  y *= 8.2629582881927490e-8f;
  return y - 87.989971088f;
}

#ifdef __SSE2__

static inline v4sf
vfastlog2 (v4sf x)
{
  union { v4sf f; v4si i; } vx = { x };
  union { v4si i; v4sf f; } mx; mx.i = (vx.i & v4sil (0x007FFFFF)) | v4sil (0x3f000000);
  v4sf y = v4si_to_v4sf (vx.i);
  y *= v4sfl (1.1920928955078125e-7f);

  const v4sf c_124_22551499 = v4sfl (124.22551499f);
  const v4sf c_1_498030302 = v4sfl (1.498030302f);
  const v4sf c_1_725877999 = v4sfl (1.72587999f);
  const v4sf c_0_3520087068 = v4sfl (0.3520887068f);

  return y - c_124_22551499
           - c_1_498030302 * mx.f 
           - c_1_725877999 / (c_0_3520087068 + mx.f);
}

static inline v4sf
vfastlog (v4sf x)
{
  const v4sf c_0_69314718 = v4sfl (0.69314718f);

  return c_0_69314718 * vfastlog2 (x);
}

static inline v4sf 
vfasterlog2 (v4sf x)
{
  union { v4sf f; v4si i; } vx = { x };
  v4sf y = v4si_to_v4sf (vx.i);
  y *= v4sfl (1.1920928955078125e-7f);

  const v4sf c_126_94269504 = v4sfl (126.94269504f);

  return y - c_126_94269504;
}

static inline v4sf
vfasterlog (v4sf x)
{
//  const v4sf c_0_69314718 = v4sfl (0.69314718f);
//
//  return c_0_69314718 * vfasterlog2 (x);

  union { v4sf f; v4si i; } vx = { x };
  v4sf y = v4si_to_v4sf (vx.i);
  y *= v4sfl (8.2629582881927490e-8f);

  const v4sf c_87_989971088 = v4sfl (87.989971088f);

  return y - c_87_989971088;
}
#endif // __SSE2__

#endif // __FAST_LOG_H_
