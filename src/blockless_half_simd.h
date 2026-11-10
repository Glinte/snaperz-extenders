#pragma once
// Finite AVX2 wavefront on up to 64 byte coordinates. Gate k of sweep t
// executes at tick 2*t+k (or 2*t+n-2-k for inverse sweeps). Gates at a
// tick are disjoint; warm-up/drain masks make the returned state synchronous.
namespace snaperz { namespace blockless { namespace half {
inline __m256i left1(__m256i v)
{
  const auto lo = _mm256_permute2x128_si256(v, v, 0x08);
  return _mm256_alignr_epi8(v, lo, 15);
}

template<bool Backward = false>
inline void power_simd(State& a, unsigned count, size_t start = 0)
{
  const size_t n = a.size() - start;
  if (n < 2 || !count) return;
  if (n > 64 || *std::max_element(a.begin()+start, a.end()) > 127)
  {
    for (unsigned t=0; t<count; ++t)
      if constexpr (Backward) inverse(a,start); else sweep(a,start);
    return;
  }
  alignas(32) uint8_t ev[32] = {}, od[32] = {};
  for (size_t i=0; i<n; ++i) (i&1 ? od : ev)[i/2] = a[start+i];
  auto e = _mm256_load_si256(reinterpret_cast<const __m256i*>(ev));
  auto o = _mm256_load_si256(reinterpret_cast<const __m256i*>(od));
  const auto zero = _mm256_setzero_si256(), one = _mm256_set1_epi8(1);
  alignas(32) uint8_t indices[32];
  for (unsigned i=0;i<32;++i) indices[i]=2*i;
  const auto idx = _mm256_load_si256(reinterpret_cast<const __m256i*>(indices));
  const unsigned last = unsigned(n)-2;
  const unsigned end = 2*(count-1)+last;
  for (unsigned tick=0; tick<=end; ++tick)
  {
    const bool odd = (tick ^ (Backward ? last : 0)) & 1;
    auto l = odd ? o : e;
    auto r = odd ? shr1(e) : o;
    auto k = _mm256_add_epi8(idx, _mm256_set1_epi8(odd));
    auto mask = _mm256_cmpgt_epi8(_mm256_set1_epi8(last+1), k);
    if (tick < last || tick > 2*(count-1))
    {
      auto distance = Backward ? _mm256_sub_epi8(_mm256_set1_epi8(last), k) : k;
      if (tick < last)
        mask = _mm256_and_si256(mask, _mm256_cmpgt_epi8(_mm256_set1_epi8(tick+1), distance));
      if (tick > 2*(count-1))
        mask = _mm256_and_si256(mask, _mm256_cmpgt_epi8(distance, _mm256_set1_epi8(tick-2*(count-1)-1)));
    }
    __m256i delta;
    if constexpr (Backward)
    {
      delta = _mm256_sub_epi8(one, _mm256_and_si256(_mm256_cmpeq_epi8(r,zero),l));
      mask = _mm256_andnot_si256(_mm256_cmpeq_epi8(l,zero),mask);
    }
    else
      delta = _mm256_add_epi8(_mm256_and_si256(_mm256_cmpeq_epi8(l,one),r), _mm256_cmpgt_epi8(l,one));
    delta = _mm256_and_si256(delta, mask);
    if (odd) { o = _mm256_add_epi8(o,delta); e = _mm256_sub_epi8(e,left1(delta)); }
    else { e = _mm256_add_epi8(e,delta); o = _mm256_sub_epi8(o,delta); }
  }
  _mm256_store_si256(reinterpret_cast<__m256i*>(ev),e);
  _mm256_store_si256(reinterpret_cast<__m256i*>(od),o);
  for (size_t i=0;i<n;++i) a[start+i]=(i&1 ? od : ev)[i/2];
}
}}}
