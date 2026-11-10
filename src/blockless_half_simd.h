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

template<bool Backward = false, bool Inject = false>
inline void power_simd(State& a, unsigned count, size_t start = 0)
{
  const size_t n = a.size() - start;
  if (!count) return;
  if (n < 2) { if constexpr (Inject) a[start] += count; return; }
  uint64_t mass = Inject ? count : 0;
  for (size_t i=start; i<a.size(); ++i) mass += a[i];
  if (n > 64 || mass > 127)
  {
    for (unsigned t=0; t<count; ++t)
      { if constexpr (Inject) ++a[start];
        if constexpr (Backward) inverse(a,start); else sweep(a,start); }
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
  const auto first = _mm256_setr_epi8(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  for (unsigned tick=0; tick<=end; ++tick)
  {
    const bool odd = (tick ^ (Backward ? last : 0)) & 1;
    if constexpr (Inject)
      if (!odd && tick <= 2*(count-1)) e = _mm256_add_epi8(e,first);
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

namespace snaperz { namespace blockless { namespace half {
// Construct D using the inverse wavefront itself. Each completed right-hand
// coordinate extends a suffix sum for its logical sweep. Particle rank t+1
// in row t is selected by that suffix crossing mass-t. This avoids storing
// all rows or serializing the inverse sweeps for particle-position queries.
inline bool retiming_simd(const State& source, State& dest)
{
  const unsigned n = source.size()-1;
  uint64_t total_mass = 0;
  for (size_t i=1;i<source.size();++i) total_mass += source[i];
  if (n > 64 || total_mass > 127) return false;
  const unsigned mass = total_mass;
  std::fill(dest.begin(),dest.end(),0);
  dest[0]=source[0];
  if (!mass) return true;
  for (unsigned i=1;i<=n;++i) if(source[i]) { ++dest[i]; break; }
  if (mass==1 || n<2) { if(n==1) dest[1]=mass; return true; }
  alignas(32) uint8_t ev[32]={},od[32]={},ix[32];
  for(unsigned i=0;i<n;++i) (i&1?od:ev)[i/2]=source[i+1];
  for(unsigned i=0;i<32;++i) ix[i]=i;
  auto e=_mm256_load_si256(reinterpret_cast<const __m256i*>(ev));
  auto o=_mm256_load_si256(reinterpret_cast<const __m256i*>(od));
  const auto idx=_mm256_load_si256(reinterpret_cast<const __m256i*>(ix));
  const auto zero=_mm256_setzero_si256(),one=_mm256_set1_epi8(1);
  auto se=zero,so=zero,ce=zero,co=zero;
  const int last=n-2, count=mass-1;
  for(int tick=0;tick<=2*(count-1)+last;++tick)
  {
    const bool odd=(tick^last)&1;
    auto k=_mm256_add_epi8(_mm256_add_epi8(idx,idx),_mm256_set1_epi8(odd));
    auto distance=_mm256_sub_epi8(_mm256_set1_epi8(last),k);
    auto mask=_mm256_cmpgt_epi8(_mm256_set1_epi8(last+1),k);
    if(tick<last) mask=_mm256_and_si256(mask,_mm256_cmpgt_epi8(_mm256_set1_epi8(tick+1),distance));
    if(tick>2*(count-1)) mask=_mm256_and_si256(mask,_mm256_cmpgt_epi8(distance,_mm256_set1_epi8(tick-2*(count-1)-1)));
    auto l=odd?o:e, r=odd?shr1(e):o;
    auto delta=_mm256_sub_epi8(one,_mm256_and_si256(_mm256_cmpeq_epi8(r,zero),l));
    delta=_mm256_and_si256(delta,_mm256_andnot_si256(_mm256_cmpeq_epi8(l,zero),mask));
    auto right=_mm256_sub_epi8(r,delta);
    auto suffix=odd?shr1(se):so;
    suffix=_mm256_andnot_si256(_mm256_cmpeq_epi8(k,_mm256_set1_epi8(last)),suffix);
    auto sum=_mm256_add_epi8(suffix,right);
    auto target=_mm256_sub_epi8(_mm256_set1_epi8(int(mass)-1-(tick-last+int(odd))/2),idx);
    auto hit=_mm256_and_si256(mask,_mm256_and_si256(_mm256_cmpgt_epi8(target,suffix),
      _mm256_andnot_si256(_mm256_cmpgt_epi8(target,sum),one)));
    if(odd) {o=_mm256_add_epi8(o,delta);e=_mm256_sub_epi8(e,left1(delta));so=sum;ce=_mm256_add_epi8(ce,left1(hit));}
    else
    {
      e=_mm256_add_epi8(e,delta);o=right;se=sum;co=_mm256_add_epi8(co,hit);
      // Gate zero completes the head too; its remaining suffix contains all
      // other coordinates. Only its lane may contribute this head particle.
      if(tick>=last && tick<=2*(count-1)+last)
        if((_mm256_movemask_epi8(_mm256_cmpgt_epi8(target,sum)) & 1)) ++dest[1];
    }
  }
  _mm256_store_si256(reinterpret_cast<__m256i*>(ev),ce);
  _mm256_store_si256(reinterpret_cast<__m256i*>(od),co);
  for(unsigned i=0;i<n;++i) dest[i+1]+=(i&1?od:ev)[i/2];
  return true;
}
}}}
