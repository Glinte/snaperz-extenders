#pragma once
// Finite AVX2/AVX512 wavefront on up to 64/128 byte coordinates. Gate k of sweep t
// executes at tick 2*t+k (or 2*t+n-2-k for inverse sweeps). Gates at a
// tick are disjoint; warm-up/drain masks make the returned state synchronous.
// Byte arithmetic is safe only when total mass (including injections) <=127;
// larger inputs retain the scalar implementation. Narrow registers are preferred
// when they fit. These helpers expect a valid start index and a nonempty state.
namespace snaperz { namespace blockless { namespace half {
struct Avx2
{
  using Vector = __m256i;
  static constexpr unsigned lanes = 32;
  static Vector load(const uint8_t* p) { return _mm256_load_si256(reinterpret_cast<const Vector*>(p)); }
  static void store(uint8_t* p, Vector v) { _mm256_store_si256(reinterpret_cast<Vector*>(p), v); }
  static Vector zero() { return _mm256_setzero_si256(); }
  static Vector set(int x) { return _mm256_set1_epi8(static_cast<char>(x)); }
  static Vector add(Vector a, Vector b) { return _mm256_add_epi8(a, b); }
  static Vector sub(Vector a, Vector b) { return _mm256_sub_epi8(a, b); }
  static Vector eq(Vector a, Vector b) { return _mm256_cmpeq_epi8(a, b); }
  static Vector gt(Vector a, Vector b) { return _mm256_cmpgt_epi8(a, b); }
  static Vector both(Vector a, Vector b) { return _mm256_and_si256(a, b); }
  static Vector without(Vector a, Vector b) { return _mm256_andnot_si256(a, b); }
  static uint64_t bits(Vector v) { return static_cast<uint32_t>(_mm256_movemask_epi8(v)); }
  static Vector right(Vector v) { return shr1(v); }
  static Vector left(Vector v)
  {
    const auto lo = _mm256_permute2x128_si256(v, v, 0x08);
    return _mm256_alignr_epi8(v, lo, 15);
  }
};

#if SNAPERZ_BLOCKLESS_HAVE_AVX512
struct Avx512
{
  using Vector = __m512i;
  static constexpr unsigned lanes = 64;
  static Vector load(const uint8_t* p) { return _mm512_load_si512(p); }
  static void store(uint8_t* p, Vector v) { _mm512_store_si512(p, v); }
  static Vector zero() { return _mm512_setzero_si512(); }
  static Vector set(int x) { return _mm512_set1_epi8(static_cast<char>(x)); }
  static Vector add(Vector a, Vector b) { return _mm512_add_epi8(a, b); }
  static Vector sub(Vector a, Vector b) { return _mm512_sub_epi8(a, b); }
  static Vector eq(Vector a, Vector b) { return _mm512_maskz_set1_epi8(_mm512_cmpeq_epi8_mask(a, b), -1); }
  static Vector gt(Vector a, Vector b) { return _mm512_maskz_set1_epi8(_mm512_cmpgt_epi8_mask(a, b), -1); }
  static Vector both(Vector a, Vector b) { return _mm512_and_si512(a, b); }
  static Vector without(Vector a, Vector b) { return _mm512_andnot_si512(a, b); }
  static uint64_t bits(Vector v) { return _mm512_movepi8_mask(v); }
  static Vector right(Vector v) { return wide::shr1(v); }
  static Vector left(Vector v)
  {
    const auto idx = _mm512_setr_epi64(0, 0, 1, 2, 3, 4, 5, 6);
    const auto prev = _mm512_maskz_permutexvar_epi64(0xFE, idx, v);
    return _mm512_or_si512(_mm512_slli_epi64(v, 8), _mm512_srli_epi64(prev, 56));
  }
};
#endif

template <class V, bool Backward = false, bool Inject = false>
inline void power_vector(State& a, unsigned count, size_t start = 0)
{
  const size_t n = a.size() - start;
  if (!count) return;
  if (n < 2)
  {
    if constexpr (Inject) a[start] += count;
    return;
  }
  uint64_t mass = Inject ? count : 0;
  for (size_t i=start; i<a.size(); ++i) mass += a[i];
  if (n > 2*V::lanes || mass > 127)
  {
    for (unsigned t = 0; t < count; ++t)
    {
      if constexpr (Inject) ++a[start];
      if constexpr (Backward) inverse(a, start); else sweep(a, start);
    }
    return;
  }
  alignas(64) uint8_t ev[V::lanes] = {}, od[V::lanes] = {};
  for (size_t i=0; i<n; ++i)
    if (i & 1) od[i/2] = a[start+i]; else ev[i/2] = a[start+i];
  auto e = V::load(ev);
  auto o = V::load(od);
  const auto zero = V::zero(), one = V::set(1);
  alignas(64) uint8_t indices[V::lanes];
  for (unsigned i=0;i<V::lanes;++i) indices[i]=2*i;
  const auto idx = V::load(indices);
  const unsigned last = unsigned(n)-2;
  const unsigned end = 2*(count-1)+last;
  alignas(64) uint8_t first_lane[V::lanes] = {1};
  const auto first = V::load(first_lane);
  for (unsigned tick=0; tick<=end; ++tick)
  {
    const bool odd = (tick ^ (Backward ? last : 0)) & 1;
    if constexpr (Inject)
      if (!odd && tick <= 2*(count-1)) e = V::add(e, first);
    auto l = odd ? o : e;
    auto r = odd ? V::right(e) : o;
    auto k = V::add(idx, V::set(odd));
    auto mask = V::gt(V::set(last+1), k);
    if (tick < last || tick > 2*(count-1))
    {
      auto distance = Backward ? V::sub(V::set(last), k) : k;
      if (tick < last)
        mask = V::both(mask, V::gt(V::set(tick+1), distance));
      if (tick > 2*(count-1))
        mask = V::both(mask, V::gt(distance, V::set(tick-2*(count-1)-1)));
    }
    typename V::Vector delta;
    if constexpr (Backward)
    {
      delta = V::sub(one, V::both(V::eq(r, zero), l));
      mask = V::without(V::eq(l, zero), mask);
    }
    else
      delta = V::add(V::both(V::eq(l, one), r), V::gt(l, one));
    delta = V::both(delta, mask);
    if (odd) { o = V::add(o, delta); e = V::sub(e, V::left(delta)); }
    else { e = V::add(e, delta); o = V::sub(o, delta); }
  }
  V::store(ev, e);
  V::store(od, o);
  for (size_t i=0;i<n;++i) a[start+i]=((i&1) ? od[i/2] : ev[i/2]);
}
// Construct D using the inverse wavefront itself. Each completed right-hand
// coordinate extends a suffix sum for its logical sweep. Particle rank t+1
// in row t is selected by that suffix crossing mass-t. This avoids storing
// all rows or serializing the inverse sweeps for particle-position queries.
template <class V>
inline bool retiming_vector(const State& source, State& dest)
{
  const unsigned n = source.size()-1;
  uint64_t total_mass = 0;
  for (size_t i=1;i<source.size();++i) total_mass += source[i];
  if (n > 2*V::lanes || total_mass > 127) return false;
  const unsigned mass = total_mass;
  std::fill(dest.begin(), dest.end(), 0);
  dest[0]=source[0];
  if (!mass) return true;
  for (unsigned i=1;i<=n;++i) if (source[i]) { ++dest[i]; break; }
  if (mass==1 || n<2) { if (n==1) dest[1]=mass; return true; }
  alignas(64) uint8_t ev[V::lanes]={}, od[V::lanes]={}, ix[V::lanes];
  for (unsigned i=0;i<n;++i)
    if (i & 1) od[i/2] = source[i+1]; else ev[i/2] = source[i+1];
  for (unsigned i=0;i<V::lanes;++i) ix[i]=i;
  auto e=V::load(ev);
  auto o=V::load(od);
  const auto idx=V::load(ix);
  const auto zero=V::zero(), one=V::set(1);
  // Suffix accumulators and output particle counts, indexed by gate parity.
  auto se = zero, so = zero, ce = zero, co = zero;
  const int last=n-2, count=mass-1;
  for (int tick=0;tick<=2*(count-1)+last;++tick)
  {
    const bool odd=(tick^last)&1;
    auto k=V::add(V::add(idx, idx), V::set(odd));
    auto distance=V::sub(V::set(last), k);
    auto mask=V::gt(V::set(last+1), k);
    if (tick<last) mask=V::both(mask, V::gt(V::set(tick+1), distance));
    if (tick>2*(count-1)) mask=V::both(mask, V::gt(distance, V::set(tick-2*(count-1)-1)));
    auto l=odd?o:e, r=odd?V::right(e):o;
    auto delta=V::sub(one, V::both(V::eq(r, zero), l));
    delta=V::both(delta, V::without(V::eq(l, zero), mask));
    auto right=V::sub(r, delta);
    auto suffix=odd?V::right(se):so;
    suffix=V::without(V::eq(k, V::set(last)), suffix);
    auto sum=V::add(suffix, right);
    // At gate k=2*i+odd the completed row is
    // t=(tick-last+odd)/2+i+1, so its particle has right-rank mass-t.
    auto target=V::sub(V::set(int(mass)-1-(tick-last+int(odd))/2), idx);
    auto hit=V::both(mask, V::both(V::gt(target, suffix),
      V::without(V::gt(target, sum), one)));
    if (odd)
    {
      o = V::add(o, delta);
      e = V::sub(e, V::left(delta));
      so = sum;
      ce = V::add(ce, V::left(hit));
    }
    else
    {
      e = V::add(e, delta);
      o = right;
      se = sum;
      co = V::add(co, hit);
      // Gate zero completes the head too; its remaining suffix contains all
      // other coordinates. Only its lane may contribute this head particle.
      if (tick>=last && tick<=2*(count-1)+last)
        if ((V::bits(V::gt(target, sum)) & 1)) ++dest[1];
    }
  }
  V::store(ev, ce);
  V::store(od, co);
  for (unsigned i=0;i<n;++i) dest[i+1]+=((i&1) ? od[i/2] : ev[i/2]);
  return true;
}
template <bool Backward = false, bool Inject = false>
inline void power_simd(State& a, unsigned count, size_t start = 0)
{
#if SNAPERZ_BLOCKLESS_HAVE_AVX512
  if (a.size()-start > 64)
    return power_vector<Avx512, Backward, Inject>(a, count, start);
#endif
  power_vector<Avx2, Backward, Inject>(a, count, start);
}

inline bool retiming_simd(const State& source, State& dest)
{
#if SNAPERZ_BLOCKLESS_HAVE_AVX512
  if (source.size() > 65) return retiming_vector<Avx512>(source, dest);
#endif
  return retiming_vector<Avx2>(source, dest);
}
}}}
