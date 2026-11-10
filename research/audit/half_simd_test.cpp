// Deterministic arbitrary-state checks against sequential gates, including
// warm-up/drain, lane boundaries, fallback bounds and complete counts.
#include "../../src/blockless_half.h"
#include <random>
#include <cstdio>
using namespace snaperz::blockless::half;
static int powers_and_counts(){std::mt19937 rng(123);
for(unsigned n=1;n<=66;++n) for(unsigned c: {0,1,2,15,31,435}) for(unsigned start: {0,1}) {
State a(n); for(unsigned i=0;i<n;++i) ++a[rng()%n];
a[0]+=n;
for(bool back:{false,true}) {auto x=a,y=a;
for(unsigned i=0;i<c;++i) {if(back) inverse(x,start);else sweep(x,start);}
if(back)power_simd<true>(y,c,start);else power_simd(y,c,start);
if(x!=y){printf("FAIL %u %u %u %d\n",n,c,start,back);return 1;}
}}
for(unsigned n=17;n<=44;++n){auto x=run<false>(n),y=run(n);if(x.result.total!=y.result.total||x.result.p!=y.result.p||x.steps!=y.steps){printf("count fail %u\n",n);return 1;}}
puts("PASS powers and counts"); return 0;}

static int reflection_and_injection(){std::mt19937 rng(321);
for(unsigned n=1;n<=65;++n)for(unsigned mass=0;mass<=128;++mass){State a(n+1),b(n+1),c(n+1);a[0]=9;b[0]=9;
for(unsigned i=0;i<mass;++i)++a[1+rng()%n];auto x=a;
for(unsigned r=1;r<=mass;++r){unsigned total=0;for(unsigned i=1;i<=n;++i){total+=x[i];if(total>=r){++b[i];break;}}if(r<mass)inverse(x,1);}
bool ok=retiming_simd(a,c);if(ok && b!=c){printf("fail n=%u mass=%u\n",n,mass);return 1;}
if(ok!=(n<=64&&mass<=127)) return 2;
for(unsigned count:{0,1,3,31}){auto u=a,v=a;for(unsigned i=0;i<count;++i){++u[0];sweep(u);}power_simd<false,true>(v,count);if(u!=v){puts("injection fail");return 3;}}
}puts("PASS retiming and injection"); return 0;}

int main() { if (powers_and_counts()) return 1; return reflection_and_injection(); }
