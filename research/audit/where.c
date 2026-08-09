#include <stdio.h>
#include <string.h>
static inline void Astep(int *x,int L){for(int i=0;i<L-1;i++){int a=x[i];if(!a)continue;if(a==1){x[i]=1+x[i+1];x[i+1]=0;}else{x[i]=a-1;x[i+1]++;}}}
static int good(int u,int v){return u==0?v%2==0:(u%2==1&&v%2==1);}
static int ell(const int*x,int L){int l=0;while(l+1<L&&good(x[l],x[l+1]))l+=2;return l;}
int main(void){const int L=69;int x[69],x0[69];for(int i=0;i<L;i++)x[i]=x0[i]=1;long n=0;
for(;;){Astep(x,L);if(x[0]==1){int e=ell(x,L);if(e<56)printf("  checkpoint #%ld (outer step %ld, offset %ld): ell=%d\n",n,n/112,n%112,e);n++;}
if(!memcmp(x,x0,sizeof x))break;}return 0;}
