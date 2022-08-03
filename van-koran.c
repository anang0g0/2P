//#include "1331.h"
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "global.h"
#include "struct.h"
//#include "chash.c"
#include "debug.c"
//#include "lu.c"
//#include "invm.c"


//#define N 256
//#define K 32
#define EN T

unsigned short vb[K][N] = {0};
unsigned short syn[K] = {0};
int a[T][T + 1] = {0};

#define O 257 // 6859 //1331 //2197,4913,6859
//#define K 4
#define P 257

// sagemath上での原始多項式
unsigned short pp[4][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}};
// {0,0,9,2}, {1,0,11,2}, {1,0,16,3}, {1,0,15,2};
// GF(11^3,13^3,17^3,19^3)
// unsigned short ff[2][7]={{1,0,0,0,0,2,0,2},{0,0,1,0,0,0,1,2}}; //GF(3^7,5^5)

// unsigned short O]={0},O]={0};
// int N =0,M=0;



// OP型からベクトル型への変換
vec o2v(OP f)
{
    vec a = {0};
    int i;

    //#pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (f.t[i].a > 0 && f.t[i].n < DEG)
            a.x[f.t[i].n] = f.t[i].a;
    }

    return a;
}

//ベクトル型からOP型への変換
OP v2o(vec a)
{
    int i, j = 0;
    OP f = {0};

    //#pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            f.t[j].n = i;
            f.t[j++].a = a.x[i];
        }
    }

    return f;
}

void op_print_raw(const OP f)
{
    puts("op_print_raw:");
    for (int i = 0; i < DEG; i++)
    {
        if (f.t[i].a > 0)
            printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
    }
}




int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return x*y%P; //((x + y - 2) % (ORD - 1)) + 1;
}

int mltn(int n, int x) {
    int ret = 1;
    while (n > 0) {
        if (n & 1) ret = mlt(ret , x) ;  // n の最下位bitが 1 ならば x^(2^i) をかける
        x = mlt(x , x);
        n >>= 1;  // n を1bit 左にずらす
    }
    return ret
    ;
}

int mlt2(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = mlt(i, x);

  return i;
}



//多項式の次数(default)
int deg(vec a)
{
    int i, n = 0;

    //#pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
            n = i;
    }

    return n;
}

//配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
    int i;
    vec v = {0};

    for (i = 0; i < n; i++)
    {
        v.x[n - 1 - i] = c[i];
    }

    return v;
}

//配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
    OP g;
    vec a;
    int i;

    memset(c, 0, sizeof(c));
    // memcpy (c, f, n);
    for (i = 0; i < n; i++)
    {
        c[i] = f[i];
        // printf("%d,",f[i]);
    }
    // exit(1);
    a = Setvec(n);

    g = v2o(a);

    return g;
}

//多項式を表示する(default)
void printpol(vec a)
{
    int i, n;

    n = deg(a);

    // printf ("baka\n");
    assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            printf("%u", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            // if (i > 0)
            printf("+");
        }
    }
    //  printf("\n");

    return;
}

//多項式の代入値
unsigned short
trace(OP f, unsigned short x)
{
    int i, d;
    unsigned short u = 0;

    d = deg(o2v(f));

    for (i = 0; i < d + 1; i++)
    {
        u += mlt(f.t[i].a, mlt2(f.t[i].n, x));
    }

    return u;
}


//多項式の代入値
unsigned short
xtrace(OP f, unsigned short x)
{
    int i, d;
    unsigned short u = 0, v = 1;

    d = deg(o2v(f));
    // printpol(o2v(f));
    // printf(" =ff\n");

    for (i = 0; i < d + 1; i++)
    {
        v = 0;
        if (f.t[i].a > 0)
        {
            v = 1;
            //for (int j = 0; j < f.t[i].n; j++)
            {
                v = mltn(f.t[i].n , x);
            }
            u += mlt(v , f.t[i].a);

            // printf("\nv=%d",v);
        }
        //u = (u ^ v);
    }
    // printf("u=%d\n",u%O);

    return u;
}


OP minus(OP f)
{
    unsigned int i, j;

    j = deg(o2v(f));
    for (i = 0; i < j + 1; i++)
    {
        if (f.t[i].a > 0)
            f.t[i].a = P - f.t[i].a;
    }

    return f;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
    int i, k;
    oterm t = {0};

    // k = deg (o2v (f));
    for (i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.t[i].a > 0)
        {
            t.n = f.t[i].n;
            t.a = f.t[i].a % P;
        }
    }

    return t;
}

// OP型を正規化する
OP conv(OP f)
{
    vec v = {0};
    OP g = {0};

    v = o2v(f);
    g = v2o(v);

    return g;
}

// 20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
    vec a = {0}, b = {0}, c = {0};
    int i, j, k, l = 0;
    OP h = {0}, f2 = {0}, g2 = {0};

    // for(i=0;i<257;i++)
    //  printf("%d %d %d %d %d\n",i,f.t[i].a,f.t[i].n,g.t[i].a,g.t[i].n);

    //  exit(1);
    f = conv(f);
    g = conv(g);

    a = o2v(f);
    // exit(1);
    b = o2v(g);

    j = deg(o2v(f));
    l = deg(o2v(g));
    printpol(o2v(f));
    printf(" =f_& %d\n", j);
    printpol(o2v(g));
    printf(" =g_& %d\n", l);
    // exit(1);

    if (j >= l)
    {
        k = j + 1;
    }
    else
    {

        k = l + 1;
    }
    // for(i=0;i<k;i++)
    // printf("%d %d\n",i,b.x[i]);
    //   exit(1);

    for (i = 0; i < k; i++)
    {
        // if(a.x[i]>b.x[i])
        c.x[i] = (a.x[i] + b.x[i]) % P;
        if (a.x[i] == b.x[i])
            c.x[i] = a.x[i] * 2 % P;
    }
    //
    h = v2o(c);
    printpol(o2v(h));
    printf(" =====in oadd\n");

    return h;
}

//多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{

    // assert (op_verify (f));
    int i, k, j;
    OP h = {0};
    vec test;
    unsigned short n;

    // f=conv(f);
    k = deg(o2v(f));
    j = 0;
    for (i = 0; i < k + 1; i++)
    {
        h.t[i].n = f.t[i].n + t.n;
        h.t[i].a = (f.t[i].a * t.a) % P;
    }

    // h=conv(h);
    // assert (op_verify (h));
    return h;
}

//多項式の掛け算
OP omul(OP f, OP g)
{
    f = conv(f);
    g = conv(g);
    // assert (op_verify (f));
    // assert (op_verify (g));
    int i, count = 0, k, l;
    oterm t = {0};
    OP h = {0}, e = {0}, r = {0};
    vec c = {0};

    k = deg(o2v(f));
    l = deg(o2v(g));
    if (l > k)
    {
        k = l;
    }

    for (i = 0; i < k + 1; i++)
    {
        t = g.t[i];
        e = oterml(f, t);
        h = oadd(h, e);
    }
    // assert (op_verify (h));
    return h;
}

OP osub(OP f, OP g)
{
    vec a = {0}, b = {0}, d = {0};
    int i, k, l, m;
    OP ans = {0};

    a = o2v(f);
    b = o2v(g);
    l = deg(a);
    m = deg(b);
    if (l >= m)
    {
        k = l;
    }
    else
    {
        k = m;
    }
    for (i = 0; i < k + 1; i++)
    {
        if (a.x[i] >= b.x[i])
        {
            d.x[i] = a.x[i] - b.x[i];
        }

        else
        {
            d.x[i] = (P - (a.x[i] + b.x[i]));
        }
        /*
        if(d.x[i]<0){
          printf("%d\n",d.x[i]);
          d.x[i]+=P;
        }
        */
    }

    ans = v2o(d);

    return ans;
}

OP confer(OP f, int a)
{
    vec r;
    int n, i;
    OP g;

    r = o2v(f);
    n = deg(r);
    for (i = 0; i < n + 1; i++)
        r.x[i] = (r.x[i] * a) % P;
    g = v2o(r);

    return g;
}

int oequ(OP f, OP g)
{
    vec v, x;
    int i, flg = 0;

    v = o2v(f);
    x = o2v(g);
    for (i = 0; i < 512; i++)
    {
        if (v.x[i] != x.x[i])
            return -1;
    }

    return 0;
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
    int i;

    for (i = 0; i < N; i++)
    {
        if ((a * i) % P == b)
            break;
    }
    return i;
}


// invert of integer
unsigned short inv(unsigned short a, unsigned short n)
{
    unsigned short d, x, s, q, r, t, gcd;
    d = n;
    x = 0;
    s = 1;
    while (a != 0)
    {
        q = d / a;
        r = d % a;
        d = a;
        a = r;
        t = (x - q * s) % P;
        x = s;
        s = t;
    }
    gcd = d % P;

    return ((x + n) % (n / d)) % P;
}

//有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
    int i;

    if (a == 0)
        return -1;

    for (i = 0; i < N; i++)
    {
        if ((a * i) % P == 1)
            return (unsigned short)i;
    }

    printf("no return \n");
    //  exit (1);
}

//有限体の元の逆数
unsigned short
iinv(unsigned short a)
{
    int i;

    if (a == 0)
        return -1;

    for (i = 0; i < P; i++)
    {
        if (mlt(a , i)  == 1)
            return (unsigned short)i;
    }

    printf("no return \n");
    //  exit (1);
}

//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
    oterm tt = {0}, s = {
                        0};

    tt = LT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a * inv(t.a, P)) % P;
        s.n = tt.n;
    }
    else
    {
        printf("debug in LTdiv\n");
        exit(1);
    }

    return s;
}

//多項式の剰余を取る
OP omod(OP f, OP g)
{
    int i = 0, j, n, k;
    OP h = {0}, e = {
                    0};
    oterm a, b = {0}, c = {0};

    n = LT(g).n;

    //  assert (("baka^\n", LT (f).n != 0));

    //  assert (("baka(A)\n", LT (g).n != 0));

    if (LT(f).n < LT(g).n)
    {
        //    exit(1);
        return f;
    }

    // printf ("in omod\n");
    // exit(1);

    k = LT(g).n;
    b = LT(g);
    OP ll;

    assert(("double baka\n", b.a != 0 && b.n != 0));
    while (LT(f).n > 0 && LT(g).n > 0)
    {

        c = LTdiv(f, b);
        h = oterml(g, c);
        printpol(o2v(f));
        printf("======f_before_omod\n");
        printpol(o2v(h));
        printf("======h_before_omod\n");
        f = osub(f, (h));
        printpol(o2v((h)));
        printf(" =====h_minus_omod\n");
        printpol(o2v(f));
        printf(" =====f_after_omod\n");
        // exit(1);
        if (deg(o2v(f)) == 0 || deg(o2v(g)) == 0)
        {
            //      printf("blake1\n");
            break;
        }

        if (c.n == 0 || b.n == 0)
            break;
        if (LT(f).a == 4 && deg(o2v(f)) == 0)
            exit(1);
    }
    printpol(o2v(f));
    printf("\n");
    // exit(1);

    return f;
}

//項の数
int terms(OP f)
{
    int i, count = 0;

    for (i = 0; i < DEG; i++)
        if (f.t[i].a > 0)
            count++;

    return count;
}

//モニック多項式にする
OP coeff(OP f)
{
    int i, j, k;
    vec a, b;
    oterm t;

    t = LT(f);
    // f = conv(f);
    k = deg(o2v(f)) + 1;
    for (i = 0; i < k; i++)
        f.t[i].a = (f.t[i].a * inv(t.a, P)) % P;

    return f;
}

OP kei(unsigned short u, OP g)
{
    vec v = {0};
    int j, i = 0;
    OP f = {0};

    v = o2v(g);
    i = deg(v);
    printpol(v);
    printf("\n");
    for (j = 0; j < i + 1; j++)
        v.x[j] = (v.x[j] * u) % P;
    f = v2o(v);
    printpol(v);
    printf("\n");
    // exit(1);

    return f;
}

//多項式の商を取る
OP odiv(OP f, OP g)
{

    f = conv(f);
    g = conv(g);
    // assert (op_verify (f));
    // assert (op_verify (g));
    int i = 0, j, n, k;
    OP h = {0}, e = {0}, tt = {0};
    oterm a, b = {0}, c = {0};

    printpol(o2v(f));
    printf("\n");
    printpol(o2v(g));
    printf("\n");
    // exit(1);

    if (LT(f).n == 0 && LT(g).a == 0)
    {
        printf("baka^\n");
        // return f;
        exit(1);
    }
    if (LT(g).a == 0)
    {
        printf("a==0\n");
        // print_trace ();
        exit(1);
    }
    if (LT(g).n == 0 && LT(g).a > 1)
        return g; // coeff (f);

    k = deg(o2v(g));
    b = LT(g);
    if (b.a == 1 && b.n == 0)
        return f;
    if (b.a == 0 && b.n == 0)
    {
        printf("baka in odiv\n");
        exit(1);
    }
    if (deg(o2v(f)) < deg(o2v(g)))
    {
        return f;
        //  a=LT(f);
    }
    OP null = {0};
    i = 0;
    while (LT(f).n > 0 && LT(g).n > 0)
    {
        c = LTdiv(f, b);
        c.a = c.a % P;
        assert(c.n < DEG);
        tt.t[i] = c;
        i++;

        h = oterml(g, c);
        f = oadd(f, minus(h));
        printpol(o2v(h));
        printf(" ===h in_odiv\n");
        printpol(o2v(f));
        printf(" ===f in_odiv\n");
        if (deg(o2v(f)) == 0 && LT(f).a > 0)
        {
            return f;
            // exit(1);
        }
        else
        {
            // return g;
        }
        if (deg(o2v(f)) == 0 || deg(o2v(g)) == 0)
        {
            printf("blake2\n");
            break;
        }
        if (oequ(f, g) == 0)
        {
            printpol(o2v(tt));
            printf("\n");
            break;
            // exit(1);
        }
        if (c.n == 0)
            break;
    }

    // tt は逆順に入ってるので入れ替える
    OP ret = {0};
    int tt_terms = terms(tt);
    for (i = 0; i < tt_terms; i++)
    {
        ret.t[i] = tt.t[tt_terms - i - 1];
    }
    ret = conv(ret);
    printpol(o2v(ret));
    printf("\n");
    //  exit(1);

    // assert (op_verify (ret));
    return ret;
}

void synd(unsigned short zz[T])
{
    int i,k,count=0;

    memset(syn,0,sizeof(syn));
    for (i = 1; i < N ; i++){
        if(zz[i]>0){
            printf("%d,",i);
            
            for(k=0;k<K;k++){
            syn[k] ^= vb[k][i];
            printf("@%d,",vb[k][i]);
            }
            for(k=0;k<K;k++)
            printf("\n%d,%d\n",i,syn[k]);
            printf("\n");
            
        }
    }
//exit(1);

}
/*
 * Author: Hiroyuki Chishiro
 * License: 2-Clause BSD
 */

OP hyoe[1331] = {0};
unsigned short sue[1331] = {0};
unsigned short sie[1331] = {0};
void tas()
{
    int i, j, k;
    unsigned short ccp[4] = {0, 0, 9, 2};
    OP g = {0}, h = {0}, f = {0}, m = {0};
    vec v = {0}, x = {0}, w = {0};

    v.x[1] = 9;
    v.x[0] = 2;
    m = v2o(v);
    printpol(v);
    printf("\n");
    g = m;
    w.x[1] = 1;
    printpol(w);
    printf("\n");

    h = v2o(w);
    // g=omul(g,h);
    printpol(o2v(g));
    printf("\n");
    //  exit(1);

    memset(v.x, 0, sizeof(v.x));
    for (i = 0; i < 2; i++)
    {
        v.x[0] = i;
        hyoe[i] = v2o(v);
    }
    memset(v.x, 0, sizeof(v.x));
    v.x[1] = 1;
    hyoe[2] = v2o(v);
    memset(v.x, 0, sizeof(v));
    v.x[2] = 1;
    hyoe[3] = v2o(v);
    x.x[3] = 1;
    f = v2o(x);
    hyoe[4] = g;
    g = omul(g, h);
    printpol(o2v(g));
    hyoe[5] = g;
    // exit(1);

    for (i = 6; i < 1331; i++)
    {
        g = omul(h, g);
        if (deg(o2v(g)) == 3)
        {
            g = oadd(g, kei(LT(g).a, m));
            g = omod(g, f);
            printpol(o2v(g));
            printf("====3\n");
            // exit(1);
            printpol(o2v(m));
            printf("\n");
        }
        //  exit(1);
        //  g=oadd(g,

        printpol(o2v(g));
        printf("\n");
        hyoe[i] = g;
        printpol(o2v(g));
        printf("\n");
    }
    // exit(1);
    for (i = 0; i < 2; i++)
        sue[i] = i;
    for (i = 2; i < 1331; i++)
        sue[i] = trace(hyoe[i], P);
    for (i = 0; i < 1331; i++)
        sie[sue[i]] = i;
}


void kakeru(MTX a,unsigned short e[T]){
int i,j,k,c[T]={0};

for(i=0;i<T;i++){
    for(j=0;j<T;j++)
    c[i]^=mlt(a.x[j][i],e[j]);
}
for(i=0;i<T;i++)
printf("%d,",c[i]);
printf("\n");

}


//#define EN K
int g2(void)
{
  unsigned short a[3][3 + 1] = {
{123,231,11,1},
{51,12,4,77},
{11,22,33,44}
  };
  
  unsigned short p, d;
  int i, j, k;



//printf("%d\n",hiku(11,121));
//exit(1);


  for (i = 0; i < EN; i++) {
    p = a[i][i];
    if(p<0)
    p+=P;
    for (j = 0; j < (EN + 1); j++) {
      a[i][j] = mlt(a[i][j] , iinv(p));
    }

    for (j = 0; j < EN; j++) {
      if (i != j) {
        d = a[j][i];
        if(d<0)
        d+=P;
        for (k = i; k < (EN + 1); k++) {
          a[j][k] = (a[j][k] - mlt(d , a[i][k]));
        }
      }
    }
  }

  for (i = 0; i < EN; i++) {
    printf("x%d = %d\n", i + 1, a[i][EN]);
  }

  return 0;
}


//連立方程式の解
#define EN 3
int renritu(MTX z)
{
  unsigned short a[EN][EN + 1] = {0
  };
  
  unsigned short p, d;
  int i, j, k;

for(i=0;i<T;i++){
    for(j=0;j<T+1;j++){
    a[i][j]=z.x[i][j];
    printf("%d,",a[i][j]);
    }
    printf("\n");
}

//printf("%d\n",hiku(11,121));
//exit(1);


  for (i = 0; i < EN; i++) {
    p = a[i][i];

    for (j = 0; j < (EN + 1); j++) {
      a[i][j] = mlt(a[i][j] , inv(p,P));
    }

    for (j = 0; j < EN; j++) {
      if (i != j) {
        d = a[j][i];

        for (k = i; k < (EN + 1); k++) {
          a[j][k] = a[j][k] ^ mlt(d , a[i][k]);
        }
      }
    }
  }

  for (i = 0; i < EN; i++) {
    printf("x%d = %d\n", i + 1, a[i][EN]);
  }

  return 0;
}


//逆行列
MTX gyaku(MTX z)
{
    //  int a[EN][EN + 1] = { 0 };
    int p, d;
    int i, j, k;
    MTX V={0};

    int a[K][K*2]={0},ken[K][K]={0},zan[K][K]={0};
    int l=0;
for(i=0;i<K;i++){
    for(j=0;j<K;j++)
    a[i][j]=z.x[i][j];
}

for(i=0;i<K;i++){
    for(j=0;j<K;j++)
    ken[i][j]=a[i][j];
}

for(i=0;i<K;i++){
    for(j=K;j<K*2;j++)
    a[i][K+i]=1;
}


for(i=0;i<K;i++){
    for(j=0;j<K*2;j++)
    printf("%d,",a[i][j]);
    printf("\n");
}
//exit(1);

 for(k=0;k<K;k++){   
    d=a[k][k];
    if(d<0)
    d=P+d;
    for(j=0;j<K*2;j++){
        a[k][j]=mlt(a[k][j],iinv(d));
    }
 
    for(i=k;i<K-1;i++){
        d=a[i+1][k];
        if(d<0)
        d=P+d;
        printf("d=%d\n",d);
        for(j=0;j<T*2;j++)
        
        for(j=0;j<K*2;j++){
        a[i+1][j]-=mlt(a[k][j],d);
        if(a[i+1][j]<0)
        a[i+1][j]+=P;
        if(a[i+1][j]>P)
        a[i+1][j]%=P;
        }
    }

    for(i=0;i<K;i++){
        for(j=0;j<K*2;j++)
        printf("%d,",a[i][j]);
        printf("\n");
    }

}
//exit(1);

 for(k=0;k<K;k++){   
    d=a[k][k];
    if(d<0)
    d+=P;
    for(j=0;j<K*2;j++){
        a[k][j]=mlt(a[k][j],iinv(d));
    }

 
    for(i=0;i<k;i++){
        d=a[i][k];
        if(d<0)
        d+=P;
        printf("d=%d\n",d);
        for(j=0;j<K*2;j++){
        a[i][j]-=mlt(a[k][j],d);
        if(a[i][j]<0)
        a[i][j]+=P;
        if(a[i][j]>P)
        a[i][j]%=P;
        }
    }

    for(i=0;i<K;i++){
        for(j=0;j<K*2;j++)
        printf("%d,",a[i][j]);
        printf("\n");
    }
for(i=0;i<K;i++){
    for(j=K;j<K*2;j++){
    zan[i][j-K]=a[i][j];
    V.x[i][j-K]=a[i][j];
    }
}
}

unsigned short zzz[K][K]={0};
//乾山
for(i=0;i<K;i++){
    for(j=0;j<K;j++){
        for(k=0;k<K;k++){
        zzz[i][j]+=mlt(zan[i][k],ken[k][j]);
        if(zzz[i][j]<0)
        zzz[i][j]+=P;
        if(zzz[i][j]>=P)
        zzz[i][j]%=P;
        }
     }

}
for(i=0;i<K;i++){
    for(j=0;j<K;j++){
    printf("%d,",zzz[i][j]);
    }
    printf("\n");
}
//exit(1);
return V;
}


//エラーロケーター
int gaus(MTX z)
{
    //  int a[EN][EN + 1] = { 0 };
    int p, d;
    int i, j, k;

    int a[3][4]={ 0
        /*
    {1],5],7],2]},
    {5],7],2],3]},
    {7],2],3],6]}
    */
        };
    int l=0;
    for(i=0;i<T;i++){
        for(j=0;j<T+1;j++){
        a[i][j]=z.x[i][j];
        printf("%d,",a[i][j]);
        }
        printf("\n");
    }
//exit(1);

 for(k=0;k<T;k++){   
    d=a[k][k];
    for(j=0;j<T+1;j++){
        a[k][j]=mlt(a[k][j],inv(d,P));
    }
 
    for(i=k;i<T-1;i++){
        d=a[i+1][k];
        printf("d=%d\n",d);
        for(j=0;j<T+1;j++)
        a[i+1][j]^=mlt(a[k][j],d);
    }

    for(i=0;i<T;i++){
        for(j=0;j<T+1;j++)
        printf("%d,",a[i][j]);
        printf("\n");
    }

}

    //exit(1);
    unsigned short t[T]={0};

    for (i = 0; i < T; i++)
    {
        printf("x%d = %d\n", i + 1, a[i][T]);
        t[i]=a[i][T];
    }

//exit(1);
    kakeru(z, t);


    return 0;
}

void van()
{
    int i, j;

    printf("van der\n");

    for (i = 0; i < N; i++)
        vb[0][i] = 1;
    //#pragma omp parallel for private(i, j)
    for (i = 1; i < K + 1; i++)
    {
        for (j = 0; j < N; j++)
        {
            // vb[i][j] = mltn(i, j])];
            vb[i][j] = mltn(i, j);
            // printf("%d,", vb[i][j]);
        }
        // printf("\n");
    }
}

void ban(int kk)
{
    int i, j, k;

    printf("van der\n");

    for (i = 0; i < N+1; i++)
        vb[0][i] = 1;
    //#pragma omp parallel for private(i, j)
    for (i = 1; i < kk+1; i++)
    {
        for (j = 0; j < 16; j++)
        {
            vb[i][j] = mlt2(i, j);
            // vb[i][j] = mltn(i, j])];
        }
        printf("\n");
    }
    for (i = 1; i < kk; i++)
    {
        for (j = 1; j < 16; j++)
            printf("%d,", vb[i][j]);
        printf("\n");
    }
}

int main()
{
    int i, j;
    OP L={0};
    unsigned short q[4]={1,0,1,3};
    unsigned short q2[3][4]={{4,12,7,8},
                            {12,7,8,11},
                            {7,8,11,13}
                            };
int kb[16][16]={
{64,0,0,0,102,0,150,1,0,147,0,12,0,200,54,0},
{17,17,17,17,189,189,231,150,150,100,100,18,18,142,112,112},
{68,136,13,26,190,97,183,170,73,125,250,44,88,250,54,108},
{71,201,70,202,132,145,213,225,62,140,137,220,121,46,42,126},
{123,241,227,171,220,87,226,144,122,149,110,188,202,103,103,129},
{97,248,63,195,134,164,34,19,95,117,188,9,45,155,157,211},
{77,179,141,9,66,145,30,6,20,112,61,177,129,220,159,101},
{182,37,251,219,201,69,11,78,247,172,99,10,54,154,6,18},
{117,143,12,96,50,141,120,193,70,208,206,139,44,207,12,96},
{13,101,106,29,231,140,81,160,201,84,206,195,149,124,224,179},
{117,101,197,241,49,247,124,25,250,176,148,84,50,91,110,139},
{41,46,31,217,61,178,83,66,236,58,131,8,88,193,177,47},
{19,212,153,226,54,117,227,253,104,22,232,8,96,103,154,246},
{32,189,242,222,12,92,160,253,149,224,20,226,14,34,188,255},
{55,23,202,152,198,208,25,213,34,182,74,152,195,221,164,182},
{211,213,247,20,188,154,71,212,248,37,206,230,235,56,43,148}
};
    //ban(T*2+1);
    //exit(1);
    
    MTX m={0};
    for(i=0;i<16;i++){
        for(j=0;j<16;j++)
        m.x[i][j]=kb[i][j];
    } 
    //g2();              
    gyaku(m);
    exit(1);

    L=setpol(q,4);
    printpol(o2v(L));
    printf("\n");
    //printf("%d %d %d %d\n",7],12],7],12]);
    //exit(1);

MTX z={0};
unsigned short x[T][T*2]={
    {1,2,3,1,0,0},
    {4,5,6,0,1,0},
    {7,7,7,0,0,1}
};
for(i=0;i<T;i++){
    for(j=0;j<T*2;j++)
    z.x[i][j]=x[i][j];
}
gyaku(z);
//exit(1);

  tas();
/*
  ban(8);
  for(i=1;i<5;i++){
    printf("^%d,",plus(vb[i][3],vb[i][4]));
  }
//  exit(1);
  printf("\n");
  printf("%d %d\n",inv(9),mlt(9],inv(9))]);
  gaus();
    printf("%d\n",hiku(11,121));
  //rec(s);
  exit(1);
*/


//    z=matinv(X,2);
  //  gaus(z);
    //exit(1);
unsigned short zz[N]={0};
zz[2]=1;
zz[3]=1;
zz[4]=1;
for(i=0;i<N;i++){
    if(zz[i]>0)
    printf("%d,",i);
}
//exit(1);

    ban(2*T);

    for (j = 1; j < K + 1; j++)
    {
        for (i = 1; i < N; i++)
        {
            printf("%d,", vb[j][i]);
        }
        printf("\n");
    }

    synd(zz);
    for (i = 0; i < K + 1; i++)
        printf("s%d,", syn[i]);
    printf("\n");

    for (i = 0; i < T; i++)
    {
        int count=0;
        for (j = i; j < T + i + 1; j++)
        {
            a[i][count++] = syn[j ];

            printf("%d,", a[i][j ]);
        }
        printf("\n");
    }
    for(i=0;i<T;i++){
        for(j=0;j<T+1;j++)
        printf("%d,",a[i][j]);
        printf("\n");
    }
//    exit(1);
    
    for (i = 0; i < T; i++)
    {
        int count=0;
        for (j = 0; j < T+1 ; j++)
        {
            z.x[i][j]=vb[i][j]; //a[i][j];
            printf("*%d,", z.x[i][j]);
        }
        printf("\n");
    }
    MTX tt={0};
    //printf("%d\n",inv2(z));
    //exit(1);
    //exit(1);
    printf("locater------------------------\n");
    gaus(z);
    printf("逆行列------------------------\n");
 for(i=0;i<T;i++){
    for(j=0;j<T;j++)
    z.x[i][j]=vb[i][j];
 }
    for(i=T;i<T*2;i++){
    z.x[i][i]=1;
}
    gyaku(z);
    printf("解------------------------\n");
    for(i=0;i<T;i++){
        for(j=0;j<T+1;j++)
        tt.x[i][j]=vb[i+2][j+3];
    }
    renritu(tt);
    


    return 0;
}