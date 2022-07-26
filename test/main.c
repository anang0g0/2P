#include "global.h"
#include "struct.h"

#include "debug.c"
#include "zech.c"
#include "util_mat.c"
#include "util_poly.c"
#include "mtx.c"



MTX MSM(MAT a,int k,int n){
  MTX b;

  b=mtx_new(k,n);
  b=mtx_copy(a.x,k,n);

  return b;
}

MAT SMS(MTX a){
  MAT V;
  
  memcpy(V.x,a.x,sizeof(V.x));

  return V;
}


//言わずもがな
int
main (void)
{
  int i, j, k, l;
  int count = 0;
  FILE *fp, *fq;
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
  //  {0};


  int flg, o1 = 0;
  OP f = { 0 }, r = {
    0
  }, w = {
    0
  }, ff = {
    0
  }, tt = {
    0
  };
  EX hh = { 0 };
  vec v;
  //  unsigned short dd[K*N] = {0},gg[K+1]={0};
  time_t t;
  OP r1 = { 0 }, r2 = {
    0
  };
  OP g1 = { 0 }, tmp = {
    0
  };
  int fail = 0;

  MAT G={0},GG={0},PP={0},SS={0};
  unsigned short P[N][N]=
    {
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0},
     {  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0},
     {  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
     {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
     {  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}
    };

unsigned short invP[N][N]=
  {
   {  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0},
   {  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
   {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0},
   {  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}
  };
  
  

unsigned short S[K][K]=
  {
   {1,0,0,0,0,0,1,1,},
   {1,1,0,0,1,0,0,0,},
   {1,0,1,1,1,0,0,1,},
   {1,1,0,1,0,0,0,1,},
   {0,0,0,1,0,1,1,1,},
   {0,1,1,0,0,0,1,1,},
   {1,1,1,0,0,1,1,0,},
   {0,0,0,0,0,0,0,1,},
  };
 
unsigned short inv_S[K][K]=
  {
   {0,1,1,0,1,0,1,0,},
   {1,1,1,1,0,1,0,0,},
   {0,0,0,1,1,0,1,0,},
   {1,0,0,0,1,1,1,1,},
   {1,1,0,1,1,1,1,0,},
   {0,1,1,0,1,1,0,1,},
   {1,1,1,0,1,0,1,1,},
   {0,0,0,0,0,0,0,1,},
  };


/*
   unsigned char S[K][K]=
      {
       {1,1,1,0,1,0,0,0,1,0,1,1},
       {1,0,0,0,0,0,0,1,0,0,0,1},
       {0,1,0,0,0,0,0,1,1,0,0,1},
       {0,1,1,1,0,1,1,0,1,1,1,0},
       {1,0,0,0,1,1,0,1,0,1,1,0},
       {1,1,1,1,0,0,1,0,0,0,1,1},
       {1,0,1,0,1,1,1,0,1,0,0,0},
       {1,1,1,1,0,1,1,0,1,0,0,1},
       {1,0,1,1,0,0,1,0,0,0,0,1},
       {0,0,1,1,1,1,1,0,1,0,1,1},
       {1,0,0,0,1,1,0,1,0,0,1,1},
       {0,0,0,1,0,0,0,1,1,0,1,1}
      };
unsigned char inv_S[K][K]=
  {
   {0,0,1,1,1,1,0,1,0,1,0,0},
   {0,0,1,0,0,0,0,0,1,1,1,0},
   {1,1,0,0,0,0,1,1,0,0,0,1},
   {0,0,1,0,0,1,0,0,1,0,0,1},
   {0,0,1,1,1,0,0,0,1,0,0,0},
   {0,1,0,1,1,1,0,0,1,1,0,0},
   {1,1,1,0,0,1,0,1,1,0,1,1},
   {0,1,0,0,0,0,1,0,1,0,1,1},
   {0,1,1,1,1,1,0,1,1,0,1,0},
   {0,0,1,1,0,1,1,1,1,1,0,1},
   {0,0,1,0,0,1,0,0,0,1,1,0},
   {0,0,1,1,1,1,1,1,1,1,1,1}
  };
*/
  //N,K
  MAT gen={0};
  MAT mat2={0},mat3={0};
  unsigned short q=10,p;

  
  //printf("%d\n",v2i(i2v(65535)));
  //exit(1);

  //公開鍵matはグローバル変数でスタック領域に取る。
//ヒープ領域は使うときはここを有効にしてください。
/*
  mat = malloc (N * sizeof (unsigned short *));
  base=malloc(sizeof(unsigned short)*K*N);
  #pragma omp parallel for
  for(i=0;i<N;i++){
  //連続したメモリ空間にデータを配置
  mat[i]=base+i*K;
  memset(mat[i],0,K);
  }
*/

#ifdef SRAND
  srand (SRAND);
#else
  const unsigned int seed = clock () + time (&t);
  printf ("srand(%u)\n", seed);
  srand (seed);
#endif

label:

  //makeS();
  mtx_test();
  mtx_test2();
  wait();
    
  do
    {
      fail = 0;
      j=0;
      k = 0;
      flg = 0;
      //memset (g, 0, sizeof (g));
      memset (ta, 0, sizeof (ta));
      //ginit ();

      for (i = 0; i < K + 1; i++)
        {
          if (g[K - 1] > 0)
            flg = 1;
          if (i % 2 == 1 && g[i] > 0 && i < K)
            k++;
        }
      if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        {
          w = setpol (g, K + 1);
	  j=1;
        }

      //w = setpol (g, K + 1);
      //oprintpol (w);

      //多項式の値が0でないことを確認
      for (i = 0; i < D; i++)
        {
          ta[i] = trace (w, i);
          if (ta[i] == 0)
            {
              printf ("trace 0 @ %d\n", i);
              fail = 1;
              break;
            }
        }

    }
  while (fail || j==0);
  
  oprintpol (w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");
  wait();
  
  
#pragma omp parallel for
  for (i = 0; i < N; i++)
    tr[i] = oinv (ta[i]);


  memset (mat, 0, sizeof (mat));

  printf ("\nすげ、オレもうイキそ・・・\n");
  //keygen(g);



  //パリティチェックを生成する。
  //パリティチェックに0の列があったら、なくなるまでやり直す。
  
  do
    {
      i = deta (g);
    }
  while (i < 0);
  
  //unsigned short gen[N][K]={0};

lab:

  matmul ();

  //makeS();
  //exit(1);
  MTX Q;
      G=genSGP();
      printf("after SG\n");
      pMAT(G,K,N,0);
      wait();
      Q=MSM(G,K,N);
      printf("test\n");
      mtx_print("MTX=",Q);
      wait();
      
      //置換の確認
      for(i=0;i<K;i++){
	for(j=0;j<N;j++){
	  for(k=0;k<N;k++)
	    gen.x[j][i]^=gf[mlt(fg[G.y[k][i]],fg[invP[k][j]])];
	}
      }
      printf("after invP\n");
      pMAT(gen,N,K,0);
      wait();

      for(j=0;j<N;j++){
	for(i=0;i<K;i++){
	  for(k=0;k<K;k++)
	    mat3.x[j][i]^=gf[mlt(fg[inv_S[i][k]],fg[gen.x[j][k]])];
	}
      }
      printf("after inv_S\n");
      pMAT(mat3,N,K,0);
      wait();
      
      printf("decode of G\n");
      pMAT(gen,N,K,0);
      //exit(1);
      
      
      printf("original\n");
      pMAT(G,K,N,0);
      //exit(1);
      
      
      

      printf("gen2mat\n");
      for(i=0;i<8;i++){
	for(j=0;j<M;j++)
	  printf("%2d,",mat[j][i]);
	printf("\n");
      }
      //exit(1);
      
      
      
      
      
      unsigned short code[N]={0},code2[N]={0},code3[N]={0};
      
      
      //decode開始
      k = 0;
      while (1)
	{
	  memset(code,0,sizeof(code));
	  memset(code2,0,sizeof(code2));
	  o1 = 0;
	  
	  count = 0;
	  
	  memset (zz, 0, sizeof (zz));
	  
	  //for(i=0;i<4;i++)
	  //zz[i]=i+1;
	  
	  
	  j = 0;
	  k = 0;
	  while (j < T)
	    {
	      l = xor128 () % D;
	      //printf("l=%d\n",l);
	      if (0 == zz[l])
		{
		  k=rand()%M;
		  if(k>0){
		    zz[l] = k;
		    j++;
		  }
		}
	    }
	  
	  for(i=0;i<N;i++)
	    printf("zz=%d %d\n",i,zz[i]);
	  wait();

	  for(i=0;i<N;i++)
	    code[i]=G.y[i][0];
	  for(i=0;i<N;i++)
	    code[i]^=zz[i];
	  unsigned short z3[N]={0};
	  memcpy(z3,zz,sizeof(z3));
	  for(i=0;i<N;i++){
	    for(j=0;j<N;j++)
	      code2[i]^=gf[mlt(fg[code[j]],fg[invP[j][i]])];
	  }

	  for (i = 0; i < D; i++)
	    {
	      if (code2[i] > 0)
		printf ("l=%d %d\n", i, code2[i]);
	    }
	  wait();
	  //exit(1);
	  
	  printf("code0\n");
	  for(i=0;i<N;i++)
	    printf("%d,",G.y[i][0]);
	  printf("\n");
	  printf("code\n");
	  for(i=0;i<N;i++)
	    printf("%d,",code[i]);
	  printf("\n");
	  printf("code2\n");
	  for(i=0;i<N;i++)
	    printf("%d,",code2[i]);
	  printf("\n");
	  //  exit(1);
	  
	  wait();

	  
	  vec ef={0},gh={0};
	          //code2
	  f = synd (code2);
	  //exit(1);
	  
	  f=conv(f);
	  printpol(o2v(f));
	  printf("\n");

	  /*
	  //exit(1);
	  ef=o2v(f);
	  for(i=0;i<K;i++){
	  for(k=0;k<K;k++)
	  gh.x[i]^=gf[mlt(fg[inv_S[k][i]],fg[ef.x[k]])];
	  }
	  f=v2o(gh);
	  f=conv(f);
	  //exit(1);
	  */

	  count = 0;
	  /*
	    count = 0;
	    for (i = 0; i < N; i++)
	    {
	    if (zz[i] > 0)
	    count++;
	    }
	    printf ("%d\n", count);
	  */
	  printpol(o2v(w));
	  printf(" ==========goppa\n");
	  printpol(o2v(f));
	  printf(" ==========synd\n");

      //exit(1);
	  
      r = decode (w, f);

      for(i=0;i<N;i++)
	printf("%d,",zz[i]);
      printf("\n");
      wait();
      for (i = 0; i < T; i++)
        {
	  if(i==0)
	    printf("i==0 %d\n",r.t[i].a);
	    
          if (r.t[i].a > 0 && count > 0)        // == r.t[i].n)
            {
              printf ("e=%d %d %d %s\n", i,r.t[i].a, r.t[i].n, "お");
              count++;
            }
          if (count == 0 && r.t[i].a > 0)
            {
              printf ("\ne=%d %d %d %s\n", i,r.t[i].a, r.t[i].n, "う");
              count++;
            }

        }
      if (count != T)
        {
          printf ("error pattarn too few %d\n", count);
	  for(i=0;i<N;i++)
	    printf("%d,",zz[i]);
	  printf("\n");
	  
	  unsigned short e[N]={0},mm[N]={0};
	  for(i=0;i<N;i++){
	    for(j=0;j<N;j++)
	      e[i]^=gf[mlt(fg[zz[j]],fg[invP[j][i]])];
	    printf("%d,",e[i]);
	  }
	  printf("\n");
          exit (1);
        }
      
      printf ("err=%dっ！！\n", count);
      if(count<T)
	{
	  
	  printf("%d baka1\n",count);
	  exit(1);
	}


      ef=o2v(r);
      for(i=0;i<N;i++){
	printf("e=%d\n",ef.x[i]);
	//code2[i]^=ef.x[i];
      }
      printf("\n\n");
      printf("code2\n");
      for(i=0;i<N;i++)
	printf("%d,",code2[i]);
      printf("\n");

      unsigned short tt[K]={0},t3[K]={0};
      for(i=0;i<K;i++)
	tt[i]=code2[i]^ef.x[i];
      printf("m=");
      for(i=0;i<K;i++){
	for(j=0;j<K;j++)
	  t3[i]^=gf[mlt(fg[tt[j]],fg[inv_S[j][i]])];
	printf("%d",t3[i]);
      }
      printf("\n");
      //exit(1);
      wait();
      
      //exit(1);
      //goto label;

    
    patta:

      
      
      //printf("パターソンアルゴリズムを実行します。何か数字を入れてください。\n");
      //wait();
      count=0;
      memset (z1, 0, sizeof (z1));
      memset(code,0, sizeof (code));
      for(i=0;i<N;i++)
	code[i]=G.y[i][2];
      j = 0;
      
      while (j < T * 2)
        {
          l = xor128 () % D;
          //printf ("l=%d\n", l);
          if (0 == z1[l])
            {
              z1[l] = 1;
	      printf("l=%d\n",l);
              j++;
            }
        }
      //wait();
      unsigned short code2[N]={0};
      printf("code\n");
      for(i=0;i<N;i++){
	code[i]^=z1[i];
	printf("%d,",code[i]);
      }
      printf("\n");
      printf("code2\n");
      for(i=0;i<N;i++){
	for(j=0;j<N;j++)
	  code2[i]^=gf[mlt(fg[code[j]],fg[invP[j][i]])];
	printf("%d,",code2[i]);
      }
      printf("\n");
      //for(i=0;i<8;i++)
      //z1[i]=1;

      //encryotion
      //test (w, z1);
      //memcpy(mat,mat2,sizeof(mat));

      f = synd (code2);
      /*
      memset(gh.x,0,sizeof(gh.x));
      ef=o2v(f);
      for(i=0;i<K;i++){
      for(k=0;k<K;k++)
	gh.x[i]^=gf[mlt(inv_S[i][k],fg[ef.x[k]])];
      }
      f=v2o(gh);
      f=conv(f);
      */
      
      //バグトラップのためのコード（省略）
      //trap(w,f);
      //バグトラップ（ここまで）

      count = 0;
      //復号化の本体
      v = pattarson (w, f);

      
      //エラー表示
      for (i = 0; i < T * 2; i++)
        {
          if (i == 0){
            printf ("error position=%d %d う\n", i, v.x[i]);
	    count++;
	  }
          if (i > 0 && v.x[i] > 0){
            printf ("error position=%d %d お\n", i, v.x[i]);
	    count++;
	  }
          if (i > 0 && v.x[i] == 0)
            {
              printf ("baka %d %d\n", i, v.x[i]);
	      printf("v.x[K-1]=%d\n",v.x[K-1]);
	      printpol(o2v(w));
	      printf(" ============goppa\n");
	      printf("{");
	      for(i=0;i<N;i++)
		printf("%d,",z1[i]);
	      printf("};\n");
	      //AA++;
	      //wait();
	      break;
	      //
              //exit (1);
            }
        }
      
      if(count==T*2){
      printf ("err=%dっ!! \n", count);
      B++;
      unsigned short e[N]={0},mm[N]={0},x[K]={0};
      for(i=0;i<K;i++)
	e[v.x[i]]=1;
      for(i=0;i<N;i++)
	printf("%d,",code2[i]);
      printf(" ========code2\n");
      printf("m'=");
      for(i=0;i<N;i++){
	mm[i]=code2[i]^e[i];
	printf("%d,",mm[i]);
      }
      printf("\n");
      //wait();
      for(i=0;i<K;i++)
	x[i]=mm[i];
      printf("m=");
      for(i=0;i<K;i++){
	for(k=0;k<K;k++)
	  gh.x[i]^=gf[mlt(fg[inv_S[k][i]],fg[x[k]])];
	printf("%d,",gh.x[i]);
      }
      printf("\n");
		  
      }
      if (count < T * 2){
        printf ("error is too few\n");
	AA++;
	memcpy(zz,z1,sizeof(zz));
	printf("{");
	for (i = 0; i < D; i++)
	  printf ("%d,", z1[i]);
	printf ("};\n");
	printpol(o2v(w));
	printf(" =========goppa\n");
	exit(1);
      }
      if(AA==1000){
	printf("B=%d",B);
	wait();
      }
      if(B>10000){
	count=0;
	printf("false=%d\n",AA);
	printf("success=%d\n",B);
	//for(i=0;i<16;i++)
	//printf("%d,",zz[i]);
	//printf("\n");
	printf("C=%d\n",C);
	printpol(o2v(w));
	printf(" =======goppa\n");
	exit(1);
      }

      //exit(1);
      //goto lab;
      //wait();
      
      break;
      }

  return 0;
}
