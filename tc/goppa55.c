#define N 64


static int gf[64] = 
  {0, 1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20, 40, 19, 38, 15, 
   30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56, 51, 37, 9, 18, 36, 11, 22, 44, 
   27, 54, 47, 29, 58, 55, 45, 25, 50, 39, 13, 26, 52, 43, 21, 42, 23, 46, 
   31, 62, 63, 61, 57, 49, 33};

/*
static int gf[512]={
0,1,2,4,8,16,32,64,128,256,305,
339,407,31,62,124,248,496,209,418,117,
234,468,153,306,341,411,7,14,28,56,
112,224,448,177,354,501,219,438,93,186,
372,473,131,262,317,331,423,127,254,508,
201,402,21,42,84,168,336,401,19,38,
76,152,304,337,403,23,46,92,184,368,
465,147,294,381,459,167,334,429,107,214,
428,105,210,420,121,242,484,249,498,213,
426,101,202,404,25,50,100,200,400,17,
34,68,136,272,273,275,279,287,271,303,
367,495,239,478,141,282,261,315,327,447,
79,158,316,329,419,119,238,476,137,274,
277,283,263,319,335,431,111,222,444,73,
146,292,377,451,183,366,493,235,470,157,
314,325,443,71,142,284,265,291,375,479,
143,286,269,299,359,511,207,414,13,26,
52,104,208,416,113,226,452,185,370,469,
155,310,349,395,39,78,156,312,321,435,
87,174,348,393,35,70,140,280,257,307,
343,415,15,30,60,120,240,480,241,482,
245,490,229,458,165,330,421,123,246,492,
233,466,149,298,357,507,199,398,45,90,
180,360,481,243,486,253,506,197,394,37,
74,148,296,353,499,215,430,109,218,436,
89,178,356,505,195,390,61,122,244,488,
225,450,181,362,485,251,502,221,442,69,
138,276,281,259,311,351,399,47,94,188,
376,449,179,358,509,203,406,29,58,116,
232,464,145,290,373,475,135,270,301,363,
487,255,510,205,410,5,10,20,40,80,
160,320,433,83,166,332,425,99,198,396,
41,82,164,328,417,115,230,460,169,338,
405,27,54,108,216,432,81,162,324,441,
67,134,268,297,355,503,223,446,77,154,
308,345,387,55,110,220,440,65,130,260,
313,323,439,95,190,380,457,163,326,445,
75,150,300,361,483,247,494,237,474,133,
266,293,379,455,191,382,461,171,342,413,
11,22,44,88,176,352,497,211,422,125,
250,500,217,434,85,170,340,409,3,6,
12,24,48,96,192,384,49,98,196,392,
33,66,132,264,289,371,471,159,318,333,
427,103,206,412,9,18,36,72,144,288,
369,467,151,302,365,491,231,462,173,346,
389,59,118,236,472,129,258,309,347,391,
63,126,252,504,193,386,53,106,212,424,
97,194,388,57,114,228,456,161,322,437,
91,182,364,489,227,454,189,378,453,187,
374,477,139,278,285,267,295,383,463,175,
350,397,43,86,172,344,385,51,102,204,
408};
*/
/*
gf[256] = 
  {0, 1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 
   152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 
   78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70, 140, 
   5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161, 95, 190, 
   97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231, 
   211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 
   175, 67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 
   31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197, 151, 51, 102, 204, 
   133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 
   168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115, 
   230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 
   255, 227, 219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 
   65, 130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 
   166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 
   9, 18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 
   22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142};
*/
/*
static int gf[64]={
0,1,2,4,8,16,32,33,35,39,47,
63,31,62,29,58,21,42,53,11,22,
44,57,19,38,45,59,23,46,61,27,
54,13,26,52,9,18,36,41,51,7,
14,28,56,17,34,37,43,55,15,30,
60,25,50,5,10,20,40,49,3,6,
12,24,48};
*/



main(){
unsigned int g,z,k[4],g2[N],g3[N];
int i,j,count=0,l,m,n;


for(z=0;z<N;z++){

k[0]=add(mltn(8,z),mlt(35,mltn(7,z)));
k[1]=add(mlt(13,mltn(6,z)),mlt(54,mltn(5,z)));
k[2]=add(mlt(32,mltn(4,z)),mlt(10,mltn(3,z)));
k[3]=add(add(mlt(mlt(51,z),z),mlt(29,z)),23);

 k[0]=add(add(k[0],k[1]),add(k[2],k[3])); 
/* k[0]=add(mltn(41,z),31); */

/* g=gf[mltn(41,z)]^gf[31]; */

g= gf[mltn(8,z)]^gf[mlt(35,mltn(7,z))]^gf[mlt(13,mltn(6,z))]^
gf[mlt(54,mltn(5,z))]^gf[mlt(32,mltn(4,z))]^gf[mlt(10,mltn(3,z))]^
gf[mlt(mlt(51,z),z)]^gf[mlt(29,z)]^gf[23];


 printf("%d ",gf[div(1,k[0])]); 
if(g!=0){
g2[count]=k[0]; /* G */
g3[count]=z; /* alpha */
count++;
}
}

printf("\n\n");

for(i=0;i<count;i++)
 printf("%d ",gf[div(1,g2[i])]); 
printf("\n");

for(i=0;i<count;i++)
 printf("%d ",gf[div(g3[i],g2[i])]); 
printf("\n");
for(j=2;j<8;j++){
for(i=0;i<count;i++)
 printf("%d ",gf[div(mltn(j,g3[i]),g2[i])]); 
printf("\n");
}
printf("%d\n",count);

for(j=0;j<8;j++){
  for(i=0;i<count;i++){
    l=div(mltn(j,g3[i],g2[i]));
    for(n=0;n<64;n++){
      if(mlt(n,l)==0)
	printf("%d@ ",n);
    }
  } 
printf("\n");
}

}



int mlt(int x,int y)
{
    if(x==0||y==0){
        return(0);
    }
    return ((x+y-2)%(N-1))+1;
}


int mltn(int n,int x)
{
int i,j;

i=x;
for(j=1;j<n;j++)
i=mlt(i,x);

return i;

}


int div(int x,int y)
{
    if(y==0){
	/* printf("y��0\n"); */
	/* exit(1); */
    }
    else if(x==0){
        return 0;
    }
    return ((x-y+(N-1))%(N-1))+1;
  }


int add(int x,int y)
{
    int z;


    if(x==0){
        return(y);
    }
    else if(y==0){
        return(x);
    }
    else if(x > y){
        z=(gf[x-y+1]-1)+(y-1);
        z=(z%(N-1))+1;
    }
    else if(x < y){
        z=(gf[y-x+1]-1)+x-1;
        z=(z%(N-1))+1;
    }
    else{
        z=0;
    }
    return(z);
}

