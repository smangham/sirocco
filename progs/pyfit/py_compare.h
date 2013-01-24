

#define LINELENGTH 200
#define NCOMMENTS   10
#define NWAVES      14000
#define MAX_MOD    10
#define BIG        1.e32
#define MAXTOT     100000
#define NREG	10
#define MAXVAR	10

#include <stdio.h>
#include <string.h>

#define MAXWORD 80



typedef struct Spec
{
  char name[LINELENGTH];
  float ang;
  float w[NWAVES];
  float f[NWAVES];
  int good[NWAVES];
  int nw, ng;
} spec_dummy, *SpecPtr;


