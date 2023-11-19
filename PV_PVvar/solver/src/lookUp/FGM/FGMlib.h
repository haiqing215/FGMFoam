/* Header file for FGMlib */

# define VAR_NAME_LENGTH 16

typedef struct {
  int Ncv;
  int *Ngrid;
  double *gridpower;
  int Nvar;
  char (*varname)[VAR_NAME_LENGTH];
  double *data;
} FGM;

FGM * readFGM(const char filename[]);

int freeFGM(FGM *fgm);

// Private

int locateKeyWord(FILE *fid, char keyword[]);

int lookupFGM_1D (FGM *fgm, double *x, double *f);

int lookupFGM_2D (FGM *fgm, double *x, double *f);

int lookupFGM_ND (FGM *fgm, double *x, double *f);

int lookupFGM_NDIS (FGM *fgm, double *x, double *f);

int NDinterp (FGM *fgm, double *x, double *f);

