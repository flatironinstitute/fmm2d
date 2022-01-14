
#include <complex.h>

void h2dformmpc_(int *nd, double complex *zk, double *rscale,  double *source,
                 int *ns, double complex *charge, double *center,
                 int *nterms, double complex *mpole);


void h2dmpmp_(int *nd, double complex *zk, double *rscale1, double *center1,
              double complex *hexp1, int *nterms1, double *rscale2,
              double *center2,  double complex *hexp2,  int *nterms2);
  

void h2dmpmp_(int *nd, double _Complex *zk, double *rscale0, double *center0,
                  double _Complex *mpole0, int *nterms0,
                  double *rscale1, double *center1,
                  double _Complex *mpole1, int *nterms1);


