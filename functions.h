#include <fitsio.h>

/* Alard.c */
void        getKernelVec();
int         fillStamp(stamp_struct *, float *, float *);
double      *kernel_vector(int, int, int, int, int *);
double      *kernel_vector_PCA(int, int, int, int, int *);
void        xy_conv_stamp(stamp_struct *, float *, int, int);
void        xy_conv_stamp_PCA(stamp_struct *, float *, int, int);
void        fitKernel(stamp_struct *, float *, float *, float *, double *, double *, double *, int *);
void        build_matrix0(stamp_struct *);
void        build_scprod0(stamp_struct *, float *);
double      check_stamps(stamp_struct *, int, float *, float *);
void        build_matrix(stamp_struct *, int, double **);
void        build_scprod(stamp_struct *, int, float *, double *);
void        getStampSig(stamp_struct *, double *, float *, double *, double *, double *);
void        getFinalStampSig(stamp_struct *, float *, float *, double *);
char        check_again(stamp_struct *, double *, float *, float *, float *, double *, double *, int *);
void        spatial_convolve(float *, float **, int, int, double *, float *, int *);
double      make_kernel(int, int, double *);
double      get_background(int, int, double *);
void        make_model(stamp_struct *, double *, float *);
int         ludcmp(double **, int, int *, double *);
void        lubksb(double **, int, int *, double *);

/* Functions.c */
int         allocateStamps(stamp_struct *, int);
void        buildStamps(int, int, int, int, int *, int *, int, int, int,
                        stamp_struct *, stamp_struct *, float *, float *,
                        float, float);
void        cutStamp(float *, float *, int, int, int, int, int, stamp_struct *);
void        buildSigMask(stamp_struct *, int, int, int *);
int         cutSStamp(stamp_struct *, float *);
double      checkPsfCenter(float *, int, int, int, int, int, int, double, float, float,
			   int, int, int, int);
int         getPsfCenters(stamp_struct *, float *, int, int, double, int, int);
int         getPsfCentersORIG(stamp_struct *, float *, int, int, double, int, int);
int         getStampStats3(float *, int, int, int, int, double *, double *, double *, double *, double *, double *, double *, int, int, int);
void        getNoiseStats3(float *, float *, double *, int *, int, int);
int         stampStats(double *, int *, long, double *, double *, double *, double *, double *, double *, double *);
int         sigma_clip(float *, int, double *, double *, int);
float      *calculateAvgNoise(float *, int *, int, int, int, int, int); /* not used? */
void        freeStampMem(stamp_struct *, int);
/*int         makeNoiseImage2(float **, float, float, float *, float, float, int, int, double *);*/
/*int         makeNoiseImage3(float *, float, float, float *, float, float, int, int);*/
float       *makeNoiseImage4(float *, float, float);
void        getKernelInfo(char *);
void        readKernel(char *, int, double **, double **, int *, int *, int *, int *, double *, double *, double *, double *, double *, int *);
void        fits_get_kernel_btbl(fitsfile *, double **, int);
void        spreadMask(int *, int);
void        makeInputMask(float *, float *, int *);
void        makeOutputMask(float *, float, float, float *, float, float, int *, int *, int *);
int         hp_fits_copy_header(fitsfile *, fitsfile *, int *);
void        hp_fits_correct_data(float *, int, float, float, int);
void        hp_fits_correct_data_int(int *, int, float, float, int);
int         hp_fits_write_subset(fitsfile *, long, long, long *,
                                 float *, int *,
                                 int, float, float,
                                 int, int, int, int, int, int);
int         hp_fits_write_subset_int(fitsfile *, long, long, long *,
				     int *, int *,
				     int, float, float,
				     int, int, int, int, int, int);
void        fset(float *, double, int, int);
void        dfset(double *, double, int, int);
void        printError(int);
double      ran1(int *);
void        quick_sort(double *,int *, int);
int         imin(int, int);
int         imax(int, int);

/* armin */
void savexy(stamp_struct *, int, long, long, int);
void loadxyfile(char *,int);

/* mysinc.c */
int    swarp_remap(float *, float *, double, double, int, int,
		   int, float *, float *, int);
double luptonD(int, double);
double luptonD_appx(int, double);
void   make_lupton_kernel(double, double *, int);
void   lanczos3(double, double *);
void   lanczos4(double, double *);
void   lanczos(double, double *, int);



/* Vargs.c */
void        vargs(int, char *[]);

/* jtwarp.c */
/*
void         jtrebin(int, int, float *, int, int, int, int,
		     double *, double *, double *, float *, float *,
		     float, float, int, float, float, float, float, float);
int          jtdotri(int, int, double *, double [], double [], double, float *, int []);
void         jtsprinkle(int, int, double *, double, double, double, double,
			double, float, float *, int []);
*/
