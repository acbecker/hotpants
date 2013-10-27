typedef struct
{
   int       x0,y0;       /* origin of stamp in region coords*/
   int       x,y;         /* center of stamp in region coords*/
   int       nx,ny;       /* size of stamp */
   int       *xss;        /* x location of test substamp centers */
   int       *yss;        /* y location of test substamp centers */
   int       nss;         /* number of detected substamps, 1 .. nss     */
   int       sscnt;       /* represents which nss to use,  0 .. nss - 1 */
   double    **vectors;   /* contains convolved image data */
   double    *krefArea;   /* contains kernel substamp data */
   double    **mat;       /* fitting matrices */
   double    *scprod;     /* kernel sum solution */
   double    sum;         /* sum of fabs, for sigma use in check_stamps */
   double    mean;
   double    median;
   double    mode;        /* sky estimate */
   double    sd;
   double    fwhm;
   double    lfwhm;
   double    chi2;        /* residual in kernel fitting */
   double    norm;        /* kernel sum */
   double    diff;        /* (norm - mean_ksum) * sqrt(sum) */
} stamp_struct;

/* GLOBAL VARS POSSIBLY SET ON COMMAND LINE */
char      *template, *image, *outim;

float     tUThresh, tUKThresh, tLThresh, tGain, tRdnoise, iUThresh, iUKThresh, iLThresh, iGain, iRdnoise;
char      *tNoiseIm, *iNoiseIm, *tMaskIm, *iMaskIm, *kernelImIn, *kernelImOut, *outMask;
float     tPedestal, iPedestal;
int       hwKernel;
float     kerFitThresh, scaleFitThresh, minFracGoodStamps;
float     kfSpreadMask1, kfSpreadMask2;
int       gdXmin, gdXmax, gdYmin, gdYmax;
int       nRegX, nRegY;
char      *regFile;
char      *regKeyWord;
int       numRegKeyWord;
int       nStampY, nStampX, useFullSS;
int       nKSStamps, hwKSStamp;
char      *sstampFile;
int       findSSC;
int       kerOrder, bgOrder;
float     statSig, kerSigReject, kerFracMask;
char      *forceConvolve, *photNormalize, *figMerit;
int       sameConv, rescaleOK;
float     fillVal, fillValNoise;
char      *effFile, *noiseImage, *sigmaImage, *convImage;
int       doSum, inclNoiseImage, inclSigmaImage, inclConvImage, noClobber;
int       doKerInfo, outShort, outNShort;
float     outBzero, outBscale, outNiBzero, outNiBscale;
int       convolveVariance;
int       usePCA, fwKernelPCA;
float     **PCA;

/* GLOBAL VARS NOT SET ON COMMAND LINE */
int       ngauss, *deg_fixe;
float     *sigma_gauss;

int       rPixX, rPixY;
int       nStamps, nS, nCompKer, nC;

int       nComp, nCompBG, nBGVectors, nCompTotal;

int       fwKernel, fwStamp, hwStamp, fwKSStamp, kcStep, *indx;
int       cmpFile;
float     *temp, *temp2;
double    *check_stack,*filter_x,*filter_y,**kernel_vec;
double    **wxy,*kernel_coeffs,*kernel,**check_mat,*check_vec;
char      version[32];

/* REGION SIZED */
int       *mRData;   /* bad input data mask */

/* armin */
/* a dummy varialbe to do some testing */
int        dummy;
/* verbose for debugging */
int        verbose;
/* cmp file stuff */
char       xyfilename[1000];
int        savexyflag;
float      *xcmp,*ycmp;
int        Ncmp;
