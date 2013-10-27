#include<stdio.h>
#include<string.h>
#include<strings.h>
#include<stdlib.h>
#include<assert.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/*
  
  Set defaults for global variables, parse any command line arguments,
  and do some basic insanity checks.

  08/20/01 acbecker@physics.bell-labs.com
  
*/

void vargs(int argc, char *argv[]) {

   int iarg, i, status=0;
   char help[8192], *scrStr;
   fitsfile *kPtr;
   int kBitpix, kNaxis;
   long kNaxes[MAXDIM];

   /* set the defaults */
   image    = NULL;
   template = NULL;
   outim    = NULL;

   tNoiseIm    = NULL;
   iNoiseIm    = NULL;
   tMaskIm     = NULL;
   iMaskIm     = NULL;
   kernelImIn  = NULL;
   kernelImOut = NULL;
   outMask     = NULL;
   
   /* nullify arrays of pointers */
   mRData        = NULL;
   
   check_stack   = NULL;
   filter_x      = NULL;
   filter_y      = NULL;
   temp          = NULL;
   kernel_vec    = NULL;
   wxy           = NULL;
   kernel_coeffs = NULL;
   kernel        = NULL;
   check_mat     = NULL;
   check_vec     = NULL;
   
   ngauss      = D_NGAUSS;
   deg_fixe    = (int *)calloc(ngauss, sizeof(int));
   deg_fixe[0] = D_DEG_GAUSS1;
   deg_fixe[1] = D_DEG_GAUSS2;
   deg_fixe[2] = D_DEG_GAUSS3;

   sigma_gauss    = (float *)calloc(ngauss, sizeof(float));
   sigma_gauss[0] = (1.0/(2.0*D_SIG_GAUSS1*D_SIG_GAUSS1));
   sigma_gauss[1] = (1.0/(2.0*D_SIG_GAUSS2*D_SIG_GAUSS2));
   sigma_gauss[2] = (1.0/(2.0*D_SIG_GAUSS3*D_SIG_GAUSS3));


   /* Things also possibly changed on the command line */
   tUThresh  = D_UTHRESH;
   tUKThresh = tUThresh;
   tLThresh  = D_LTHRESH;
   tGain     = D_GAIN;
   tRdnoise  = D_RDNOISE;

   iUThresh  = D_UTHRESH;
   iUKThresh = iUThresh;
   iLThresh  = D_LTHRESH;
   iGain     = D_GAIN;
   iRdnoise  = D_RDNOISE;

   outBzero    = D_OUTBZERO;
   outBscale   = D_OUTBSCALE;
   outNiBzero  = D_OUTBZERO;
   outNiBscale = D_OUTBSCALE;

   tPedestal       = D_PEDESTAL;
   iPedestal       = D_PEDESTAL;
   hwKernel        = D_HWKERNEL;
   fwKernelPCA     = 2 * D_HWKERNEL + 1;
   kerFitThresh    = D_KFITTHRESH;
   minFracGoodStamps = D_NFITTHRESH;
   scaleFitThresh  = D_SFITTHRESH;
   
   kfSpreadMask1   = D_INMASKFSPREAD;
   kfSpreadMask2   = D_OUMASKFSPREAD;   
   nRegX           = D_NREGIONS;
   nRegY           = D_NREGIONS;
   nStampY         = D_NSTAMPS;
   nStampX         = D_NSTAMPS;
   useFullSS       = D_USEFULLSS;
   cmpFile         = D_CMPFILE;
   findSSC         = D_FINDSSCENTERS;
   nKSStamps       = D_NKSSTAMPS;
   hwKSStamp       = D_HWKSSTAMP;
   kerOrder        = D_KORDER;
   bgOrder         = D_BGORDER;
   statSig         = D_STATSIG;
   kerSigReject    = D_KSIGREJECT;
   kerFracMask     = D_KFRACMASK;
   rescaleOK       = D_RESCALEOK;
   photNormalize   = D_NORMALIZE;
   forceConvolve   = D_CONVOLVE;
   figMerit        = D_FIGMERIT;
   fillVal         = D_FILL;
   fillValNoise    = D_FILLNOISE;
   doSum           = D_OUTSUM;
   inclNoiseImage  = D_NOILAYER;
   inclSigmaImage  = D_SIGLAYER;
   inclConvImage   = D_CONVLAYER;
   noClobber       = D_NOCLOBBER;
   doKerInfo       = D_KINFO;
   outShort        = D_OUTSHORT;
   outNShort       = D_OUTSHORT;

   savexyflag      = D_SAVEXY;
   verbose         = D_VERBOSE;

   convolveVariance= D_CONVVAR;
   usePCA          = D_USEPCA;
   
   sprintf(help, "Usage:  hotpants [options]\n");
   
   sprintf(help, "%sVersion %s\n", help, version);

   sprintf(help, "%sRequired options:\n", help);
   sprintf(help, "%s   [-inim fitsfile]  : comparison image to be differenced\n", help);
   sprintf(help, "%s   [-tmplim fitsfile]: template image\n", help);
   sprintf(help, "%s   [-outim fitsfile] : output difference image\n\n", help);

   sprintf(help, "%sAdditional options:\n", help);
   sprintf(help, "%s   [-tu tuthresh]    : upper valid data count, template (%.0f)\n", help, D_UTHRESH);
   sprintf(help, "%s   [-tuk tucthresh]  : upper valid data count for kernel, template (tuthresh)\n", help);
   sprintf(help, "%s   [-tl tlthresh]    : lower valid data count, template (%.0f)\n", help, D_LTHRESH);
   sprintf(help, "%s   [-tg tgain]       : gain in template (%.0f)\n", help, D_GAIN);
   sprintf(help, "%s   [-tr trdnoise]    : e- readnoise in template (%.0f)\n", help, D_RDNOISE);
   sprintf(help, "%s   [-tp tpedestal]   : ADU pedestal in template (%.0f)\n", help, D_PEDESTAL);
   sprintf(help, "%s   [-tni fitsfile]   : input template noise array (undef)\n", help);
   sprintf(help, "%s   [-tmi fitsfile]   : input template mask image (undef)\n", help);

   sprintf(help, "%s   [-iu iuthresh]    : upper valid data count, image (%.0f)\n", help, D_UTHRESH);
   sprintf(help, "%s   [-iuk iucthresh]  : upper valid data count for kernel, image (iuthresh)\n", help);
   sprintf(help, "%s   [-il ilthresh]    : lower valid data count, image (%.0f)\n", help, D_LTHRESH);
   sprintf(help, "%s   [-ig igain]       : gain in image (%.0f)\n", help, D_GAIN);
   sprintf(help, "%s   [-ir irdnoise]    : e- readnoise in image (%.0f)\n", help, D_RDNOISE);
   sprintf(help, "%s   [-ip ipedestal]   : ADU pedestal in image (%.0f)\n", help, D_PEDESTAL);
   sprintf(help, "%s   [-ini fitsfile]   : input image noise array (undef)\n", help);
   sprintf(help, "%s   [-imi fitsfile]   : input image mask image (undef)\n\n", help);

   sprintf(help, "%s   [-ki fitsfile]    : use kernel table in image header (undef)\n", help);
   sprintf(help, "%s   [-r rkernel]      : convolution kernel half width (%d)\n", help, hwKernel);
   sprintf(help, "%s   [-kcs step]       : size of step for spatial convolution (2*rkernel + 1)\n", help);
   sprintf(help, "%s   [-ft fitthresh]   : RMS threshold for good centroid in kernel fit (%.1f)\n", help, kerFitThresh);
   sprintf(help, "%s   [-sft scale]      : scale fitthresh by this fraction if... (%.1f)\n", help, scaleFitThresh);
   sprintf(help, "%s   [-nft fraction]   : this fraction of stamps are not filled (%.1f)\n", help, minFracGoodStamps);
   
   sprintf(help, "%s   [-mins spread]    : Fraction of kernel half width to spread input mask (%.1f)\n", help, kfSpreadMask1);
   sprintf(help, "%s   [-mous spread]    : Ditto output mask, negative = no diffim masking (%.1f)\n", help, kfSpreadMask2);
   sprintf(help, "%s   [-omi  fitsfile]  : Output bad pixel mask (undef)\n", help);
   
   sprintf(help, "%s   [-gd xmin xmax ymin ymax]\n", help);
   sprintf(help, "%s                     : only use subsection of full image (full image)\n\n", help);
   
   sprintf(help, "%s   [-nrx xregion]    : number of image regions in x dimension (%d)\n", help, nRegX);
   sprintf(help, "%s   [-nry yregion]    : number of image regions in y dimension (%d)\n", help, nRegY);
   sprintf(help, "%s   -- OR --\n", help);
   sprintf(help, "%s   [-rf regionfile]  : ascii file with image regions 'xmin:xmax,ymin:ymax'\n", help);
   sprintf(help, "%s   -- OR --\n", help);
   sprintf(help, "%s   [-rkw keyword num]: header 'keyword[0->(num-1)]' indicates valid regions\n\n", help);

   sprintf(help, "%s   [-nsx xstamp]     : number of each region's stamps in x dimension (%d)\n", help, nStampX);
   sprintf(help, "%s   [-nsy ystamp]     : number of each region's stamps in y dimension (%d)\n", help, nStampY);
   sprintf(help, "%s   -- OR --\n", help);
   sprintf(help, "%s   [-ssf stampfile]  : ascii file indicating substamp centers 'x y'\n", help);
   sprintf(help, "%s   -- OR --\n", help);
   sprintf(help, "%s   [-cmp cmpfile]    : .cmp file indicating substamp centers 'x y'\n\n", help);

   sprintf(help, "%s   [-afssc find]     : autofind stamp centers so #=-nss when -ssf,-cmp (%d)\n", help, findSSC);
   sprintf(help, "%s   [-nss substamps]  : number of centroids to use for each stamp (%d)\n", help, nKSStamps);
   sprintf(help, "%s   [-rss radius]     : half width substamp to extract around each centroid (%d)\n\n", help, hwKSStamp);

   sprintf(help, "%s   [-savexy file]    : save positions of stamps for convolution kernel (undef)\n", help);
   
   sprintf(help, "%s   [-c  toconvolve]  : force convolution on (t)emplate or (i)mage (undef)\n", help);
   sprintf(help, "%s   [-n  normalize]   : normalize to (t)emplate, (i)mage, or (u)nconvolved (%s)\n", help, photNormalize);
   sprintf(help, "%s   [-fom figmerit]   : (v)ariance, (s)igma or (h)istogram convolution merit (%s)\n", help, figMerit);
   sprintf(help, "%s   [-sconv]          : all regions convolved in same direction (%d)\n", help, sameConv);
   
   sprintf(help, "%s   [-ko kernelorder] : spatial order of kernel variation within region (%d)\n", help, kerOrder);
   sprintf(help, "%s   [-bgo bgorder]    : spatial order of background variation within region (%d)\n", help, bgOrder);
   sprintf(help, "%s   [-ssig statsig]   : threshold for sigma clipping statistics  (%.1f)\n", help, statSig);
   sprintf(help, "%s   [-ks badkernelsig]: high sigma rejection for bad stamps in kernel fit (%.1f)\n", help, kerSigReject);
   sprintf(help, "%s   [-kfm kerfracmask]: fraction of abs(kernel) sum for ok pixel (%.3f)\n", help, kerFracMask);
   sprintf(help, "%s   [-okn]            : rescale noise for 'ok' pixels (%d)\n", help, rescaleOK);
   sprintf(help, "%s   [-fi fill]        : value for invalid (bad) pixels (%.1e)\n", help, fillVal);
   sprintf(help, "%s   [-fin fill]       : noise image only fillvalue (%.1e)\n", help, fillValNoise);
   sprintf(help, "%s   [-convvar]        : convolve variance not noise (%d)\n\n", help, convolveVariance);
   
   sprintf(help, "%s   [-oni fitsfile]   : output noise image (undef)\n", help);
   sprintf(help, "%s   [-ond fitsfile]   : output noise scaled difference image (undef)\n", help);
   sprintf(help, "%s   [-nim]            : add noise image as layer to sub image (%d)\n", help, inclNoiseImage);
   sprintf(help, "%s   [-ndm]            : add noise-scaled sub image as layer to sub image (%d)\n\n", help, inclSigmaImage);
   sprintf(help, "%s   [-oci fitsfile]   : output convolved image (undef)\n", help);
   sprintf(help, "%s   [-cim]            : add convolved image as layer to sub image (%d)\n\n", help, inclConvImage);
   
   sprintf(help, "%s   [-allm]           : output all possible image layers\n\n", help);

   sprintf(help, "%s   [-nc]             : do not clobber output image (%d)\n", help, noClobber);
   sprintf(help, "%s   [-hki]            : print extensive kernel info to output image header (%d)\n\n", help, doKerInfo);
   sprintf(help, "%s   [-oki fitsfile]   : new fitsfile with kernel info (under)\n\n", help);
   sprintf(help, "%s   [-sht]            : output images 16 bitpix int, vs -32 bitpix float (%d)\n", help, outShort);
   sprintf(help, "%s   [-obs bscale]     : if -sht, output image BSCALE, overrides -inim (%.1f)\n", help, D_OUTBSCALE);
   sprintf(help, "%s   [-obz bzero]      : if -sht, output image BZERO , overrides -inim (%.1f)\n", help, D_OUTBZERO);
   sprintf(help, "%s   [-nsht]           : output noise image 16 bitpix int, vs -32 bitpix float (%d)\n", help, outNShort);
   sprintf(help, "%s   [-nbs bscale]     : noise image only BSCALE, overrides -obs (%.1f)\n", help, D_OUTBSCALE);
   sprintf(help, "%s   [-nbz bzero]      : noise image only BZERO,  overrides -obz (%.1f)\n\n", help, D_OUTBZERO);

   sprintf(help, "%s   [-ng  ngauss degree0 sigma0 .. degreeN sigmaN]\n", help);
   sprintf(help, "%s                     : ngauss = number of gaussians which compose kernel (%d)\n", help, ngauss);
   sprintf(help, "%s                     : degree = degree of polynomial associated with gaussian #\n", help);
   sprintf(help, "%s                                (%d %d %d)\n", help, D_DEG_GAUSS1, D_DEG_GAUSS2, D_DEG_GAUSS3);
   sprintf(help, "%s                     : sigma  = width of gaussian #\n", help);
   sprintf(help, "%s                                (%.2f %.2f %.2f)\n", help, D_SIG_GAUSS1, D_SIG_GAUSS2, D_SIG_GAUSS3);
   sprintf(help, "%s                     : N = 0 .. ngauss - 1\n\n", help);
   sprintf(help, "%s                     : (%d %d %.2f %d %.2f %d %.2f\n", help, ngauss, D_DEG_GAUSS1, D_SIG_GAUSS1, D_DEG_GAUSS2, D_SIG_GAUSS2, D_DEG_GAUSS3, D_SIG_GAUSS3);

   sprintf(help, "%s   [-pca nk k0.fits ... n(k-1).fits]\n", help);
   sprintf(help, "%s                     : nk      = number of input basis functions\n", help);
   sprintf(help, "%s                     : k?.fits = name of fitsfile holding basis function\n", help);
   sprintf(help, "%s                     : Since this uses input basis functions, it will fix :\n", help);
   sprintf(help, "%s                     :    hwKernel \n", help);
   sprintf(help, "%s                     :    \n", help);
   
   sprintf(help, "%s   [-v] verbosity    : level of verbosity, 0-2 (%d)\n", help, verbose);

   sprintf(help, "%s NOTE: Fits header params will be added to the difference image\n", help);
   sprintf(help, "%s       COMMAND             (what was called on the command line)\n", help);
   sprintf(help, "%s       NREGION             (number of regions in image)\n", help);
   sprintf(help, "%s       PHOTNORM            (to which system the difference image is normalized)\n", help);
   sprintf(help, "%s       TARGET              (image which was differenced)\n", help);
   sprintf(help, "%s       TEMPLATE            (template for the difference imaging)\n", help);
   sprintf(help, "%s       DIFFIM              (output difference image)\n", help);
   sprintf(help, "%s       MASKVAL             (value for masked pixels)\n", help);
   sprintf(help, "%s       REGION\?\?            (IRAF-format limits for each region in the image)\n", help);
   sprintf(help, "%s       CONVOL\?\?            (which image was convolved for each region)\n", help);
   sprintf(help, "%s       KSUM\?\?              (sum of the convolution kernel for each region)", help);

   /* read in command options */
   if (argc < 2) {
      fprintf(stderr, "%s\n", help);
      exit(1);
   }
   
   for (iarg=1; iarg < argc; iarg++) {
      if (argv[iarg][0]=='-') {
         
         if (strcasecmp(argv[iarg]+1,"inim")==0) {
            image = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"tmplim")==0) {
            template = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"outim")==0) {
            outim = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"tu")==0) {
            sscanf(argv[++iarg], "%f", &tUThresh);
         } else if (strcasecmp(argv[iarg]+1,"tuk")==0) {
            sscanf(argv[++iarg], "%f", &tUKThresh);
         } else if (strcasecmp(argv[iarg]+1,"tl")==0) {
            sscanf(argv[++iarg], "%f", &tLThresh);
         } else if (strcasecmp(argv[iarg]+1,"tg")==0) {
            sscanf(argv[++iarg], "%f", &tGain);
         } else if (strcasecmp(argv[iarg]+1,"tr")==0) {
            sscanf(argv[++iarg], "%f", &tRdnoise);
         } else if (strcasecmp(argv[iarg]+1,"tp")==0) {
            sscanf(argv[++iarg], "%f", &tPedestal);
         } else if (strcasecmp(argv[iarg]+1,"tni")==0) {
            tNoiseIm = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"tmi")==0) {
            tMaskIm = argv[++iarg];
            
         } else if (strcasecmp(argv[iarg]+1,"iu")==0) {
            sscanf(argv[++iarg], "%f", &iUThresh);
         } else if (strcasecmp(argv[iarg]+1,"iuk")==0) {
            sscanf(argv[++iarg], "%f", &iUKThresh);
         } else if (strcasecmp(argv[iarg]+1,"il")==0) {
            sscanf(argv[++iarg], "%f", &iLThresh);
         } else if (strcasecmp(argv[iarg]+1,"ig")==0) {
            sscanf(argv[++iarg], "%f", &iGain);
         } else if (strcasecmp(argv[iarg]+1,"ir")==0) {
            sscanf(argv[++iarg], "%f", &iRdnoise);
         } else if (strcasecmp(argv[iarg]+1,"ip")==0) {
            sscanf(argv[++iarg], "%f", &iPedestal);
         } else if (strcasecmp(argv[iarg]+1,"ini")==0) {
            iNoiseIm = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"imi")==0) {
            iMaskIm = argv[++iarg];
	 } else if (strcasecmp(argv[iarg]+1,"ki")==0) {
            kernelImIn = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"r")==0) {
            sscanf(argv[++iarg], "%d", &hwKernel);
         } else if (strcasecmp(argv[iarg]+1,"kcs")==0) {
            sscanf(argv[++iarg], "%d", &kcStep);
         } else if (strcasecmp(argv[iarg]+1,"ft")==0) {
            sscanf(argv[++iarg], "%f", &kerFitThresh);
         } else if (strcasecmp(argv[iarg]+1,"sft")==0) {
            sscanf(argv[++iarg], "%f", &scaleFitThresh);
         } else if (strcasecmp(argv[iarg]+1,"nft")==0) {
            sscanf(argv[++iarg], "%f", &minFracGoodStamps);
         } else if (strcasecmp(argv[iarg]+1,"mins")==0) {
            sscanf(argv[++iarg], "%f", &kfSpreadMask1);
         } else if (strcasecmp(argv[iarg]+1,"mous")==0) {
            sscanf(argv[++iarg], "%f", &kfSpreadMask2);
         } else if (strcasecmp(argv[iarg]+1,"omi")==0) {
            outMask = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"gd")==0) {
            sscanf(argv[++iarg], "%d", &gdXmin);
            sscanf(argv[++iarg], "%d", &gdXmax);
            sscanf(argv[++iarg], "%d", &gdYmin);
            sscanf(argv[++iarg], "%d", &gdYmax);
         } else if (strcasecmp(argv[iarg]+1,"nrx")==0) {
            sscanf(argv[++iarg], "%d", &nRegX);
         } else if (strcasecmp(argv[iarg]+1,"nry")==0) {
            sscanf(argv[++iarg], "%d", &nRegY);
         } else if (strcasecmp(argv[iarg]+1,"rf")==0) {
            regFile = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"rkw")==0) {
            regKeyWord = argv[++iarg];
            sscanf(argv[++iarg], "%d", &numRegKeyWord);
         } else if (strcasecmp(argv[iarg]+1,"nsx")==0) {
            sscanf(argv[++iarg], "%d", &nStampX);
         } else if (strcasecmp(argv[iarg]+1,"nsy")==0) {
            sscanf(argv[++iarg], "%d", &nStampY);
         } else if (strcasecmp(argv[iarg]+1,"ssf")==0) {
            sstampFile = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"cmp")==0) {
            sstampFile = argv[++iarg];
	    cmpFile = 1;
         } else if (strcasecmp(argv[iarg]+1,"afssc")==0) {
            sscanf(argv[++iarg], "%d", &findSSC);
         } else if (strcasecmp(argv[iarg]+1,"nss")==0) {
            sscanf(argv[++iarg], "%d", &nKSStamps);
         } else if (strcasecmp(argv[iarg]+1,"rss")==0) {
            sscanf(argv[++iarg], "%d", &hwKSStamp);
         } else if (strcasecmp(argv[iarg]+1,"ko")==0) {
            sscanf(argv[++iarg], "%d", &kerOrder);
         } else if (strcasecmp(argv[iarg]+1,"bgo")==0) {
            sscanf(argv[++iarg], "%d", &bgOrder);
         } else if (strcasecmp(argv[iarg]+1,"ssig")==0) {
            sscanf(argv[++iarg], "%f", &statSig);
         } else if (strcasecmp(argv[iarg]+1,"ks")==0) {
            sscanf(argv[++iarg], "%f", &kerSigReject);
         } else if (strcasecmp(argv[iarg]+1,"kfm")==0) {
            sscanf(argv[++iarg], "%f", &kerFracMask);
         } else if (strcasecmp(argv[iarg]+1,"okn")==0) {
            rescaleOK = 1;
         } else if (strcasecmp(argv[iarg]+1,"c")==0) {
            forceConvolve = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"n")==0) {
            photNormalize = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"fom")==0) {
            figMerit = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"sconv")==0) {
            sameConv = 1;
         } else if (strcasecmp(argv[iarg]+1,"fi")==0) {
            sscanf(argv[++iarg], "%f", &fillVal);
         } else if (strcasecmp(argv[iarg]+1,"ef")==0) {
            effFile = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"sum")==0) {
            doSum = 1;
         } else if (strcasecmp(argv[iarg]+1,"oni")==0) {
            noiseImage = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"ond")==0) {
            sigmaImage = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"nim")==0) {
            inclNoiseImage = 1;
         } else if (strcasecmp(argv[iarg]+1,"ndm")==0) {
            inclSigmaImage = 1;
         } else if (strcasecmp(argv[iarg]+1,"oci")==0) {
            convImage = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"cim")==0) {
            inclConvImage = 1;
         } else if (strcasecmp(argv[iarg]+1,"nc")==0) {
            noClobber = 1;
         } else if (strcasecmp(argv[iarg]+1,"hki")==0) {
            doKerInfo = 1;
         } else if (strcasecmp(argv[iarg]+1,"oki")==0) {
            kernelImOut = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"sht")==0) {
            outShort = 1;
         } else if (strcasecmp(argv[iarg]+1,"obs")==0) {
	    sscanf(argv[++iarg], "%f", &outBscale);
         } else if (strcasecmp(argv[iarg]+1,"obz")==0) {
            sscanf(argv[++iarg], "%f", &outBzero);
         } else if (strcasecmp(argv[iarg]+1,"nsht")==0) {
            outNShort = 1;
         } else if (strcasecmp(argv[iarg]+1,"nbs")==0) {
	    sscanf(argv[++iarg], "%f", &outNiBscale);
         } else if (strcasecmp(argv[iarg]+1,"nbz")==0) {
            sscanf(argv[++iarg], "%f", &outNiBzero);
         } else if (strcasecmp(argv[iarg]+1,"fin")==0) {
            sscanf(argv[++iarg], "%f", &fillValNoise);
         } else if (strcasecmp(argv[iarg]+1,"convvar")==0) {
            convolveVariance = 1;
         } else if (strcasecmp(argv[iarg]+1,"allm")==0) {
            inclConvImage  = 1;
            inclSigmaImage = 1;
            inclNoiseImage = 1;
         } else if (strcasecmp(argv[iarg]+1,"v")==0) {
            sscanf(argv[++iarg], "%d", &verbose);
         } else if (strcasecmp(argv[iarg]+1,"savexy")==0) {
	    savexyflag = 1;
	    strcpy(xyfilename,argv[++iarg]);
	 } else if (strcasecmp(argv[iarg]+1,"ng")==0) {
	    sscanf(argv[++iarg], "%d", &ngauss);

	    deg_fixe    = (int *)realloc(deg_fixe,      ngauss*sizeof(int));
	    sigma_gauss = (float *)realloc(sigma_gauss, ngauss*sizeof(float));
	    
	    for (i = 0; i < ngauss; i++) {
	       sscanf(argv[++iarg], "%d", &deg_fixe[i]);
	       sscanf(argv[++iarg], "%f", &sigma_gauss[i]);
	       sigma_gauss[i] = (1.0/(2.0*sigma_gauss[i]*sigma_gauss[i]));
	    }
	 } else if (strcasecmp(argv[iarg]+1,"pca")==0) {
	    usePCA = 1;
	    sscanf(argv[++iarg], "%d", &ngauss); /* effectively ngauss... */

	    PCA = (float **)malloc(ngauss*sizeof(float *));
	    for (i = 0; i < ngauss; i++) {
	       scrStr = argv[++iarg];
	       
	       if (fits_open_file(&kPtr, scrStr, 0, &status) ||
		   fits_get_img_param(kPtr, MAXDIM, &kBitpix, &kNaxis, kNaxes, &status))
		  printError(status);

	       assert (kNaxes[0] == kNaxes[1]); /* square */
	       if (fwKernelPCA == (2*D_HWKERNEL+1))
		  fwKernelPCA = kNaxes[0];
	       else
		  assert(fwKernelPCA == kNaxes[0]);
	       
	       fprintf(stderr, "   Reading in basis image %s\n", scrStr);
	       
	       PCA[i] = (float *)malloc(kNaxes[0]*kNaxes[1]*sizeof(float));
	       if (fits_read_img_flt(kPtr, 0, 1, kNaxes[0]*kNaxes[1], 0, PCA[i], 0, &status) ||
		   fits_close_file(kPtr, &status))
		  printError(status);
	    }
	    
	    /* hack to not break other parts of the code */
	    deg_fixe    = (int *)realloc(deg_fixe,      ngauss*sizeof(int));
	    sigma_gauss = (float *)realloc(sigma_gauss, ngauss*sizeof(float));
	    for (i = 0; i < ngauss; i++) {
	       deg_fixe[i]    =  0;
	       sigma_gauss[i] = -1;
	    }
	    
         } else {
	    fprintf(stderr, "Unknown option : %s\n", argv[iarg]);
	    exit(1);
         }
      } else {
         fprintf(stderr, "Unexpected string encountered on command line : %s\n", argv[iarg]);
         exit(1);
      }
   }

   /* ngauss */
   /*
   for (i = 0; i < ngauss; i++) {
      fprintf(stderr, "%d %d %f\n", i, deg_fixe[i], sigma_gauss[i]);
   }
   */
   
   /* check that the requested options make sense */
   if ( !(image) ) {
      fprintf(stderr, "FATAL ERROR inim : required command line option\n");
      exit(1);
   }  
   if ( !(template) ) {
      fprintf(stderr, "FATAL ERROR tmplim : required command line option\n");
      exit(1);
   }  
   if ( !(outim) ) {
      fprintf(stderr, "FATAL ERROR outim : required command line option\n");
      exit(1);
   }  

   /* template thresholds and gain */
   if (tUThresh && (tUThresh <= 0)) {
      fprintf(stderr, "FATAL ERROR tu (%f) : cannot be negative (or zero)\n", tUThresh);
      exit(1);
   }
   if (tUThresh && tLThresh && (tUThresh < tLThresh)) {
      fprintf(stderr, "FATAL ERROR tu (%f) : cannot be smaller than tl (%f)\n", tUThresh, tLThresh);
      exit(1);
   }
   if (tGain && (tGain <= 0)) {
      fprintf(stderr, "FATAL ERROR tg (%f) : cannot be negative (or zero)\n", tGain);
      exit(1);
   }
   if (tRdnoise && (tRdnoise <= 0)) {
      fprintf(stderr, "FATAL ERROR tr (%f) : cannot be negative (or zero)\n", tRdnoise);
      exit(1);
   }
   
   /* image thresholds and gain */
   if (iUThresh && (iUThresh <= 0)) {
      fprintf(stderr, "FATAL ERROR iu (%f) : cannot be negative (or zero)\n", iUThresh);
      exit(1);
   }
   if (iUThresh && iLThresh && (iUThresh < iLThresh)) {
      fprintf(stderr, "FATAL ERROR iu (%f) : cannot be smaller than il (%f)\n", iUThresh, iLThresh);
      exit(1);
   }
   if (iGain && (iGain <= 0)) {
      fprintf(stderr, "FATAL ERROR ig (%f) : cannot be negative (or zero)\n", iGain);
      exit(1);
   }
   if (iRdnoise && (iRdnoise <= 0)) {
      fprintf(stderr, "FATAL ERROR ir (%f) : cannot be negative (or zero)\n", iRdnoise);
      exit(1);
   }
   
   /* kernel fit pars */
   if (hwKernel <= 0) {
      fprintf(stderr, "FATAL ERROR r (%d) : cannot be negative (or zero)\n", hwKernel);
      exit(1);
   }
   if (hwKernel == 1) {
      fprintf(stderr, "WARNING r (%d) : r = 1 is unstable\n", hwKernel);
   }
   if (kcStep && (kcStep <= 0)) {
      fprintf(stderr, "FATAL ERROR kcs (%d) : cannot be negative (or zero)\n", kcStep);
      exit(1);
   }
   if (kerFitThresh < 0) {
      fprintf(stderr, "FATAL ERROR ft (%f) : cannot be negative\n", kerFitThresh);
      exit(1);
   }
   if ((scaleFitThresh > 1) || (scaleFitThresh < 0)) {
      fprintf(stderr, "FATAL ERROR sft (%f) : cannot be negative or greater than 1\n", scaleFitThresh);
      exit(1);
   }
   if ((minFracGoodStamps > 1) || (minFracGoodStamps < 0) ) {
      fprintf(stderr, "FATAL ERROR nft (%f) : cannot be negative or greater than 1\n", minFracGoodStamps);
      exit(1);
   }
   if (kerSigReject < 0) {
      fprintf(stderr, "FATAL ERROR ks (%f) : cannot be negative\n", kerSigReject);
      exit(1);
   }

   /* masking pars */
   if (kfSpreadMask1 < 0) {
      fprintf(stderr, "FATAL ERROR mins (%f) : cannot be negative\n", kfSpreadMask1);
      exit(1);
   }

   /* asking for subsections of image */
   if ((gdXmin < 0) || (gdXmax < 0) || (gdXmin > gdXmax) || (gdYmin < 0) || (gdYmax < 0) || (gdYmin > gdYmax)) {
      fprintf(stderr, "FATAL ERROR gd : (%d,%d) < x,y < (%d,%d)\n", gdXmin, gdYmin, gdXmax, gdYmax);
      exit(1);
   }

   /* changing default region/stamp size */
   if (regFile && sstampFile) {
      fprintf(stderr, "FATAL ERROR rf (%s) and ssf (%s) : conflicting requests\n", regFile, sstampFile);
      exit(1);
   }

   /* region keyword */
   if (regKeyWord && regFile) {
      fprintf(stderr, "FATAL ERROR rkw (%s) and rf (%s) : conflicting requests\n", regKeyWord, regFile);
      exit(1);
   }
   if (regKeyWord) {
      if (numRegKeyWord < 1) {
         fprintf(stderr, "FATAL ERROR rkw (%s) : number of regions invalid (%d)\n", regKeyWord, numRegKeyWord);
         exit(1);
      }
   }
   
   /* substamp properties */
   if (nKSStamps < 0) {
      fprintf(stderr, "FATAL ERROR nss (%d) : cannot be negative\n", nKSStamps);
      exit(1);
   }
   if (hwKSStamp < 0) {
      fprintf(stderr, "FATAL ERROR nss (%d) : cannot be negative\n", hwKSStamp);
      exit(1);
   }
   if (hwKSStamp < hwKernel) {
      fprintf(stderr, "FATAL ERROR rss (%d) < r (%d) : conflicting requests\n", hwKSStamp, hwKernel);
      exit(1);
   }

   /* defining convolution or normalization */
   if (forceConvolve &&
       ( !((strncmp(forceConvolve, "i", 1)==0) || (strncmp(forceConvolve, "t", 1)==0) || (strncmp(forceConvolve, "b", 1)==0)) )) {
      fprintf(stderr, "FATAL ERROR c (%s) : invalid char ('i' or 't' or 'b')\n", forceConvolve);
      exit(1);
   }
   if (photNormalize && ( !((strncmp(photNormalize, "i", 1)==0) || (strncmp(photNormalize, "t", 1)==0) || (strncmp(photNormalize, "u", 1)==0)) )) {
      fprintf(stderr, "FATAL ERROR n (%s) : invalid char ('i', 't', or 'u')\n", photNormalize);
      exit(1);
   }

   /* good figure of merit request */
   if (figMerit && ( !((strncmp(figMerit, "s", 1)==0) || (strncmp(figMerit, "h", 1)==0) || (strncmp(figMerit, "v", 1)==0)) )) {
      fprintf(stderr, "FATAL ERROR fom (%s) : invalid char ('v', 's' or 'h')\n", figMerit);
      exit(1);
   }

   /* fit orders */
   if (kerOrder < 0) {
      fprintf(stderr, "FATAL ERROR ko (%d) : cannot be negative\n", kerOrder);
      exit(1);
   }
   if (bgOrder < 0) {
      fprintf(stderr, "FATAL ERROR bgo (%d) : cannot be negative\n", bgOrder);
      exit(1);
   }

   /* bscale, bzero */
   if (outBscale && (outBscale == 0.)) {
      fprintf(stderr, "FATAL ERROR obs (%f) : cannot be zero\n", outBscale);
      exit(1);
   }      
   
   if (outNiBscale && (outNiBscale == 0.)) {
      fprintf(stderr, "FATAL ERROR nbs (%f) : cannot be zero\n", outNiBscale);
      exit(1);
   }      

   /* PCA */
   if (usePCA) {
      hwKernel = fwKernelPCA / 2;
      fprintf(stderr, "ATTN : size of input kernel images fixes hwKernel to %d\n", hwKernel);
   }
   
   return;
}
