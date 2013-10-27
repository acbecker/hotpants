#include<assert.h>
#include<stdio.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<fitsio.h>

#define D_FILLVAL 0
#define D_MASK    0x8000
#define MAXDIM    2

void   printError(int);

int main(int argc, char **argv)
{

   int i, iarg, nPix, anynul, status, mask;
   char help[4096];
   
   fitsfile *iPtr;
   int iBitpix, iNaxis;
   long iNaxes[MAXDIM];
   float *iRData=NULL;
   
   fitsfile *mPtr;
   int mBitpix, mNaxis;
   long mNaxes[MAXDIM];
   int *mRData=NULL;

   char  *inim=NULL, *maskim=NULL, *outim=NULL;
   float fillVal=D_FILLVAL;

   /* START */
   status  = 0;
   mask    = D_MASK;
   fillVal = D_FILLVAL;
   
   sprintf(help, "Usage : maskim [options]\n");
   sprintf(help, "%sRequired options:\n", help);
   sprintf(help, "%s   [-inim fitsfile]  : input image\n", help);
   sprintf(help, "%s   [-maskim fitsfile]: corresponding mask image\n", help);
   sprintf(help, "%s   [-outim fitsfile] : output masked image\n\n", help);
   sprintf(help, "%sOptions:\n", help);
   sprintf(help, "%s   [-fi fillval]    : Value of filled pixels (%.3f)\n", help, fillVal);
   sprintf(help, "%s   [-mask hex]      : Mask for bad pixels (0x%x)\n", help, mask);

   /* read in command options. j counts # of required args given */
   for (iarg=1; iarg < argc; iarg++) {
      if (argv[iarg][0]=='-') {
         if (strcasecmp(argv[iarg]+1,"inim")==0) {
            inim = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"maskim")==0) {
            maskim = argv[++iarg];
         } else if (strcasecmp(argv[iarg]+1,"outim")==0) {
            outim = argv[++iarg];
	 } else if (strcasecmp(argv[iarg]+1,"fi")==0) {
	    sscanf(argv[++iarg], "%f", &fillVal);
         } else if (strcasecmp(argv[iarg]+1,"mask")==0) {
	    sscanf(argv[++iarg], "%d", &mask);
	 } else {
	    fprintf(stderr, "Unknown option : %s\n", argv[iarg]);
	    exit(1);
         }
      } else {
         fprintf(stderr, "Unexpected string encountered on command line : %s\n", argv[iarg]);
         exit(1);
      }
   }
   /* read in command options */
   if (argc < 2) {
      fprintf(stderr, "%s\n", help);
      exit(1);
   }

   /* insanity checking */
   if ( !(inim) ) {
      fprintf(stderr, "FATAL ERROR inim is required command line option\n");
      exit(1);
   }  
   if ( !(maskim) ) {
      fprintf(stderr, "FATAL ERROR maskim is required command line option\n");
      exit(1);
   }  
   if ( !(outim) ) {
      fprintf(stderr, "FATAL ERROR outim is required command line option\n");
      exit(1);
   }  

   /* open up, get bitpix, # dimensions, image size */
   if ( fits_open_file(&iPtr, inim, 0, &status) ||
        fits_get_img_param(iPtr, 2, &iBitpix, &iNaxis, iNaxes, &status) ) {
      printError(status);
   }
   if ( fits_open_file(&mPtr, maskim, 0, &status) ||
        fits_get_img_param(mPtr, 2, &mBitpix, &mNaxis, mNaxes, &status) ) {
      printError(status);
   }

   assert(iNaxes[0] == mNaxes[0]);
   assert(iNaxes[1] == mNaxes[1]);

   nPix = iNaxes[0]*iNaxes[1];
   
   iRData = (float *)realloc(iRData, nPix*sizeof(float));
   mRData = (int *)realloc(mRData, nPix*sizeof(int));
   if (iRData == NULL || mRData == NULL) {
      fprintf(stderr, "Cannot Allocate Standard Data Arrays\n"); 
      exit (1);
   }
   memset(iRData,   0.0, nPix*sizeof(float));
   memset(mRData,   0x0, nPix*sizeof(int));

   if (fits_read_img(iPtr, TFLOAT, 1, nPix, 0, iRData, &anynul, &status) ||
       fits_read_img(mPtr, TINT,   1, nPix, 0, mRData, &anynul, &status) ||
       fits_close_file(mPtr, &status))
      printError(status);

   /* do the masking */
   for (i = nPix; i--; )
      if (mRData[i] & mask)
	 iRData[i] = fillVal;
   
   /* reuse mptr and help */
   sprintf(help, "!%s", outim);
   if (fits_create_template(&mPtr, help, inim, &status) ||
       fits_write_img(mPtr, TFLOAT, 1, nPix, iRData, &status) ||
       fits_close_file(mPtr, &status) ||
       fits_close_file(iPtr, &status))
      printError(status);      

   free(iRData);
   free(mRData);

   return 0;
}

void printError(int status) {
   /*****************************************************/
   /* Print out cfitsio error messages and exit program */
   /*****************************************************/
   if (status) {
      fits_report_error(stderr, status); /* print error report */
      exit( status );    /* terminate the program, returning error status */
   }
   return;
}
