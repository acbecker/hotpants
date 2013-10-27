#include<stdio.h>
#include<string.h>
#include<strings.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<fitsio.h>

#define max(x,y) x>y?x:y
#define min(x,y) x<y?x:y

/* globals */
int       ngauss=3, *deg_fixe=NULL, dofullImage=0;
float     *sigma_gauss=NULL;
int       fwKernel,nCompKer,nComp,nBGVectors,nCompTotal,kerOrder,bgOrder,nR;
double    *filter_x,*filter_y,**kernel_vec,*kernel_coeffs,*kernel,kSumIm;
char      *inConv = NULL;
long      oNaxes[2];
char      **regions = NULL;
float     xConv=1e10, yConv=1e10;
int       rPixX, rPixY;

/* funcs */
void   getKernelInfo(char *);
void   getKernel(char *, double **, int, int);
void   fits_get_kernel_btbl(fitsfile *, double **, int);
     
void   getKernelVec();
double *kernel_vector(int, int, int, int, int *);
double make_kernel(int, int, double *);
void   spatial_convolve(float *, int, int, float, float, double *, float *);
double get_background(int, int, double *);
void   printError(int);

int main(int argc, char **argv)
{
   
   int    iarg, i, j, xsize, ysize, ndelta;
   int    status=0;
   float  numKerFW=2;
   double kSum;
   char   *diffim = NULL, *outConv = NULL;
   float  *delta = NULL, *dconv = NULL;
   double *kernelSol = NULL;
   long   cNaxes[2];
   char   scrStr[256], help[4096];
   fitsfile *fPtr;

   sprintf(help, "Usage : extractkern [options] diffimage outimage\n");
   sprintf(help, "%sOptions:\n", help);
   sprintf(help, "%s   [-xy x y]        : convolve kernel with delta function at x,y\n", help);
   sprintf(help, "%s   [-nkw numkwidth] : # kernel widths for outimage size (%.1f)\n", help, numKerFW);
   sprintf(help, "%s   [-a]             : sample entire diffimage size with delta functions\n", help);
   sprintf(help, "%s   [-im image]      : convolve fitsfile instead of delta function\n\n", help);
   sprintf(help, "%s   To be used in conjuntion with the diffimage produced by hotpants\n", help);   
   sprintf(help, "%s      using the -hki option.  [-xy] convolves a delta function at\n", help);
   sprintf(help, "%s      the image position x, y with the spatially varying kernel\n", help);
   sprintf(help, "%s      used in the hotpants convolution.  Provides a visual realization \n", help);
   sprintf(help, "%s      of the kernel at that position, and can be useful for cosmic ray\n", help);
   sprintf(help, "%s      discrimination.  Also, if used with the [-im] option, one may\n", help);
   sprintf(help, "%s      reconstruct the entire convolved image to avoid storing it on disk.\n", help);
   
   /* read in command options. j counts # of required args given */
   for (iarg=1, j=0; iarg < argc; iarg++) {
      if (argv[iarg][0]=='-') {
         if (strcasecmp(argv[iarg]+1,"xy")==0) {
            sscanf(argv[++iarg], "%f", &xConv);
            sscanf(argv[++iarg], "%f", &yConv);
         } else if (strcasecmp(argv[iarg]+1,"nkw")==0) {
            sscanf(argv[++iarg], "%f", &numKerFW);
         } else if (strcasecmp(argv[iarg]+1,"a")==0) {
            dofullImage = 1;
	    xConv = 0;
	    yConv = 0;
         } else if (strcasecmp(argv[iarg]+1,"im")==0) {
            inConv = argv[++iarg];
	    xConv = 0;
	    yConv = 0;
         } else {
            fprintf(stderr, "Unknown option %s\n", argv[iarg]);
            exit(1);
         }
      } else {
         diffim  = argv[iarg++];
         outConv = argv[iarg++];
      }
   }
   if (iarg < 2) {
      /* not enough command line images...*/
      fprintf(stderr, "%s\n", help);
      exit(1);
   }

   /* insanity checking */
   if ( ((xConv == 1e10) || (yConv == 1e10)) && !(dofullImage) && !(inConv) ) {
      fprintf(stderr, "Sorry, I do not know what to do, exiting...\n");
      exit(1);
   }

   deg_fixe    = (int *)calloc(ngauss, sizeof(int));
   sigma_gauss = (float *)calloc(ngauss, sizeof(float));
   
   /* fill up global kernel info */
   getKernelInfo(diffim);
   
   /* set array and comp sizes */
   nCompKer = 0;
   for (i = 0; i < ngauss; i++)
      nCompKer += ((deg_fixe[i] + 1) * (deg_fixe[i] + 2)) / 2;
   
   nComp       = ((kerOrder + 1) * (kerOrder + 2)) / 2;
   nBGVectors  = ((bgOrder + 1) * (bgOrder + 2)) / 2;
   nCompTotal  = nCompKer * nComp + nBGVectors;
   kernelSol   = (double *)realloc(kernelSol, (nCompTotal+1)*sizeof(double));
   
   /* allocate more kernel vectors */
   if ( !( filter_x      = (double *)malloc(nCompKer*fwKernel*sizeof(double))) ||
        !( filter_y      = (double *)malloc(nCompKer*fwKernel*sizeof(double))) ||
        !( kernel        = (double *)malloc(fwKernel*fwKernel*sizeof(double))) ||
        !( kernel_vec    = (double **)malloc(nCompKer*sizeof(double *))) ||
        !( kernel_coeffs = (double *)malloc(nCompKer*sizeof(double))) )
      exit(3);
   
   /* Set output image size */
   if (! (dofullImage) ) {
      xsize = numKerFW * fwKernel;
      ysize = numKerFW * fwKernel;
      
      oNaxes[0] = xsize;
      oNaxes[1] = ysize;
   }
   else {
      xsize = oNaxes[0];
      ysize = oNaxes[1];
   }
   
   /* you can input a shape to be convolved, theoretically... */
   if (inConv) {
      if ( fits_open_file(&fPtr, inConv, 0, &status)  ||
           fits_get_img_param(fPtr, 2, NULL, NULL, cNaxes, &status) )
         printError(status);
      
      delta = (float *)realloc(delta, cNaxes[0]*cNaxes[1]*sizeof(float));
      if ( fits_read_img_flt(fPtr, 1, 1, cNaxes[0]*cNaxes[1], 0, delta, 0, &status) ||
           fits_close_file(fPtr, &status) )
         printError(status);

      xsize = cNaxes[0];
      ysize = cNaxes[1];
      oNaxes[0] = cNaxes[0];
      oNaxes[1] = cNaxes[1];
   }
   else {
      delta = (float *)realloc(delta, xsize*ysize*sizeof(float));
      memset(delta, 0, xsize*ysize*sizeof(float));
      
      ndelta = 0;
      for (j = fwKernel-1; j < (ysize-1); j += fwKernel) {
         for (i = fwKernel-1; i < (xsize-1); i += fwKernel) {
            /* delta function in middle */
            delta[i + xsize*j] = 1.;
            ndelta += 1;
         }
      }
   }
   
   /* output convolved image */
   dconv = (float *)realloc(dconv, xsize*ysize*sizeof(float));
   memset(dconv, 0, xsize*ysize*sizeof(float));

   /* fill weight matrices kernel_vec, calls kernel_vector */
   getKernelVec();

   /* read kernel solution */
   getKernel(diffim, &kernelSol, xsize, ysize);
   
   /* do the convolution, fills kernel and calls make_kernel */
   spatial_convolve(delta, xsize, ysize, xConv, yConv, kernelSol, dconv);

   /* clobber output image */
   sprintf(scrStr, "!%s", outConv);
   /* create and open new empty output FITS file, using input image as template.*/
   if ( fits_create_file(&fPtr, scrStr, &status)                            ||
        fits_create_img(fPtr, FLOAT_IMG, 2, oNaxes, &status)                ||
        fits_write_img_flt(fPtr, 1, 1, xsize*ysize, dconv, &status) ||
        fits_close_file(fPtr, &status) )
      printError(status);
   
   /* sanity check - add up all pixels in the image */
   kSum = 0;
   for (i = 0; i < xsize*ysize; i++) {
      kSum += dconv[i];
   }
   fprintf(stderr, " Actual sum of pixels in convolved image : %.6f\n", kSum);
   fprintf(stderr, " Kernel Sum from input image             : %.6f\n", kSumIm);
   
   
   for (i = 0; i < 10; i++)
      if (regions[i]) free(regions[i]);
   if (regions)       free(regions);
   if (kernelSol)     free(kernelSol);
   if (delta)         free(delta);
   if (dconv)         free(dconv);
   if (filter_x)      free(filter_x);
   if (filter_y)      free(filter_y);
   if (kernel)        free(kernel);
   if (kernel_vec)    free(kernel_vec);
   if (kernel_coeffs) free(kernel_coeffs);

   return 1;
}


/* ********************************** */
/* from functions.c */
/* ********************************** */

void getKernelInfo(char *kimage) {
   /*****************************************************
     Get all 1-time info from kernel fits header, overriding defaults
      and command line options.
   *****************************************************/
   
   fitsfile *kPtr;
   int i, existsTable, status = 0;
   char hKeyword[1024];

   /* open the input kernel image */
   if ( fits_open_file(&kPtr, kimage, 0, &status) )
      printError(status);

   /* required keyword in primary HDU */
   if ( fits_read_key_log(kPtr, "KERINFO", &existsTable, NULL, &status) )
      printError(status);
   
   if (!(existsTable)) {
      fits_close_file(kPtr, &status);
      fprintf(stderr, "This image does not appear to contain a kernel table, exiting...\n");
      exit(1);
   }

   /* move to binary kernel table... */
   if ( fits_get_num_hdus(kPtr, &existsTable, &status) ||
        fits_movabs_hdu(kPtr, existsTable, NULL, &status) ||
        fits_read_key(kPtr, TINT,    "NGAUSS",  &ngauss, NULL, &status) ||
        fits_read_key(kPtr, TINT,    "FWKERN",  &fwKernel, NULL, &status) ||
        fits_read_key(kPtr, TINT,    "CKORDER", &kerOrder, NULL, &status) ||
        fits_read_key(kPtr, TINT,    "BGORDER", &bgOrder, NULL, &status) )
      printError(status);

   deg_fixe    = (int *)realloc(deg_fixe,      ngauss*sizeof(int));
   sigma_gauss = (float *)realloc(sigma_gauss, ngauss*sizeof(float));
   
   /* read kernel gaussian info */
   for (i = 0; i < ngauss; i++) {
      sprintf(hKeyword, "DGAUSS%d", i+1);
      if (fits_read_key(kPtr, TINT, hKeyword, &deg_fixe[i], NULL, &status))
         printError(status);
      sprintf(hKeyword, "SGAUSS%d", i+1);
      if (fits_read_key(kPtr, TFLOAT, hKeyword, &sigma_gauss[i], NULL, &status))
         printError(status);
   }
   
   if (fits_close_file(kPtr, &status) )
      printError(status);
   return;
}

void getKernel(char *kimage, double **kerSol, int xsize, int ysize) {
   
   /* read in kernel image for region */

   fitsfile *fPtr;
   int status = 0, i;
   char hKeyword[1024];
   int  rXMin, rXMax, rYMin, rYMax;
   int  rXBMin, rXBMax, rYBMin, rYBMax;
   int hwKernel = fwKernel/2;
   
   /* open the input kernel image */
   if ( fits_open_file(&fPtr, kimage, 0, &status) )
      printError(status);

   /* get its size if needed */
   if ( dofullImage || inConv )
      if ( fits_get_img_param(fPtr, 2, NULL, NULL, oNaxes, &status) )
         printError(status);

   /* grab all regions in the primary image header, up to 10 in number */
   regions = (char **)malloc(10*sizeof(char *));
   for (i = 0; i < 10; i++)
      regions[i] = (char *)malloc(80*sizeof(char));
   if ( fits_read_keys_str(fPtr, "REGION", 0, 9, regions, &nR, &status ) )
      printError(status);
   for (i = 0; i < nR; i++) {
      if (sscanf(regions[i], "[%d:%d,%d:%d]", &rXMin, &rXMax, &rYMin, &rYMax) != 4) {
         fprintf(stderr, "Problem with region %d (%s), exiting...\n", i, regions[i]);
         exit(1);
      }
      if ( ((xConv>-1) && (yConv>-1)) && ((xConv+1) >= rXMin &&
					  (xConv+1) <= rXMax &&
					  (yConv+1) >= rYMin &&
					  (yConv+1) <= rYMax)) {
         /* we got our guy, lets roll! */
         sprintf(hKeyword, "KSUM%02d", i);
         if (fits_read_key(fPtr, TDOUBLE, hKeyword, &kSumIm, NULL, &status))
            printError(status);

	 /* set region sizes */
         rXBMin = max(0, rXMin - hwKernel);
         rYBMin = max(0, rYMin - hwKernel); 
         rXBMax = min(xsize, rXMax + hwKernel);
         rYBMax = min(ysize, rYMax + hwKernel);
	 rPixX  = rXBMax - rXBMin + 1;
	 rPixY  = rYBMax - rYBMin + 1;

         break;
      }
      else {
	 rXMin = rXMax = rYMin = rYMax = 0;
      }
   }
   if ( (rXMin == 0) && (rXMax == 0) && (rYMin == 0) && (rYMax == 0) ) {
      /* we did not get our guy, shucks... */
      fprintf(stderr, "Unable to locate appropriate region for %.2f, %.2f, exiting...\n", xConv, yConv);
      exit(2);
   }

   fits_get_kernel_btbl(fPtr, &(*kerSol), i);
   return;
}

void fits_get_kernel_btbl(fitsfile *kPtr, double **kernelSol, int nRegion) {
   int status=0, existsTable;
   
   /* move to binary kernel table... */
   if ( fits_get_num_hdus(kPtr, &existsTable, &status) ||
        fits_movabs_hdu(kPtr, existsTable, NULL, &status) )
      printError(status);

   memset(*kernelSol, 0, (nCompTotal+1)*sizeof(double));
   if (fits_read_col(kPtr, TDOUBLE, nRegion+1, 1, 1, (nCompTotal+1), 0, *kernelSol, 0, &status))
      printError(status);

   return;
}

/* ********************************** */
/* from alard.c */
/* ********************************** */

void getKernelVec() {
   /*****************************************************
    * Fills kernel_vec with kernel weight filter, called only once
    *****************************************************/
   int ig, idegx, idegy, nvec;
   int ren;
   
   nvec = 0;
   for (ig = 0; ig < ngauss; ig++) {
      for (idegx = 0; idegx <= deg_fixe[ig]; idegx++) {
         for (idegy = 0; idegy <= deg_fixe[ig]-idegx; idegy++) {
            /* stores kernel weight mask for each order */
            kernel_vec[nvec] = kernel_vector(nvec, idegx, idegy, ig, &ren);
            nvec++;
         }
      }
   }
}

double *kernel_vector(int n, int deg_x, int deg_y, int ig, int *ren) {
   /*****************************************************
    * Creates kernel sized entry for kernel_vec for each kernel degree 
    *   Mask of filter_x * filter_y, filter = exp(-x**2 sig) * x^deg 
    *   Subtract off kernel_vec[0] if n > 0
    * NOTE: this does not use any image
    ******************************************************/

   double    *vector=NULL,*kernel0=NULL;
   int       i,j,k,dx,dy,ix;
   double    sum_x,sum_y,x,qe;
   
   vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
   dx = (deg_x / 2) * 2 - deg_x;
   dy = (deg_y / 2) * 2 - deg_y;
   sum_x = sum_y = 0.0;
   *ren = 0;
   
   for (ix = 0; ix < fwKernel; ix++) {
      x            = (double)(ix - fwKernel/2);
      k            = ix+n*fwKernel;
      qe           = exp(-x * x * sigma_gauss[ig]);
      filter_x[k]  = qe * pow(x, deg_x);
      filter_y[k]  = qe * pow(x, deg_y);
      sum_x       += filter_x[k];
      sum_y       += filter_y[k];
   }
   
   if (n > 0)
      kernel0 = kernel_vec[0];
   
   sum_x = 1. / sum_x;
   sum_y = 1. / sum_y;
   
   if (dx == 0 && dy == 0) {
      for (ix = 0; ix < fwKernel; ix++) {
         filter_x[ix+n*fwKernel] *= sum_x;
         filter_y[ix+n*fwKernel] *= sum_y;
      }
      
      for (i = 0; i < fwKernel; i++) {
         for (j = 0; j < fwKernel; j++) {
            vector[i+fwKernel*j] = filter_x[i+n*fwKernel] * filter_y[j+n*fwKernel];
         }
      }
      
      if (n > 0) {
         for (i = 0; i < fwKernel * fwKernel; i++) {
            vector[i] -= kernel0[i];
         }
         *ren = 1;
      }
   } else {
      for (i = 0; i < fwKernel; i++) {
         for (j = 0; j < fwKernel; j++) {
            vector[i+fwKernel*j] = filter_x[i+n*fwKernel] * filter_y[j+n*fwKernel];
         }
      }
   }
   return vector;
}

double make_kernel(int xi, int yi, double *kernelSol) {
   /*****************************************************
    * Create the appropriate kernel at xi, yi, return sum
    *****************************************************/
   
   int    i1,k,ix,iy,i;
   double ax,ay,sum_kernel;
   double xf, yf;
   
   k = 2;
   /* RANGE FROM -1 to 1 */
   xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
   yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
   
   for (i1 = 1; i1 < nCompKer; i1++) {
      kernel_coeffs[i1] = 0.0;
      ax = 1.0;
      for (ix = 0; ix <= kerOrder; ix++) {
         ay = 1.0;
         for (iy = 0; iy <= kerOrder - ix; iy++) {
            kernel_coeffs[i1] += kernelSol[k++] * ax * ay;
            ay *= yf;
         }
         ax *= xf;
      }
   }
   kernel_coeffs[0] = kernelSol[1]; 
   
   for (i = 0; i < fwKernel * fwKernel; i++)
      kernel[i] = 0.0;
   
   sum_kernel = 0.0;
   for (i = 0; i < fwKernel * fwKernel; i++) {
      for (i1 = 0; i1 < nCompKer; i1++) {
         kernel[i] += kernel_coeffs[i1] * kernel_vec[i1][i];
      }
      sum_kernel += kernel[i];    
   }
   return sum_kernel;
}


void spatial_convolve(float *image, int xSize, int ySize, float xConv, float yConv, double *kernelSol, float *outim) {
   /*****************************************************
    * Take image and convolve it using the kernelSol every kernel width
    *****************************************************/
  
   int       i1,j1,i2,j2,nsteps_x,nsteps_y,i,j,i0,j0,ic,jc,ik,jk;
   double    q;
   int       x0, y0, hwKernel = fwKernel/2;
   
   nsteps_x = ceil((double)(xSize)/(double)fwKernel);
   nsteps_y = ceil((double)(ySize)/(double)fwKernel);

   x0 = max(0, (int)(xConv - xSize/2));
   y0 = max(0, (int)(yConv - ySize/2));

   for (j1 = 0; j1 < nsteps_y; j1++) {
      j0 = j1 * fwKernel + hwKernel;
      
      for(i1 = 0; i1 < nsteps_x; i1++) {
         i0 = i1 * fwKernel + hwKernel;
	 
         make_kernel(x0 + i0 + hwKernel, y0 + j0 + hwKernel, kernelSol);
	 
         for (j2 = 0; j2 < fwKernel; j2++) {
            j = j0 + j2;
            if (j >= (ySize - hwKernel)) break;
            
            for (i2 = 0; i2 < fwKernel; i2++) {
               i = i0 + i2;
               if (i >= (xSize - hwKernel)) break;

               q    = 0;
               for (jc = j - hwKernel; jc <= j + hwKernel; jc++) {
                  jk = j - jc + hwKernel;
                  
                  for (ic = i - hwKernel; ic <= i + hwKernel; ic++) {
                     ik = i - ic + hwKernel;

		     q += image[ic+xSize*jc] * kernel[ik+jk*fwKernel];
		  }
	       }
	       outim[i+xSize*j] = q;
	    }
	 }
      }
   }
   return;
}

double get_background(int xi, int yi, double *kernelSol) {
   /*****************************************************
    * Return background value at xi, yi
    *****************************************************/
   
   double  background,ax,ay,xf,yf;
   int     i,j,k;
   int     ncompBG;
   
   ncompBG = (nCompKer - 1) * ( ((kerOrder + 1) * (kerOrder + 2)) / 2 ) + 1;
   
   background = 0.0;
   k          = 1;
   /* RANGE FROM -1 to 1 */
   xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
   yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
   
   ax=1.0;
   for (i = 0; i <= bgOrder; i++) {
      ay = 1.0; 
      for (j = 0; j <= bgOrder - i; j++) {
         background += kernelSol[ncompBG+k++] * ax * ay;
         ay *= yf;
      }
      ax *= xf;
   }
   return background;
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
