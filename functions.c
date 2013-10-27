#include<stdio.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<fitsio.h>
#include<ctype.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/*
  
Some of these subroutines appear originally in code created by Gary
Bernstein for psfmatch, but have been modified and/or rewritten for
the current software package.  In particular, the determination of
image noise statistics has been taken directly from the psfmatch
code.

08/20/01 acbecker@lucent.com

*/

/* armin */ 
void loadxyfile(char *filename, int cmpfileflag){      
    FILE *xyfile;
    int Nalloc,c;
    char line[SCRLEN];
    xyfile   = fopen(filename, "r");
    fprintf(stderr, "WARNING : INPUT FORMAT HARDCODED : X=1, Y=2; 1-indexed coordinates\n");
    if (cmpfileflag) {
        for (;;){
            c=getc(xyfile);
            if (c== EOF) break;	 
            if (c=='\n') break;	  
        }
    }
    Nalloc=100;
    if ( !(xcmp = (float *)malloc(Nalloc*sizeof(float))) ||
         !(ycmp = (float *)malloc(Nalloc*sizeof(float))) )
        exit(1);
    Ncmp=0;
    while(fgets(line, SCRLEN, xyfile) != NULL) {
        /* reallocate if necessary */
        if (Ncmp==Nalloc){
            Nalloc+=100;
            xcmp  = (float *)realloc(xcmp, Nalloc *sizeof(float));
            ycmp  = (float *)realloc(ycmp, Nalloc *sizeof(float));
        }
        
        /* skip poorly formatted lines... */
        if (isalpha(line[0])) continue;
        if (sscanf(line, "%f %f", &xcmp[Ncmp], &ycmp[Ncmp]) == 2) {
            Ncmp++;
        }
    }
    fclose(xyfile);
    
    /* take care of 1-indexing, our arrays are 0-indexed */
    for (c = 0; c < Ncmp; c++) {
        xcmp[c] -= 1;
        ycmp[c] -= 1;
    }
}

/*armin*/
void savexy(stamp_struct *stamps, int nStamps, long xmin, long ymin, int regioncounter){
    int  sscnt,istamp;
    FILE *xyfileused,*xyfileall,*xyfileskipped;
    char xyfilenameall[1000];
    char xyfilenameskipped[1000];
    
    sprintf(xyfilenameall,"%s.all",xyfilename);
    sprintf(xyfilenameskipped,"%s.skipped",xyfilename);
    
    if (regioncounter==0) {
        xyfileused = fopen(xyfilename, "w");
        xyfileall  = fopen(xyfilenameall, "w");
        xyfileskipped = fopen(xyfilenameskipped, "w");
    } else {
        xyfileused = fopen(xyfilename, "a");
        xyfileall  = fopen(xyfilenameall, "a");
        xyfileskipped = fopen(xyfilenameskipped, "a");
    }
    
    fprintf(xyfileall,    "#%4s %4s\n", "X", "Y");    
    fprintf(xyfileused,   "#%4s %4s\n", "X", "Y");    
    fprintf(xyfileskipped,"#%4s %4s\n", "X", "Y");
    
    /* note single pixel offset for ds9 display compared to 0-index array */
    for (istamp = 0; istamp < nStamps; istamp++) {
        for (sscnt = 0; sscnt < stamps[istamp].nss; sscnt++) {
            /*fprintf(xyfileall,        " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
            fprintf(xyfileall,        "image;box(%4ld,%4ld,%d,%d) # color=yellow\n",
                    stamps[istamp].xss[sscnt] + xmin + 1,
                    stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);    
            
            if (sscnt == stamps[istamp].sscnt)
                /*fprintf(xyfileused,    " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
                fprintf(xyfileused,    "image;box(%4ld,%4ld,%d,%d) # color=green\n",
                        stamps[istamp].xss[sscnt] + xmin + 1,
                        stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);
            
            if (sscnt < stamps[istamp].sscnt)
                /*fprintf(xyfileskipped, " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
                fprintf(xyfileskipped, "image;box(%4ld,%4ld,%d,%d) # color=red\n",
                        stamps[istamp].xss[sscnt] + xmin + 1,
                        stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);    
        }
    }
    
    fclose(xyfileskipped);
    fclose(xyfileall);
    fclose(xyfileused);
}

int allocateStamps(stamp_struct *stamps, int nStamps) {                                                                     
    int i,j;
    int nbgVectors;
    
    nbgVectors = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    if (stamps) {
        for (i = 0; i < nStamps; i++) {
            
            /* **************************** */
            if(!(stamps[i].vectors = (double **)calloc((nCompKer+nbgVectors), sizeof(double *))))
                return 1;
            
            for (j = 0; j < nCompKer + nbgVectors; j++) {
                if(!(stamps[i].vectors[j] = (double *)calloc(fwKSStamp*fwKSStamp, sizeof(double))))
                    return 1;
            }
            
            /* **************************** */
            
            if (!(stamps[i].krefArea      = (double *)calloc(fwKSStamp*fwKSStamp, sizeof(double))))
                return 1;
            
            /* **************************** */
            
            if (!(stamps[i].mat    = (double **)calloc(nC, sizeof(double *))))
                return 1;
            
            for (j = 0; j < nC; j++)
                if ( !(stamps[i].mat[j] = (double *)calloc(nC, sizeof(double))) )
                    return 1;
            
            /* **************************** */
            
            if (!(stamps[i].xss = (int *)calloc(nKSStamps, sizeof(int))))
                return 1;
            
            /* **************************** */
            
            if (!(stamps[i].yss = (int *)calloc(nKSStamps, sizeof(int))))
                return 1;
            
            /* **************************** */
            
            if (!(stamps[i].scprod = (double *)calloc(nC, sizeof(double))))
                return 1;
            
            /* **************************** */
            
            stamps[i].x0 = stamps[i].y0 = stamps[i].x = stamps[i].y = 0;
            stamps[i].nss = stamps[i].sscnt = 0;
            stamps[i].nx = stamps[i].ny = 0;
            stamps[i].sum = stamps[i].mean = stamps[i].median = 0;
            stamps[i].mode = stamps[i].sd = stamps[i].fwhm = 0;
            stamps[i].lfwhm = stamps[i].chi2 = 0;
            stamps[i].norm = stamps[i].diff = 0;
        }
    }
    return 0;
}


void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int *niS, int *ntS,
                 int getCenters, int rXBMin, int rYBMin, stamp_struct *ciStamps, stamp_struct *ctStamps,
                 float *iRData, float *tRData, float hardX, float hardY) {
    
    int sPixX, sPixY;
    int xmax, ymax, k, l, nr2;
    int nss;
    /*int bbitt1=0x100, bbitt2=0x200, bbiti1=0x400, bbiti2=0x800;*/
    int bbitt1=FLAG_T_BAD, bbitt2=FLAG_T_SKIP, bbiti1=FLAG_I_BAD, bbiti2=FLAG_I_SKIP;
    float *refArea=NULL;
    double check;
    
    if (verbose >= 1)
        fprintf(stderr, "    Stamp in region : %d:%d,%d:%d\n",
                sXMin, sXMax, sYMin, sYMax);
    
    /* global vars */
    sPixX  = sXMax - sXMin + 1;
    sPixY  = sYMax - sYMin + 1;
    
    if (!(strncmp(forceConvolve, "i", 1)==0)) {
        
        if (ctStamps[*ntS].nss == 0) {
            refArea = (float *)calloc(sPixX*sPixY, sizeof(float));
            
            /* temp : store the whole stamp in refArea */
            cutStamp(tRData, refArea, rPixX, sXMin-rXBMin, sYMin-rYBMin,
                     sXMax-rXBMin, sYMax-rYBMin, &ctStamps[*ntS]);
            
            if ( !( getStampStats3(refArea, ctStamps[*ntS].x0, ctStamps[*ntS].y0, sPixX, sPixY,
                                   &ctStamps[*ntS].sum, &ctStamps[*ntS].mean, &ctStamps[*ntS].median,
                                   &ctStamps[*ntS].mode, &ctStamps[*ntS].sd, &ctStamps[*ntS].fwhm,
                                   &ctStamps[*ntS].lfwhm, 0x0, 0xffff, 3))) {
                /* pointless */
                if (verbose >= 1)
                    fprintf(stderr, "    Tmpl  xs : %4i ys : %4i  (sky,dsky = %.1f,%.1f)\n",
                            ctStamps[*ntS].x, ctStamps[*ntS].y, ctStamps[*ntS].mode, ctStamps[*ntS].fwhm);
            }
            free(refArea);
        }
    }
    
    if (!(strncmp(forceConvolve, "t", 1)==0)) {
        
        if (ciStamps[*niS].nss == 0) {
            refArea = (float *)calloc(sPixX*sPixY, sizeof(float));
            
            /* store the whole of the stamp in .refArea */
            cutStamp(iRData, refArea, rPixX, sXMin-rXBMin, sYMin-rYBMin,
                     sXMax-rXBMin, sYMax-rYBMin, &ciStamps[*niS]);
            
            if ( !( getStampStats3(refArea, ciStamps[*niS].x0, ciStamps[*niS].y0, sPixX, sPixY,
                                   &ciStamps[*niS].sum, &ciStamps[*niS].mean, &ciStamps[*niS].median,
                                   &ciStamps[*niS].mode, &ciStamps[*niS].sd, &ciStamps[*niS].fwhm,
                                   &ciStamps[*niS].lfwhm, 0x0, 0xffff, 3))) {
                
                /* pointless */
                /* buildSigMask(&ciStamps[*niS], sPixX, sPixY, misRData); */
                if (verbose >= 1)
                    fprintf(stderr, "    Image xs : %4i ys : %4i  (sky,dsky = %.1f,%.1f)\n",
                            ciStamps[*niS].x, ciStamps[*niS].y, ciStamps[*niS].mode, ciStamps[*niS].fwhm);
            }
            free(refArea);
        }
    }
    
    if (!(strncmp(forceConvolve, "i", 1)==0)) {
        nss = ctStamps[*ntS].nss;
        if (getCenters) {
            /* get potential centers for the kernel fit */
            getPsfCenters(&ctStamps[*ntS], tRData, sPixX, sPixY, tUKThresh, bbitt1, bbitt2);
            if (verbose >= 1)
                fprintf(stderr, "    Tmpl     : scnt = %2i nss = %2i\n",
                        ctStamps[*ntS].sscnt, ctStamps[*ntS].nss);
            
        } else {
            /* don't increment ntS inside subroutine, BUT MAKE SURE YOU DO IT OUTSIDE! */
            if (nss < nKSStamps) {
                
                if (hardX)
                    xmax = (int)(hardX);
                else
                    xmax = sXMin + (int)(fwStamp/2);
                
                if (hardY)
                    ymax = (int)(hardY);
                else
                    ymax = sYMin + (int)(fwStamp/2);
                
                check = checkPsfCenter(tRData, xmax-ctStamps[*ntS].x0, ymax-ctStamps[*ntS].y0, sPixX, sPixY,
                                       ctStamps[*ntS].x0, ctStamps[*ntS].y0, tUKThresh,
                                       ctStamps[*ntS].mode, 1. / ctStamps[*ntS].fwhm,
                                       0, 0, bbitt1 | bbitt2 | 0xbf, bbitt1);
                
                /* its ok? */
                if (check != 0.) {
                    /* globally mask out the region around this guy */
                    for (l = ymax-hwKSStamp; l <= ymax+hwKSStamp; l++) {
                        /*yr2 = l + ctStamps[*ntS].y0;*/
                        
                        for (k = xmax-hwKSStamp; k <= xmax+hwKSStamp; k++) {
                            /*xr2 = k + ctStamps[*ntS].x0;*/
                            nr2 = l+rPixX*k;
                            
                            /*if ((k > 0) && (k < sPixX) && (l > 0) && (l < sPixY))*/
                            if (nr2 >= 0 && nr2 < rPixX*rPixY)
                                mRData[nr2] |= bbitt2;
                        }
                    }
                    
                    ctStamps[*ntS].xss[nss] = xmax;
                    ctStamps[*ntS].yss[nss] = ymax;	    
                    ctStamps[*ntS].nss += 1;
                    if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d\n", nss, xmax, ymax);
                }
            }
        }
    }
    
    if (!(strncmp(forceConvolve, "t", 1)==0)) {
        nss = ciStamps[*niS].nss;
        
        if (getCenters) {
            /* get potential centers for the kernel fit */
            getPsfCenters(&ciStamps[*niS], iRData, sPixX, sPixY, iUKThresh, bbiti1, bbiti2);
            if (verbose >= 1)
                fprintf(stderr, "    Image    : scnt = %2i nss = %2i\n",
                        ciStamps[*niS].sscnt, ciStamps[*niS].nss);
        } else {
            /* don't increment niS inside subroutine, BUT MAKE SURE YOU DO IT OUTSIDE! */
            if (nss < nKSStamps) {
                if (hardX)
                    xmax = (int)(hardX);
                else
                    xmax = sXMin + (int)(fwStamp/2);
                
                if (hardY)
                    ymax = (int)(hardY);
                else
                    ymax = sYMin + (int)(fwStamp/2);
                
                check = checkPsfCenter(iRData, xmax-ciStamps[*niS].x0, ymax-ciStamps[*niS].y0, sPixX, sPixY,
                                       ciStamps[*niS].x0, ciStamps[*niS].y0, iUKThresh,
                                       ciStamps[*niS].mode, 1. / ciStamps[*niS].fwhm,
                                       0, 0, bbiti1 | bbiti2 | 0xbf, bbiti1);
                
                /* its ok? */
                if (check != 0.) {
                    /* globally mask out the region around this guy */
                    for (l = ymax-hwKSStamp; l <= ymax+hwKSStamp; l++) {
                        /*yr2 = l + ciStamps[*niS].y0;*/
                        
                        for (k = xmax-hwKSStamp; k <= xmax+hwKSStamp; k++) {
                            /*xr2 = k + ciStamps[*niS].x0;*/
                            /*nr2 = xr2+rPixX*yr2;*/
                            nr2 = l+rPixX*k;
                            
                            /*if ((k > 0) && (k < sPixX) && (l > 0) && (l < sPixY))*/
                            if (nr2 >= 0 && nr2 < rPixX*rPixY)
                                mRData[nr2] |= bbiti2;
                        }
                    }
                    
                    ciStamps[*niS].xss[nss] = xmax;
                    ciStamps[*niS].yss[nss] = ymax;	    
                    ciStamps[*niS].nss += 1;
                    if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d\n", nss, xmax, ymax);
                    
                }
            }
        }
    }
    
    return;
}

void cutStamp(float *data, float *refArea, int dxLen, int xMin, int yMin,
              int xMax, int yMax, stamp_struct *stamp) {
    
    int i, j;
    int x, y;
    int sxLen;
    
    sxLen = xMax - xMin + 1;
    /* NOTE, use j <= yMax here, not j < yMax, include yMax point */
    for (j = yMin; j <= yMax; j++) {
        y = j - yMin;
        
        for (i = xMin; i <= xMax; i++) {
            x = i - xMin;
            
            refArea[x+y*sxLen] = data[i+j*dxLen];
            /*fprintf(stderr, "%d %d %d %d %f\n", x, y, i, j, data[i+j*dxLen]); */
        }
    }
    stamp->x0 = xMin;
    stamp->y0 = yMin;
    stamp->x  = xMin + (xMax - xMin) / 2;
    stamp->y  = yMin + (yMax - yMin) / 2;
    
    return;
}

int cutSStamp(stamp_struct *stamp, float *iData) {
    /*****************************************************
     * NOTE: The Centers here are in the regions' coords
     *    To grab them from the stamp, adjust using stamp->x0,y0
     *****************************************************/
    
    int i, j, k, nss, sscnt;
    int x, y, dy, xStamp, yStamp;
    float dpt;
    double sum = 0.;
    
    /*stamp->krefArea = (double *)realloc(stamp->krefArea, fwKSStamp*fwKSStamp*sizeof(double));*/
    dfset(stamp->krefArea, fillVal, fwKSStamp, fwKSStamp);
    
    nss    = stamp->nss;
    sscnt  = stamp->sscnt;
    xStamp = stamp->xss[sscnt] - stamp->x0;
    yStamp = stamp->yss[sscnt] - stamp->y0;
    
    /* have gone through all the good substamps, reject this stamp */
    if (sscnt >= nss) {
        if (verbose >= 2)
            fprintf(stderr, "    xs : %4i ys : %4i sig: %6.3f sscnt: %4i nss: %4i ******** REJECT stamp (out of substamps)\n",
                    stamp->x, stamp->y, stamp->chi2, sscnt, nss);
        else
            if (verbose >= 1)
                fprintf(stderr, "        Reject stamp\n");
        return 1;
    }
    /*
      fprintf(stderr, "    xs : %4i ys : %4i substamp (sscnt=%d, nss=%d) xss: %4i yss: %4i\n",
      stamp->x, stamp->y, sscnt, nss, stamp->xss[sscnt], stamp->yss[sscnt]);
    */
    if (verbose >= 1)
        fprintf(stderr, "    xss : %4i yss : %4i\n", stamp->xss[sscnt], stamp->yss[sscnt]);
    
    for (j = yStamp - hwKSStamp; j <= yStamp + hwKSStamp; j++) {
        y  = j - (yStamp - hwKSStamp);
        dy = j + stamp->y0;
        
        for (i = xStamp - hwKSStamp; i <= xStamp + hwKSStamp; i++) {
            x = i - (xStamp - hwKSStamp);
            
            k   = i+stamp->x0 + rPixX*dy;
            dpt = iData[k];
            
            stamp->krefArea[x+y*fwKSStamp] = dpt;
            sum += (mRData[k] & FLAG_INPUT_ISBAD) ? 0 : fabs(dpt);
        }
    }
    
    stamp->sum = sum;
    return 0;
}

double checkPsfCenter(float *iData, int imax, int jmax, int xLen, int yLen,
                      int sx0, int sy0,
                      double hiThresh, float sky, float invdsky,
                      int xbuffer, int ybuffer, int bbit, int bbit1) {
    
    /* note the imax and jmax are in the stamp coordinate system */
    
    int brk, l, k;
    int xr2, yr2, nr2;
    double dmax2, dpt2;
    
    /* (5). after centroiding, re-check for fill/sat values and overlap (a little inefficient) */
    /*      zero tolerance for bad pixels! */
    brk   = 0;
    
    /* since we have zero tolerance for bad pixels, the sum
       of all the pixels in this box can be compared to the
       sum of all the pixels in the other boxes.  rank on
       this!  well, sum of high sigma pixels anyways...*/
    dmax2 = 0.;
    
    for (l = jmax-hwKSStamp; l <= jmax+hwKSStamp; l++) {
        if ((l < ybuffer) || (l >= yLen-ybuffer)) 
            continue; /* continue l loop */
        
        yr2 = l + sy0;
        
        for (k = imax-hwKSStamp; k <= imax+hwKSStamp; k++) {
            if ((k < xbuffer) || (k >= xLen-xbuffer)) 
                continue; /* continue k loop */
            
            xr2 = k + sx0;
            nr2 = xr2+rPixX*yr2;
            
            if (mRData[nr2] & bbit) {
                brk = 1;
                dmax2 = 0;
                break; /* exit k loop */
            }
            
            dpt2   = iData[nr2];
            
            if (dpt2 >= hiThresh) {
                mRData[nr2] |= bbit1;
                brk = 1;
                dmax2 = 0;
                break; /* exit k loop */
            }
            
            if (( (dpt2 - sky) * invdsky) > kerFitThresh)
                dmax2 += dpt2;
        }
        if (brk == 1)
            break; /* exit l loop */
    }
    
    return dmax2;
}

int getPsfCenters(stamp_struct *stamp, float *iData, int xLen, int yLen, double hiThresh, int bbit1, int bbit2) {
    /*****************************************************
     * Find the X highest independent maxima in the stamp
     *   Subject to some stringent cuts
     *   NOTE : THERE ARE SOME HARDWIRED THINGS IN HERE
     *****************************************************/
    
    int i, j, k, l, nr, nr2, imax, jmax, xr, yr, xr2, yr2, sy0, sx0, xbuffer, ybuffer;
    double dmax, dmax2, dpt, dpt2, loPsf, floor;
    int *xloc, *yloc, *qs, pcnt, brk, bbit, fcnt;
    double *peaks;
    double sky, invdsky;
    
    float dfrac = 0.9;
    
    if (stamp->nss >= nKSStamps) {
        fprintf(stderr,"    no need for automatic substamp search...\n");
        return(0);
    }
    
    /* 0xff, but 0x40 is OK = 0xbf */
    bbit    = bbit1 | bbit2 | 0xbf;
    sky     = stamp->mode;
    invdsky = 1. / stamp->fwhm;
    
    /* the highest peak in the image, period, in case all else fails */
    /* default to center of stamp */
    sx0  = stamp->x0;
    sy0  = stamp->y0;
    
    /*
      we can go with no buffer here, but, we need to allocate extra size in the
      stamp->refArea to allow this to happen.  too much room.  removed
      refArea from stamp structure, just use actual array.
    */
    xbuffer = ybuffer = 0;
    
    /* as low as we will go */
    /* was a 2 here before version 4.1.6 */
    floor = sky + kerFitThresh * stamp->fwhm;
    
    qs    = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    xloc  = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    yloc  = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    peaks = (double *)calloc(xLen*yLen/(hwKSStamp), sizeof(double));
    
    brk  = 0;
    pcnt = 0;
    fcnt = 2 * nKSStamps; /* maximum number of stamps to find */
    
    while (pcnt < fcnt) {
        
        /* NOTE : we keep the previous iteration's matches if they were good */
        
        loPsf = sky + (hiThresh - sky) * dfrac;
        loPsf = (loPsf > floor) ? loPsf : floor;   /* last ditch attempt to get some usable pixels */
        
        /* (0). ignore near the edges */
        for (j = ybuffer; j < yLen-ybuffer; j++) {
            yr  = j + sy0;
            
            for (i = xbuffer; i < xLen-xbuffer; i++) {
                xr = i + sx0;
                nr = xr+rPixX*yr;
                
                /* (1). pixel already included in another stamp, or exceeds hiThresh */
                if (mRData[nr] & bbit)
                    continue;
                
                dpt = iData[nr];
                
                /* (2a). skip masked out value */
                if (dpt >= hiThresh) {
                    mRData[nr] |= bbit1;
                    continue;
                }
                
                /* (2b). don't want low sigma point here, sometimes thats all that passes first test */
                if (( (dpt - sky) * invdsky) < kerFitThresh)
                    continue;
                
                /* finally, a good candidate! */
                if (dpt > loPsf) {
                    
                    dmax  = dpt;
                    imax  = i;
                    jmax  = j;
                    
                    /* center this candidate on the peak flux */
                    for (l = j-hwKSStamp; l <= j+hwKSStamp; l++) {
                        yr2 = l + sy0;
                        
                        if ((l < ybuffer) || (l >= yLen-ybuffer))
                            continue; /* continue l loop */
                        
                        for (k = i-hwKSStamp; k <= i+hwKSStamp; k++) {
                            xr2 = k + sx0;
                            nr2 = xr2+rPixX*yr2;
                            
                            if ((k < xbuffer) || (k >= xLen-xbuffer))
                                continue; /* continue k loop */
                            
                            /* (3). no fill/sat values within checked region, and not overlapping other stamp */
                            if (mRData[nr2] & bbit)
                                continue; /* continue k loop */
                            
                            dpt2 = iData[nr2];
                            
                            /* (4a). no fill/sat values within checked region, and not overlapping other stamp */
                            if (dpt2 >= hiThresh) {
                                mRData[nr2] |= bbit1;
                                continue; /* continue k loop */
                            }
                            
                            /* (4b). not as strong a problem, just don't want low sigma peak */
                            if (( (dpt2 - sky) * invdsky) < kerFitThresh)
                                continue;
                            
                            /* record position and amp of good point centroid in i,j,dmax */
                            if (dpt2 > dmax) {
                                dmax = dpt2;
                                imax = k;
                                jmax = l;
                            }
                        }
                    }
                    
                    /* (5). after centroiding, re-check for fill/sat values and overlap (a little inefficient) */
                    /*      zero tolerance for bad pixels! */
                    
                    /* since we have zero tolerance for bad pixels, the sum
                       of all the pixels in this box can be compared to the
                       sum of all the pixels in the other boxes.  rank on
                       this!  well, sum of high sigma pixels anyways...*/
                    
                    dmax2 = checkPsfCenter(iData, imax, jmax, xLen, yLen, sx0, sy0,
                                           hiThresh, sky, invdsky, xbuffer, ybuffer, bbit, bbit1);
                    
                    if (dmax2 == 0.)
                        continue; /* continue i loop */
                    
                    /* made it! - a valid peak */
                    xloc[pcnt]    = imax;
                    yloc[pcnt]    = jmax;
                    peaks[pcnt++] = dmax2;
                    
                    /* globally mask out the region around this guy */
                    for (l = jmax-hwKSStamp; l <= jmax+hwKSStamp; l++) {
                        yr2 = l + sy0;
                        
                        for (k = imax-hwKSStamp; k <= imax+hwKSStamp; k++) {
                            xr2 = k + sx0;
                            nr2 = xr2+rPixX*yr2;
                            
                            if ((k > 0) && (k < xLen) && (l > 0) && (l < yLen))
                                mRData[nr2] |= bbit2;
                        }
                    }
                    
                    /* found enough stamps, get out! */
                    if (pcnt >= fcnt)
                        brk = 2;
                }
                if (brk == 2)
                    break; /* exit x loop */
            }
            if (brk == 2)
                break; /* exit y loop */
        }
        
        if (loPsf == floor)
            break;        /* have hit the floor, get out  */
        dfrac -= 0.2;
        
    }
    
    if (pcnt+stamp->nss < nKSStamps) {
        if (verbose >= 2) fprintf(stderr, "    ...only found %d good substamps by autosearch\n", pcnt);
    }
    else {
        if (verbose >= 2) fprintf(stderr, "    ...found %d good substamps by autosearch\n", pcnt);
    }
    
    /* coords are in region's pixels */
    
    /* no good full stamps */
    if (pcnt == 0) {
        /* changed in v4.1.6, don't accept center pixel, could be bad, worse than not having one! */
        if (verbose >= 2) fprintf(stderr, "    NO good pixels, skipping...\n");
        
        free(qs);
        free(xloc);
        free(yloc);
        free(peaks);
        return 1;
    }
    else {
        quick_sort (peaks, qs, pcnt);
        if (verbose >= 2) fprintf(stderr, "    Adding %d substamps found by autosearch\n", imin(pcnt, nKSStamps-stamp->nss));
        for (i = stamp->nss, j = 0; j < pcnt && i < nKSStamps; i++,j++) {
            stamp->xss[i] = xloc[qs[pcnt-j-1]] + sx0;
            stamp->yss[i] = yloc[qs[pcnt-j-1]] + sy0;
            if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d,%8.1f\n", i, stamp->xss[i], stamp->yss[i], peaks[qs[pcnt-j-1]]); 
            stamp->nss++;
        }
    }
    
    free(qs);
    free(xloc);
    free(yloc);
    free(peaks);
    return 0;
}


void getNoiseStats3(float *data, float *noise, double *nnorm, int *nncount, int umask, int smask) {
    /*
      good pixels : umask=0,      smask=0xffff
      ok   pixels : umask=0xff,   smask=0x8000
      bad  pixels : umask=0x8000, smask=0
    */
    
    double nsum=0;
    int i, n=0, mdat;
    float ddat=0, ndat=0;
    
    
    for (i = rPixX*rPixY; i--; ) {
        ddat = data[i];
        mdat = mRData[i];
        
        /*
          fprintf(stderr, "CAW %d %f %f %f %d %d %d %d\n", n, nsum, ddat, ndat, mdat,
	      ((umask > 0) && (!(mdat & umask))),
	      ((smask > 0) &&   (mdat & smask)),
	      (ddat == fillVal) || (fabs(ddat) <= ZEROVAL));
        */
        
        if (((umask > 0) && (!(mdat & umask))) ||
            ((smask > 0) &&   (mdat & smask))  ||
            (fabs(ddat) <= ZEROVAL))
            continue;
        
        ndat  = 1. / noise[i];
        
        n    += 1;
        nsum += ddat*ddat * ndat*ndat;
    }
    if (n > 1) {
        *nncount = n;
        *nnorm   = nsum / (float)n;
    }
    else {
        *nncount = n;
        *nnorm   = MAXVAL;
    }
    return;
}

int getStampStats3(float *data,
                   int x0Reg, int y0Reg, int nPixX, int nPixY, 
                   double *sum, double *mean, double *median,
                   double *mode, double *sd, double *fwhm, double *lfwhm,
                   int umask, int smask, int maxiter) {
    /*****************************************************
     * Given an input image, return stats on the pixel distribution
     *    This should be a masked distribution
     *****************************************************/
    /*
      good pixels : umask=0,      smask=0xffff
      ok   pixels : umask=0xff,   smask=0x8000
      bad  pixels : umask=0x8000, smask=0
      
      this version 3 does sigma clipping 
    */
    
    /* this came primarily from Gary Bernstein */
    
    extern int flcomp();
    
    double   bin1,binsize,maxdens,moden;
    double   sumx,sumxx,isd;
    double   lower,upper,mode_bin, rdat;
    int      mdat,npts;
    int      bins[256],i,j,xr,yr;
    int      index,imax=0,tries,ilower,iupper, repeat;
    double   ssum;
    int      goodcnt;
    int      idum;
    float    *sdat;
    double   *work;
    
    int      nstat=100;
    float    ufstat=0.9;
    float    mfstat=0.5;
    
    
    npts = nPixX * nPixY;
    if (npts < nstat)
        return 4;
    
    /* fprintf(stderr, "DOINK!  %d %d %f\n", x0Reg, y0Reg, data[0]); */
    if ( !(sdat = (float *)calloc(npts, sizeof(float))) ||
         !(work = (double *)calloc(nstat, sizeof(double)))) {
        return (1);
    }
    
    idum   = -666;  /* initialize random number generator with the devil's seed */
    tries  = 0;     /* attempts at the histogram */
    
    goodcnt = 0;
    /* pull 100 random enough values from array to estimate required bin sizes */
    /* only do as many calls as there are points in the section! */
    /* NOTE : ignore anything with fillVal or zero */
    for (i = 0; (i < nstat) && (goodcnt < npts); i++, goodcnt++) {
        xr = (int)floor(ran1(&idum)*nPixX);
        yr = (int)floor(ran1(&idum)*nPixY);
        
        /* here data is size rPixX, rPixY */
        rdat = data[xr+yr*nPixX];
        /* region is rPixX, rPixY */
        mdat = mRData[(xr+x0Reg)+(yr+y0Reg)*rPixX];
        
        if (((umask > 0) && (!(mdat & umask))) ||
            ((smask > 0) &&   (mdat & smask))  ||
            (fabs(rdat) <= ZEROVAL))
            i--;
        else
            work[i] = rdat;
    }
    qsort(work, i, sizeof(double), flcomp);
    npts = i;
    
    binsize = (work[(int)(ufstat*npts)] - work[(int)(mfstat*npts)]) / (float)nstat; /* good estimate */
    bin1    = work[(int)(mfstat*npts)] - 128. * binsize;    /* we use 256 bins */
    
    /*** DO ONLY ONCE! ***/
    goodcnt = 0;
    for (j = 0; j < nPixY; j++) {
        for (i = 0; i < nPixX; i++) {
            rdat = data[i+j*nPixX];
            mdat = mRData[(i+x0Reg)+(j+y0Reg)*rPixX];
            
            /* fprintf(stderr, "BIGT %d %d %f : %d %d %d\n", i, j, rdat, i+x0Reg, j+y0Reg, mdat); */
            /* looks like it works.  that is, the image pixel values
               i+x0Reg, j+y0Reg are the same as printed pixel values
               i,j,rdat.  mask should work similarly */
            
            if (((umask > 0) && (!(mdat & umask))) ||
                ((smask > 0) &&   (mdat & smask))  ||
                (fabs(rdat) <= ZEROVAL))
                continue;
            
            if (rdat*0.0 != 0.0) {
                mRData[(i+x0Reg)+(j+y0Reg)*rPixX] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                /* fprintf(stderr, "OUCH %d %d %f %d\n", i+x0Reg, j+y0Reg, rdat, mdat); */
                continue;
            }
            
            /* good pixels pass, so sigma clipping part here */
            sdat[goodcnt++] = rdat;
        }
    }
    if (sigma_clip(sdat, goodcnt, mean, sd, maxiter)) {
        free(sdat);
        free(work);
        return 5;
    }
    free(sdat);
    
    /* save some speed */
    isd = 1. / (*sd);
    
    /*** DO ONLY ONCE! ***/
    
    repeat = 1;
    while (repeat) {
        
        if (tries >= 5) {
            /* too many attempts here - print message and exit*/
            /* DD fprintf(stderr, "     WARNING: 5 failed iterations in getStampStats2\n"); */
            free(work);
            return 1;
        }
        
        for (i=0; i<256; bins[i++]=0);
        
        /* rezero sums if repeating */
        ssum = sumx = sumxx = 0.;
        goodcnt = 0;
        
        for (j = 0; j < nPixY; j++) {
            for (i = 0; i < nPixX; i++) {
                rdat = data[i+j*nPixX];
                mdat = mRData[(i+x0Reg)+(j+y0Reg)*rPixX];
                
                if (((umask > 0) && (!(mdat & umask))) ||
                    ((smask > 0) &&   (mdat & smask))  ||
                    (fabs(rdat) <= ZEROVAL))
                    continue;
                
                if (rdat*0.0 != 0.0) {
                    mRData[(i+x0Reg)+(j+y0Reg)*rPixX] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                    continue;
                }
                
                /* final sigma cut */
                /* reject both high and low here */
                if ((fabs(rdat - (*mean)) * isd) > statSig)
                    continue;
                
                index = floor( (rdat-bin1)/binsize ) + 1;
                index = (index < 0 ? 0 : index);
                index = (index > 255 ? 255 : index);
                
                /* NOTE : the zero index is way overweighted since it contains everything
                   below bin1.
                */
                bins[index]++;
                
                /* ssum is the sum of absolute value of pixels in the image */
                ssum += fabs(rdat);
                /*
                  This is the number of pixels we are working with, not npts which
                  includes any zeroed or masked pixels.
                */
                goodcnt += 1;
                
            }
        }
        
        /* get out of dodge! */
        if (goodcnt == 0) {
            /* DD fprintf(stderr, "     WARNING: no stampStats2 data\n"); */
            *mode = *median = work[(int)(mfstat*npts)];
            *fwhm = *lfwhm = 0.;
            free(work);
            return 2;
        }
        
        /* quit here if the bins are degenerate */
        if (binsize == 0.) {
            /* DD fprintf(stderr, "     WARNING: no variation in stampStats2 data\n"); */
            *mode = *median = work[(int)(mfstat*npts)];
            *fwhm = *lfwhm = 0.;
            free(work);
            return 3;
        }
        
        /* find the mode - find narrowest region which holds ~10% of points*/
        sumx = maxdens = 0.;
        for (ilower = iupper = 1; iupper < 255; sumx -= bins[ilower++]) {
            while ( (sumx < goodcnt/10.) && (iupper < 255) ) 
                sumx += bins[iupper++];
            
            if (sumx / (iupper-ilower) > maxdens) {
                maxdens = sumx / (iupper-ilower);
                imax = ilower;
            }
        }
        /* if it never got assigned... */
        if (imax < 0 || imax > 255)
            imax = 0;
        
        
        /* try to interpolate between bins somewhat by finding weighted
           mean of the bins in the peak*/
        /* NOTE : need <= here and above in case goodcnt is small (< 10) */
        sumxx = sumx = 0.;
        for (i = imax; (sumx < goodcnt/10.) && (i < 255); i++) {
            sumx  += bins[i];
            sumxx += i*bins[i];
        }
        mode_bin = sumxx / sumx + 0.5; /*add 0.5 to give middle of bin*/
        *mode = bin1 + binsize * (mode_bin - 1.);
        
        /*find the percentile of the mode*/
        imax = floor(mode_bin);
        for (i = 0, sumx = 0.; i < imax; sumx += bins[i++]);
        sumx += bins[imax] * (mode_bin-imax); /*interpolate fractional bin*/
        sumx /= goodcnt;
        moden=sumx;
        
        /* find the region around mode containing half the "noise" points,
           e.g. assume that the mode is 50th percentile of the noise */
        lower = goodcnt * 0.25;
        upper = goodcnt * 0.75; 
        sumx = 0.;
        for (i = 0; sumx < lower; sumx += bins[i++]);
        lower = i - (sumx - lower) / bins[i-1];
        for ( ; sumx < upper; sumx += bins[i++]);
        upper = i - (sumx - upper) / bins[i-1];
        
        /*now a few checks to make sure the histogram bins were chosen well*/
        if ( (lower < 1.) || (upper > 255.) ) {
            /*make the bins wider, about same center*/
            bin1 -= 128. * binsize;
            binsize *= 2.;
            tries++;
            repeat = 1;
        } else if ( (upper-lower) < 40. ) {
            /*make bins narrower for better precision*/
            binsize /= 3.;
            bin1 = *mode - 128. * binsize;
            tries++;
            repeat=1;;
        } else
            repeat = 0;
        
        /* end of re-histogramming loop */
    }
    
    *sum = ssum;
    
    /* calculate the noise sd based on this distribution width
       the numerical constant converts the width to sd based on a
       gaussian noise distribution
    */
    *fwhm = binsize * (upper - lower) / 1.35;
    
    /*find the median*/
    for (i = 0, sumx = 0; sumx < goodcnt/2.; sumx += bins[i++]);
    *median = i - (sumx - goodcnt/2.) / bins[i-1];
    
    /* lower-quartile result - scale to sigma: */
    *lfwhm = binsize * (*median - lower) * 2. / 1.35;
    
    *median = bin1 + binsize*(*median-1.);
    
    free(work);
    return 0;
}

int sigma_clip(float *data, int count, double *mean, double *stdev, int maxiter) {
    int cnt, ncnt, i, iter;
    char *smask;
    double istdev;
    float d;
    
    if (count == 0) {
        *mean  = 0;
        *stdev = MAXVAL;
        return 1;
    }
    
    smask = (char *)calloc(count, sizeof(char));
    /*for (i=0; i<count; i++) smask[i] = 0;*/
    
    cnt  = 0;
    ncnt = count;
    iter = 0;
    
    while ((ncnt != cnt) && (iter < maxiter)) {
        cnt = ncnt;
        
        *mean  = 0;
        *stdev = 0;
        for (i=0; i<count; i++) {
            if (!(smask[i])) {
                d       = data[i];
                *mean  += d;
                *stdev += d*d;
            }
        }
        
        if (ncnt > 0) 
            *mean /= ncnt;
        else {
            *mean  = 0;
            *stdev = MAXVAL;
            free(smask);
            return 2;
        }
        
        if (ncnt > 1) {
            *stdev = *stdev - ncnt * (*mean) * (*mean);
            *stdev = sqrt(*stdev / (double)(ncnt - 1));
        }
        else {
            *stdev = MAXVAL;
            free(smask);
            return 3;
        }
        
        ncnt   = 0;
        istdev = 1. / (*stdev);
        for (i=0; i<count; i++) {
            if (!(smask[i])) {
                /* reject high and low outliers */
                if ((fabs(data[i] - (*mean)) * istdev) > statSig) {
                    smask[i] = 1;
                }
                else {
                    ncnt++;
                }
            }
        }
        iter += 1;
        /*fprintf(stderr, "%d %d %f %f\n", cnt, ncnt, (*mean), (*stdev));*/
    }
    /*fprintf(stderr, "\n");*/
    free(smask);
    return 0;
}

float *calculateAvgNoise(float *image, int *mask, int nx, int ny, int size, int maskval, int doavg) {
    /* if avg = 0, take stdev of pixel values (e.g in diffim) */
    /* if avg = 1, take mean of pixel values (e.g in noiseim)   */
    int i, j, ii, jj, cnt;
    float *data, *outim;
    double mean, stdev;
    
    
    if ( !(data = (float *) calloc(size * size, sizeof(float))))
        return (NULL);
    if ( !(outim = (float *) calloc(nx * ny, sizeof(float))))
        return (NULL);
    
    for (j = ny; j--; ) {
        for (i = nx; i--; ) {
            
            /* take your average! */
            cnt = 0;
            for (jj = j - size; jj <= j + size; jj++) {
                if ((jj < 0) || jj >= ny)
                    continue;
                
                for (ii = i - size; ii <= i + size; ii++) {
                    if ((ii < 0) || ii >= nx)
                        continue;
                    
                    if (!(mask[ii+jj*nx] & maskval))
                        data[cnt++] = image[ii+jj*nx];
                }
            }
            /* we'll do no clipping */
            sigma_clip(data, cnt, &mean, &stdev, 0);
            outim[i+j*nx] = doavg ? mean : stdev;
        }
    }
    free(data);
    return outim;
}


void freeStampMem(stamp_struct *stamps, int nStamps) {
    /*****************************************************
     * Free ctStamps allocation when ciStamps are used, vice versa
     *****************************************************/
    int i, j;
    if (stamps) {
        for (i = 0; i < nStamps; i++) {
            for(j = 0; j < nCompKer + nBGVectors; j++) 
                if (stamps[i].vectors[j]) free(stamps[i].vectors[j]);
            if (stamps[i].vectors) free(stamps[i].vectors);
            
            for (j = 0; j < nC; j++) 
                if (stamps[i].mat[j]) free(stamps[i].mat[j]);
            if (stamps[i].mat) free(stamps[i].mat);
            
            if (stamps[i].krefArea) free(stamps[i].krefArea);
            if (stamps[i].scprod) free(stamps[i].scprod);
            if (stamps[i].xss) free(stamps[i].xss);
            if (stamps[i].yss) free(stamps[i].yss);
        }
    }
}

float *makeNoiseImage4(float *iData, float invGain, float quad) {
    
    int    i;
    double qquad;
    float  *nData=NULL;
    
    if ( !(nData = (float *)calloc(rPixX*rPixY, sizeof(float))))
        return NULL;
    
    qquad = quad * quad;
    
    for (i = rPixX*rPixY; i--; ) 
        nData[i] = fabs(iData[i])*invGain + qquad;
    
    return nData;
}

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
         fits_read_key(kPtr, TINT,    "NGAUSS", &ngauss, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "FWKERN", &fwKernel, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "CKORDER", &kerOrder, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "BGORDER", &bgOrder, NULL, &status) )
        printError(status);
    
    deg_fixe    = (int *)realloc(deg_fixe,      ngauss*sizeof(int));
    sigma_gauss = (float *)realloc(sigma_gauss, ngauss*sizeof(float));
    
    /* this took a while to figure out! */
    photNormalize = (char *)malloc(1*sizeof(char));
    
    sprintf(hKeyword, "PHOTNORM");
    if (fits_read_key(kPtr, TSTRING, hKeyword, photNormalize, NULL, &status))
        printError(status);
    
    /* read kernel gaussian info */
    for (i = 0; i < ngauss; i++) {
        sprintf(hKeyword, "DGAUSS%d", i+1);
        if (fits_read_key(kPtr, TINT, hKeyword, &deg_fixe[i], NULL, &status))
            printError(status);
        sprintf(hKeyword, "SGAUSS%d", i+1);
        if (fits_read_key(kPtr, TFLOAT, hKeyword, &sigma_gauss[i], NULL, &status))
            printError(status);
        
        /* important! */
        sigma_gauss[i] = (1.0/(2.0*sigma_gauss[i]*sigma_gauss[i]));
    }
    
    if (fits_close_file(kPtr, &status) )
        printError(status);
    
    return;
}


void readKernel(char *kimage, int nRegion, double **tKerSol, double **iKerSol,
                int *rXMin, int *rXMax, int *rYMin, int *rYMax,
                double *meansigSubstamps, double *scatterSubstamps,
                double *meansigSubstampsF, double *scatterSubstampsF,
                double *diffrat, int *NskippedSubstamps) {
    /* read in kernel image for region */
    
    fitsfile *kPtr;
    int status = 0;
    char hKeyword[1024], hInfo[1024];
    
    /* open the input kernel image */
    if ( fits_open_file(&kPtr, kimage, 0, &status) )
        printError(status);
    
    /* grab stuff for this region */
    sprintf(hKeyword, "REGION%02d", nRegion);
    if (fits_read_key(kPtr, TSTRING, hKeyword, &hInfo, NULL, &status))
        printError(status);
    
    /* get extent of region */
    if (sscanf(hInfo, "[%d:%d,%d:%d]", rXMin, rXMax, rYMin, rYMax) != 4) {
        fprintf(stderr, "Problem with region %d (%s), exiting...\n", nRegion, hInfo);
        exit(1);
    }
    /* fits indexing starts at 1, code at 0 */
    *rXMin -= 1;
    *rXMax -= 1;
    *rYMin -= 1;
    *rYMax -= 1;
    
    /* which way to convolve */
    sprintf(hKeyword, "CONVOL%02d", nRegion);
    if (fits_read_key(kPtr, TSTRING, hKeyword, &hInfo, NULL, &status))
        printError(status);
    
    /* copy quality control stuff: mean sigma, scatter, # substamps skipped */
    sprintf(hKeyword, "SSSIG%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, meansigSubstamps, NULL, &status))
        printError(status);
    
    sprintf(hKeyword, "SSSCAT%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, scatterSubstamps, NULL, &status))
        printError(status);
    
    sprintf(hKeyword, "FSIG%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, meansigSubstampsF, NULL, &status))
        printError(status);
    
    sprintf(hKeyword, "FSCAT%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, scatterSubstampsF, NULL, &status))
        printError(status);
    
    /* sometimes does not exist */
    sprintf(hKeyword, "NSCALO%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, diffrat, NULL, &status)) {
        *diffrat = 1;
        status = 0;
    }
    
    if (strncmp(hInfo, "TEMPLATE", 8)==0) {
        forceConvolve = "t";
        *tKerSol = (double *)realloc(*tKerSol, (nCompTotal+1)*sizeof(double));
        
        fits_get_kernel_btbl(kPtr, &(*tKerSol), nRegion);
    }
    else if (strncmp(hInfo, "IMAGE", 5)==0) {
        forceConvolve = "i";
        *iKerSol = (double *)realloc(*iKerSol, (nCompTotal+1)*sizeof(double));
        
        fits_get_kernel_btbl(kPtr, &(*iKerSol), nRegion);
    }
    
    if (fits_close_file(kPtr, &status) )
        printError(status);
    
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
    
    /*
      fprintf(stderr, "OK %d, %f %f %f %f %f %f\n", nRegion, *kernelSol[0], *kernelSol[1],
      *kernelSol[192-1], *kernelSol[293-1],
      *kernelSol[294-1], *kernelSol[298-1]);
      */
    
    return;
}


void spreadMask(int *mData, int width) {
    /*****************************************************/
    /* spread mask by width size in all directions */
    /* it works, but is it too inefficient? */
    /*****************************************************/
    
    int i, j, k, l, ii, jj, w2;
    
    /* nothing to do! */
    if (width <= 0) {
        return;
    }
    w2 = width/2;
    for (j = 0; j < rPixY; j++) {
        for (i = 0; i < rPixX; i++) {
            if (mData[i+rPixX*j] & FLAG_INPUT_ISBAD) {
                for (k = -w2; k <= w2; k++) {
                    ii = i + k;
                    if (ii < 0 || ii >= rPixX)
                        continue;
                    
                    for (l = -w2; l <= w2; l++) {
                        jj = j + l;
                        if (jj < 0 || jj >= rPixY)
                            continue;
                        
                        mData[ii+rPixX*jj] |= FLAG_OK_CONV * (!(mData[ii+rPixX*jj] & FLAG_INPUT_ISBAD));
                    }
                }
            }
        }
    }
    return;
}

void makeInputMask(float *tData, float *iData, int *mData) {
    /*****************************************************/
    /* Construct input mask image, expand by fwKernel  */
    /*****************************************************/
    
    int i;
    
    for (i = rPixX*rPixY; i--; ){
        mData[i] |= (FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (tData[i] == fillVal  || iData[i] == fillVal);
        mData[i] |= (FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL)  * (tData[i] >= tUThresh || iData[i] >= iUThresh);
        mData[i] |= (FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL)  * (tData[i] <= tLThresh || iData[i] <= iLThresh);
    }
    
    spreadMask(mData, (int)(hwKernel*kfSpreadMask1));
    
    /* mask has value 0 for good pixels in the difference image */
    return;
}

int hp_fits_copy_header(fitsfile *iPtr, fitsfile *oPtr, int *status) {
#define FSTRNCMP(a,b,n)  ((a)[0]<(b)[0]?-1:(a)[0]>(b)[0]?1:strncmp((a),(b),(n)))   
    int nkeys, i;
    long naxis;
    char card[SCRLEN];
    
    if (fits_get_hdrspace(iPtr, &nkeys, NULL, status) ||
        fits_read_key_lng(iPtr, "NAXIS", &naxis, NULL, status) )
        return *status;
    
    /* copy remaining keywords, excluding NAXIS?, EXTEND, and reference COMMENT keywords */
    for (i = 4 + naxis; i <= nkeys; i++) {
        if (fits_read_record(iPtr, i, card, status))
            break;
        
        if (FSTRNCMP(card, "EXTEND  ", 8) &&
            FSTRNCMP(card, "COMMENT   FITS (Flexible Image Transport System) format is", 58) &&
            FSTRNCMP(card, "COMMENT   and Astrophysics', volume 376, page 3", 47) ) {
            if (fits_write_record(oPtr, card, status))
                return *status;
        }
    }
    return *status;
}

void hp_fits_correct_data(float *data, int npix, float bZero, float bScale, int makeShort) {
    
    float maxval=1e30, minval=-1e30;
    int   i;
    float *dptr;
    int   *mptr;
    
    /*
      BUYER BEWARE : the photometric rescaling of the kernel can take
      an innocent amount of flux and drive it higher than the allowed
      short maximum value.  we need to check for this here
      
      %%% AND WATCH OUT FOR BSCALE MADNESS... %%%
      
    */
    
    dptr = data;
    mptr = mRData;
    
    if (makeShort) {
        maxval =  32767. * bScale + bZero;
        minval = -32768. * bScale + bZero;
        
        for (i = 0; i < npix; i++, dptr++, mptr++) {
            if (*dptr > maxval) {
                *dptr  = maxval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
            else if (*dptr < minval) {
                *dptr  = minval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
        }
    }
    /* extra check for NaN, probably not necessary */
    
    /*
      dptr = data;
      for (i = 0; i < npix; i++, dptr++) {
      if (*dptr*0 != 0)
      *dptr = fillVal;
      }
    */
    
    return;
}

void hp_fits_correct_data_int(int *data, int npix, float bZero, float bScale, int makeShort) {
    
    float maxval=1e30, minval=-1e30;
    int   i;
    int  *dptr, *mptr;
    
    /*
      BUYER BEWARE : the photometric rescaling of the kernel can take
      an innocent amount of flux and drive it higher than the allowed
      short maximum value.  we need to check for this here
      
      %%% AND WATCH OUT FOR BSCALE MADNESS... %%%
      
    */
    
    dptr = data;
    mptr = mRData;
    
    if (makeShort) {
        maxval =  32767. * bScale + bZero;
        minval = -32768. * bScale + bZero;
        
        for (i = 0; i < npix; i++, dptr++) {
            if (*dptr > maxval) {
                *dptr  = maxval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
            else if (*dptr < minval) {
                *dptr  = minval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
        }
    }
    return;
}

int hp_fits_write_subset(fitsfile *fptr, long group, long naxis, long *naxes,
                         float *data, int *status, int makeShort,
                         float bZero, float bScale,
                         int fpixelX, int fpixelY, int lpixelX, int lpixelY, int xArrayLo, int yArrayLo) {
    
    int   pixX, pixY, y, x;
    float *dptr;
    long  fpixel[2], lpixel[2];
    
    hp_fits_correct_data(data, (rPixX*rPixY), bZero, bScale, makeShort);
    
    pixX = lpixelX - fpixelX + 1;
    pixY = lpixelY - fpixelY + 1;
    
    fpixel[0] = fpixelX;
    lpixel[0] = lpixelX;
    
    dptr = data;
    /* get to the first row */
    for (y = 0; y < yArrayLo; y++)
        for (x = 0; x < rPixX; x++, dptr++);
    
    for (y = 0; y < pixY; y++) {
        
        fpixel[1] = fpixelY + y;
        lpixel[1] = fpixel[1];
        
        /* get to first pixel in row to write */
        for (x = 0; x < xArrayLo; x++, dptr++);
        
        /* this works! */
        fits_write_subset_flt(fptr, group, naxis, naxes, fpixel, lpixel, dptr, status);
        
        /* clear the row, as fits_write_subset does not increment dptr */
        for (x = xArrayLo; x < rPixX; x++, dptr++);
    }
    return *status;
}

int hp_fits_write_subset_int(fitsfile *fptr, long group, long naxis, long *naxes,
                             int *data, int *status, int makeShort,
                             float bZero, float bScale,
                             int fpixelX, int fpixelY, int lpixelX, int lpixelY, int xArrayLo, int yArrayLo) {
    
    int   pixX, pixY, y, x;
    int  *dptr;
    long  fpixel[2], lpixel[2];
    
    hp_fits_correct_data_int(data, (rPixX*rPixY), bZero, bScale, makeShort);
    
    pixX = lpixelX - fpixelX + 1;
    pixY = lpixelY - fpixelY + 1;
    
    fpixel[0] = fpixelX;
    lpixel[0] = lpixelX;
    
    dptr = data;
    /* get to the first row */
    for (y = 0; y < yArrayLo; y++)
        for (x = 0; x < rPixX; x++, dptr++);
    
    for (y = 0; y < pixY; y++) {
        
        fpixel[1] = fpixelY + y;
        lpixel[1] = fpixel[1];
        
        /* get to first pixel in row to write */
        for (x = 0; x < xArrayLo; x++, dptr++);
        
        /* this works! */
        fits_write_subset_int(fptr, group, naxis, naxes, fpixel, lpixel, dptr, status);
        
        /* clear the row, as fits_write_subset does not increment dptr */
        for (x = xArrayLo; x < rPixX; x++, dptr++);
    }
    return *status;
}

void fset(float *data, double value, int nPixX, int nPixY) {
    int    i;
    float *d;
    d = data;
    for (i = nPixX*nPixY; i--; )
        *(d++) = value;
}

void dfset(double *data, double value, int nPixX, int nPixY) {
    int     i;
    double *d;
    d = data;
    for (i = nPixX*nPixY; i--; )
        *(d++) = value;
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

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(int *idum) {
    static long ix1,ix2,ix3;
    static double r[98];
    double temp;
    static int iff=0;
    int j;
    /* void nrerror(char *error_text); */
    
    if (*idum < 0 || iff == 0) {
        iff=1;
        ix1=(IC1-(*idum)) % M1;
        ix1=(IA1*ix1+IC1) % M1;
        ix2=ix1 % M2;
        ix1=(IA1*ix1+IC1) % M1;
        ix3=ix1 % M3;
        for (j=1;j<=97;j++) {
            ix1=(IA1*ix1+IC1) % M1;
            ix2=(IA2*ix2+IC2) % M2;
            r[j]=(ix1+ix2*RM2)*RM1;
        }
        *idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    /* if (j > 97 || j < 1) nrerror("RAN1: This cannot happen."); */
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

void quick_sort (double *list, int *index, int n) {
    
    int i;
    void quick_sort_1();
    
    for (i = 0; i < n; i++) index [i] = i;
    quick_sort_1 (list, index, 0, n-1);
    
    return;
}



void quick_sort_1(double *list, int *index, int left_end, int right_end) {
    int i,j,temp;
    double chosen;
    
    
    chosen = list[index[(left_end + right_end)/2]];
    i = left_end-1;
    j = right_end+1;
    
    for(;;) {
        while(list [index[++i]] < chosen);
        while(list [index[--j]] > chosen);
        if (i < j){
            temp=index [j];
            index [j] = index [i];
            index [i] = temp;}
        else if (i == j) {
            ++i; break;}
        else break;
    } 
    
    if (left_end < j)  quick_sort_1 (list, index, left_end, j);
    if (i < right_end) quick_sort_1 (list, index, i, right_end);
    
    return;
}

/**** comparison call for the qsort ****/
int flcomp(x,y)
     double *x,*y;
{
    if (*x>*y) return(1);
    else if (*x==*y) return(0);
    else return(-1);
}

int imin(int a, int b) {
    if (a < b) { return a; }
    else { return b; }
}   
int imax(int a, int b) {
    if (a < b) { return b; }
    else { return a; }
}   

