#include<stdio.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<fitsio.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/*
  
Several of these subroutines appear originally in code created by
Cristophe Alard for ISIS, but have been modified and/or rewritten
for the current software package.  In particular, the construction
of the least squares matrices have been taken directly from the ISIS
code.

08/20/01 acbecker@physics.bell-labs.com

*/

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

int fillStamp(stamp_struct *stamp, float *imConv, float *imRef) {
    /*****************************************************
     * Fills stamp->vectors with convolved images, and 
     *   pixel indices multiplied by each other for background fit 
     *****************************************************/
    
    int       ren = 0;
    int       i,j,xi,yi,dx,dy,idegx,idegy,di,dj,nv,ig,nvec;
    double    ax,ay,xf,yf;
    double *im;
    float     rPixX2, rPixY2;
    
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
    if (verbose >= 1)
        fprintf(stderr, "    xs  : %4i ys  : %4i sig: %6.3f sscnt: %4i nss: %4i \n",
                stamp->x, stamp->y, stamp->chi2, stamp->sscnt, stamp->nss);
    if (stamp->sscnt >= stamp->nss) {
        /* have gone through all the good substamps, reject this stamp */
        /*if (verbose >= 2) fprintf(stderr, "    ******** REJECT stamp (out of substamps)\n");*/
        if (verbose >= 1)
            fprintf(stderr, "        Reject stamp\n");
        return 1;
    }
    
    nvec = 0;
    for (ig = 0; ig < ngauss; ig++) {
        for (idegx = 0; idegx <= deg_fixe[ig]; idegx++) {
            for (idegy = 0; idegy <= deg_fixe[ig]-idegx; idegy++) {
                
                ren = 0;
                dx = (idegx / 2) * 2 - idegx;
                dy = (idegy / 2) * 2 - idegy;
                if (dx == 0 && dy == 0 && nvec > 0)
                    ren = 1;
                
                /* fill stamp->vectors[nvec] with convolved image */
                /* image is convolved with functional form of kernel, fit later for amplitude */
                xy_conv_stamp(stamp, imConv, nvec, ren);
                ++nvec;
            }
        }
    }
    
    /* get the krefArea data */
    if (cutSStamp(stamp, imRef))
        return 1;
    
    /* fill stamp->vectors[nvec+++] with x^(bg) * y^(bg) for background fit */
    xi = stamp->xss[stamp->sscnt];
    yi = stamp->yss[stamp->sscnt];
    di = xi - hwKSStamp;
    dj = yi - hwKSStamp;
    for (i = xi - hwKSStamp; i <= xi + hwKSStamp; i++) {
        xf = (i - rPixX2) / rPixX2;
        
        for (j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
            /* fprintf(stderr, "%d %d %d %d %d %d\n", k, xi, yi,i, j, fwKSStamp); */
            yf = (j - rPixY2) / rPixY2; 
            
            ax = 1.0;
            nv = nvec;
            for (idegx = 0; idegx <= bgOrder; idegx++) {
                ay = 1.0; 
                for (idegy = 0; idegy <= bgOrder - idegx; idegy++) {
                    im = stamp->vectors[nv];
                    im[i-di+fwKSStamp*(j-dj)] = ax * ay;
                    ay *= yf;
                    ++nv;
                }
                ax *= xf;
            }
        }
    }
    
    /* build stamp->mat from stamp->vectors */
    build_matrix0(stamp);
    /* build stamp->scprod from stamp->vectors and imRef */
    build_scprod0(stamp, imRef);
    
    return 0;
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
    
    if (usePCA) {
        return kernel_vector_PCA(n, deg_x, deg_y, ig, ren);
    }
    
    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
    dx = (deg_x / 2) * 2 - deg_x;
    dy = (deg_y / 2) * 2 - deg_y;
    sum_x = sum_y = 0.0;
    *ren = 0;
    
    for (ix = 0; ix < fwKernel; ix++) {
        x            = (double)(ix - hwKernel);
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

double *kernel_vector_PCA(int n, int deg_x, int deg_y, int ig, int *ren) {
    /*****************************************************
     * Creates kernel sized entry for kernel_vec for each kernel degree 
     *   Mask of filter_x * filter_y, filter = exp(-x**2 sig) * x^deg 
     *   Subtract off kernel_vec[0] if n > 0
     * NOTE: this does not use any image
     ******************************************************/
    double    *vector=NULL,*kernel0=NULL;
    int i,j;
    
    vector = (double *)malloc(fwKernel*fwKernel*sizeof(double));
    
    for (i = 0; i < fwKernel; i++) {
        for (j = 0; j < fwKernel; j++) {
            vector[i+fwKernel*j] = PCA[n][i+fwKernel*j];
        }
    }
    
    if (n > 0)
        kernel0 = kernel_vec[0];
    
    if (n > 0) {
        for (i = 0; i < fwKernel * fwKernel; i++) {
            vector[i] -= kernel0[i];
        }
        *ren = 1;
    }
    
    return vector;
}

void xy_conv_stamp(stamp_struct *stamp, float *image, int n, int ren) {
    /*****************************************************
     * Called for each degree of convolution, ngauss by deg_gauss 
     * Each convolution is stored in stamp->vectors[n], imc here
     ******************************************************/
    int       i,j,xc,yc,xij,sub_width,xi,yi;
    double    *v0,*imc;
    
    
    if (usePCA) {
        xy_conv_stamp_PCA(stamp, image, n, ren);
        return;
    }
    
    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    imc = stamp->vectors[n];
    
    sub_width = fwKSStamp + fwKernel - 1;
    
    /* pull area to convolve out of full reference image region */
    /* convolve with y filter */
    for(i = xi - hwKSStamp - hwKernel; i <= xi + hwKSStamp + hwKernel; i++) {
        for(j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
            xij = i - xi + sub_width / 2 + sub_width * (j - yi + hwKSStamp);
            temp[xij] = 0.0;
            for(yc = -hwKernel; yc <= hwKernel; yc++) {
                temp[xij] += image[i+rPixX*(j+yc)] * filter_y[hwKernel-yc+n*fwKernel];
            }
        }
    }
    
    /* convolve with x filter */
    for(j = -hwKSStamp; j <= hwKSStamp; j++) {
        for(i = -hwKSStamp; i <= hwKSStamp;i++) {  
            xij = i + hwKSStamp + fwKSStamp * (j + hwKSStamp);
            imc[xij] = 0.0;
            for(xc = -hwKernel; xc <= hwKernel; xc++) {
                imc[xij] += temp[i+xc+sub_width/2+sub_width*(j+hwKSStamp)] * filter_x[hwKernel-xc+n*fwKernel];
            }
        }
    }
    
    if (ren) {
        v0 = stamp->vectors[0];
        for(i = 0; i < fwKSStamp * fwKSStamp; i++) imc[i] -= v0[i];
    }
    
    return;
}

void xy_conv_stamp_PCA(stamp_struct *stamp, float *image, int n, int ren) {
    /*****************************************************
     * Called for each degree of convolution, ngauss by deg_gauss 
     * Each convolution is stored in stamp->vectors[n], imc here
     ******************************************************/
    
    int       i,j,xc,yc,xij,xi,yi;
    double    *v0,*imc;
    
    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    imc = stamp->vectors[n];
    
    /* pull area to convolve out of full reference image region */
    for(j = yi - hwKSStamp; j <= yi + hwKSStamp; j++) {
        for(i = xi - hwKSStamp; i <= xi + hwKSStamp; i++) {
            xij      = i - (xi - hwKSStamp) + fwKSStamp * (j - (yi - hwKSStamp));
            imc[xij] = 0.;
            
            for(yc = -hwKernel; yc <= hwKernel; yc++) {
                for(xc = -hwKernel; xc <= hwKernel; xc++) {
                    imc[xij] += image[(i+xc)+rPixX*(j+yc)] * PCA[n][(xc+hwKernel) + fwKernel*(yc+hwKernel)];
                }
            }
        }
    }
    
    if (ren) {
        v0 = stamp->vectors[0];
        for(i = 0; i < fwKSStamp * fwKSStamp; i++) imc[i] -= v0[i];
    }
    
    return;
}

void fitKernel(stamp_struct *stamps, float *imRef, float *imConv, float *imNoise, double *kernelSol, 
               double *meansigSubstamps, double *scatterSubstamps, int *NskippedSubstamps) {
    /*****************************************************
     * Complete fit for kernel solution
     *****************************************************/
    
    double d, **matrix;
    char check;
    int i,mat_size;
    int ncomp1, ncomp2, ncomp, nbg_vec;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    mat_size   = ncomp1 * ncomp2 + nbg_vec + 1;
    
    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",
                              mat_size, ncomp2, ncomp1, nbg_vec);
    
    /* allocate fitting matrix */
    matrix = (double **)malloc((mat_size + 1)*sizeof(double *));
    for (i = 0; i <= mat_size; i++) 
        matrix[i] = (double *)malloc((mat_size + 1)*sizeof(double));
    
    /* allocate weight matrix */
    wxy = (double **)malloc(nS*sizeof(double *));
    for (i = 0; i < nS; i++)
        wxy[i] = (double *)malloc(ncomp2*sizeof(double));
    
    
    
    if (verbose>=2) fprintf(stderr, " Expanding Matrix For Full Fit\n");
    build_matrix(stamps, nS, matrix);
    build_scprod(stamps, nS, imRef, kernelSol);
    
    ludcmp(matrix, mat_size, indx, &d);
    lubksb(matrix, mat_size, indx, kernelSol);
    
    if (verbose>=2) fprintf(stderr, " Checking again\n");
    check = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps);
    
    while(check) {
        
        fprintf(stderr, "\n Re-Expanding Matrix\n");
        build_matrix(stamps, nS, matrix);
        build_scprod(stamps, nS, imRef, kernelSol);
        
        ludcmp(matrix, mat_size, indx, &d);
        lubksb(matrix, mat_size, indx, kernelSol);
        
        fprintf(stderr, " Checking again\n");          
        check = check_again(stamps, kernelSol, imConv, imRef, imNoise, meansigSubstamps, scatterSubstamps, NskippedSubstamps); 
    }
    fprintf(stderr, " Sigma clipping of bad stamps converged, kernel determined\n");
    
    for (i = 0; i <= mat_size; i++)
        free(matrix[i]);
    for (i = 0; i < nS; i++)
        free(wxy[i]);
    free(matrix);
    free(wxy);
    
    return;
}


void build_matrix0(stamp_struct *stamp) {
    /*****************************************************
     * Build least squares matrix for each stamp
     *****************************************************/
    
    int       i,j,pixStamp,k,i1,ivecbg=0;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double    **vec;
    
    ncomp1   = nCompKer;
    ncomp2   = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp    = ncomp1 * ncomp2;
    nbg_vec  = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    pixStamp = fwKSStamp * fwKSStamp;
    
    vec      = stamp->vectors;
    
    /* loop over the convolved images created by xy_conv_stamp() */
    /* each level represents ngauss and deg_gauss */
    for (i = 0; i < ncomp1; i++) {
        for (j = 0; j <= i; j++) {
            q = 0.0;
            /* integrate W_m1 and W_m2 (sum over all pixels) */
            for (k = 0; k < pixStamp; k++) 
                q += vec[i][k] * vec[j][k];
            
            /* Q from Eqn 3. in Alard */
            stamp->mat[i+1][j+1] = q;
        }
    }
    
    for (i1 = 0; i1 < ncomp1; i1++) { 
        ivecbg = ncomp1;
        
        p0 = 0.0;
        /* integrate convolved images and first order background (equals 1 everywhere!)*/
        for (k = 0; k < pixStamp; k++)
            p0 += vec[i1][k] * vec[ivecbg][k];
        stamp->mat[ncomp1+1][i1+1] = p0;
    }  
    
    /* integrate first order background with itself */
    /* NOTE : DON'T MASK K HERE - BACKGROUND! */
    for (k = 0, q = 0.0; k < pixStamp; k++)
        q += vec[ivecbg][k] * vec[ncomp1][k];
    stamp->mat[ncomp1+1][ncomp1+1] = q;
    
    return;
}

void build_scprod0(stamp_struct *stamp, float *image) {
    /*****************************************************
     * Build the right side of each stamp's least squares matrix
     *    stamp.scprod = degree of kernel fit + 1 bg term
     *****************************************************/
    
    int       xc,yc,xi,yi,i1,k;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double **vec;
    
    ncomp1  = nCompKer;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    vec = stamp->vectors;
    xi  = stamp->xss[stamp->sscnt];
    yi  = stamp->yss[stamp->sscnt];
    
    /* Do eqn 4. in Alard */
    
    /* Multiply each order's convolved image with reference image */
    for (i1 = 0; i1 < ncomp1; i1++) {
        p0 = 0.0;
        for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
            for (yc = -hwKSStamp; yc <= hwKSStamp; yc++) {
                k   = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
                p0 += vec[i1][k] * image[xc+xi+rPixX*(yc+yi)];
            }
        }
        stamp->scprod[i1+1] = p0;
    }
    
    /* Multiply first order background model with reference image */
    q = 0.0;
    for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
        for (yc = -hwKSStamp; yc <= hwKSStamp; yc++) {
            k  = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
            q += vec[ncomp1][k] * image[xc+xi+rPixX*(yc+yi)];	 
        }
    }
    stamp->scprod[ncomp1+1] = q;
    
    return;
}

double check_stamps(stamp_struct *stamps, int nS, float *imRef, float *imNoise) {
    /*****************************************************
     * Fit each stamp independently, reject significant outliers
     *    Next fit good stamps globally
     *    Returns a merit statistic, smaller for better fits
     *****************************************************/
    
    int    nComps,i,im,jm,mcnt1,mcnt2,mcnt3;
    double d,sum=0,kmean,kstdev;
    double merit1,merit2,merit3,sig1,sig2,sig3;
    float *m1,*m2,*m3,*ks;
    int    xc, yc, nks;
    
    double **matrix;
    int mat_size;
    int ncomp1, ncomp2, ncomp, nbg_vec;
    
    int ntestStamps;
    double       *testKerSol = NULL;
    stamp_struct *testStamps = NULL;
    
    /* kernel sum */
    ks  = (float *)calloc(nS, sizeof(float));
    nks = 0;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    mat_size   = ncomp1 * ncomp2 + nbg_vec + 1;
    
    if (verbose>=2) fprintf(stderr, " Mat_size0: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n"
                            , mat_size, ncomp2, ncomp1, nbg_vec);
    
    /* for inital fit */
    nComps      = nCompKer + 1;
    
    for (i = 0; i < nS; i++) {
        
        xc    = stamps[i].xss[stamps[i].sscnt];
        yc    = stamps[i].yss[stamps[i].sscnt];
        
        /* extract check_mat to solve one particular stamp */
        for (im = 1; im <= nComps; im++) {
            check_vec[im] = stamps[i].scprod[im];
            
            for (jm = 1; jm <= im; jm++) {
                check_mat[im][jm] = stamps[i].mat[im][jm];
                check_mat[jm][im] = check_mat[im][jm];
            }
        }
        
        /* fit stamp, the constant kernel coefficients end up in check_vec */
        ludcmp(check_mat,nComps,indx,&d);
        lubksb(check_mat,nComps,indx,check_vec);
        
        /* find kernel sum */
        sum = check_vec[1];
        check_stack[i] = sum;
        stamps[i].norm = sum;
        ks[nks++]      = sum;
        
        if (verbose >= 2) fprintf(stderr, "    # %d    xss: %4i yss: %4i  ksum: %f\n", i,
                                  stamps[i].xss[stamps[i].sscnt],
                                  stamps[i].yss[stamps[i].sscnt], sum);
    }
    
    sigma_clip(ks, nks, &kmean, &kstdev, 10);
    
    fprintf(stderr, "    %.1f sigma clipped mean ksum : %.3f, stdev : %.3f, n : %i\n",
            kerSigReject, kmean, kstdev, nks);
    
    /* so we need some way to reject bad stamps here in the first test,
       we decided to use kernel sum.  is there a better way?  part of
       the trick is that if some things are variable, you get different
       kernel sums, but the subtraction itself should come out ok. */
    
    /* stamps.diff : delta ksum in sigma */
    
    /* here we want to reject high sigma points on the HIGH and LOW
       side, since we want things with the same normalization */
    for (i = 0; i < nS; i++) {
        stamps[i].diff = fabs((stamps[i].norm - kmean) / kstdev);
    }
    
    /*****************************************************
     * Global fit for kernel solution
     *****************************************************/
    
    /* do only if necessary */
    if ((strncmp(forceConvolve, "b", 1)==0)) {
        
        /* allocate fitting matrix */
        matrix = (double **)calloc((mat_size + 1), sizeof(double *));
        for (i = 0; i <= mat_size; i++) 
            matrix[i] = (double *)calloc((mat_size + 1), sizeof(double));
        
        /* allocate weight matrix */
        wxy = (double **)calloc(nS, sizeof(double *));
        for (i = 0; i < nS; i++)
            wxy[i] = (double *)calloc(ncomp2, sizeof(double));
        
        /* first find out how many good stamps to allocate */
        ntestStamps = 0;
        for (i = 0; i < nS; i++)
            if (stamps[i].diff < kerSigReject) {
                ntestStamps++;
            }
            else {
                if (verbose >= 2) fprintf(stderr, "    # %d    skipping xss: %4i yss: %4i ksum: %f sigma: %f\n", i,
                                          stamps[i].xss[stamps[i].sscnt],
                                          stamps[i].yss[stamps[i].sscnt],
                                          stamps[i].norm, stamps[i].diff);
            }
        
        /* then allocate test stamp structure */
        if(!(testStamps = (stamp_struct *)calloc(ntestStamps, sizeof(stamp_struct)))) {
            printf("Cannot Allocate Test Stamp List\n"); 
            exit (1);
        }
        testKerSol = (double *)calloc((nCompTotal+1), sizeof(double));
        
        /* and point test stamp structure to good stamps */
        ntestStamps = 0;
        for (i = 0; i < nS; i++)
            if (stamps[i].diff < kerSigReject)
                testStamps[ntestStamps++] = stamps[i];
        
        /* finally do fit */
        if (verbose >= 2) fprintf(stderr, " Expanding Test Matrix For Fit\n");
        build_matrix(testStamps, ntestStamps, matrix);
        build_scprod(testStamps, ntestStamps, imRef, testKerSol);
        ludcmp(matrix, mat_size, indx, &d);
        lubksb(matrix, mat_size, indx, testKerSol);
        
        /* get the kernel sum to normalize figures of merit! */
        kmean = make_kernel(0, 0, testKerSol);
        
        /* determine figure of merit from good stamps */
        
        /* average of sum (diff**2 / value), ~variance */
        m1 = (float *)calloc(ntestStamps, sizeof(float));
        
        /* standard deviation of pixel distribution */
        m2 = (float *)calloc(ntestStamps, sizeof(float));
        
        /* noise sd based on histogram distribution width */
        m3 = (float *)calloc(ntestStamps, sizeof(float));
        
        mcnt1 = 0;
        mcnt2 = 0;
        mcnt3 = 0;
        for (i = 0; i < ntestStamps; i++) {
            
            getStampSig(&testStamps[i], testKerSol, imNoise, &sig1, &sig2, &sig3);
            
            if ((sig1 != -1) && (sig1 <= MAXVAL)) {
                m1[mcnt1++] = sig1;
            }
            if ((sig2 != -1) && (sig2 <= MAXVAL)) {
                m2[mcnt2++] = sig2;
            }
            if ((sig3 != -1) && (sig3 <= MAXVAL)) {
                m3[mcnt3++] = sig3;
            }
        }
        sigma_clip(m1, mcnt1, &merit1, &sig1, 10);
        sigma_clip(m2, mcnt2, &merit2, &sig2, 10);
        sigma_clip(m3, mcnt3, &merit3, &sig3, 10);
        
        /* normalize by kernel sum */
        merit1 /= kmean;
        merit2 /= kmean;
        merit3 /= kmean;
        
        /* clean up this mess */
        if (testKerSol)	 free(testKerSol);
        if (testStamps)	 free(testStamps);
        for (i = 0; i <= mat_size; i++)
            free(matrix[i]);
        for (i = 0; i < nS; i++)
            free(wxy[i]);
        free(matrix);
        free(wxy);
        
        free(m1);
        free(m2);
        free(m3);
        free(ks);
        
        /* average value of figures of merit across stamps */
        fprintf(stderr, "    <var_merit> = %.3f, <sd_merit> = %.3f, <hist_merit> = %.3f\n", merit1, merit2, merit3);
        
        /* return what is asked for if possible, if not use backup */
        if (strncmp(figMerit, "v", 1)==0) {
            if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt2 > 0) {
                return merit2;
            }
            else if (mcnt3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "s", 1)==0) {
            if (mcnt2 > 0) {
                return merit2;
            }
            else if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt3 > 0) {
                return merit3;
            }
            else {
                return 666;
            }
        }
        else if (strncmp(figMerit, "h", 1)==0) {      
            if (mcnt3 > 0) {
                return merit3;
            }
            else if (mcnt1 > 0) {
                return merit1;
            }
            else if (mcnt2 > 0) {
                return merit2;
            }
            else {
                return 666;
            }
        }
    }
    else
        return 0;
    
    return 0;
}

void build_matrix_new(stamp_struct *stamps, int nS, double **matrix) {
    /*****************************************************
     * Build overall matrix including spatial variations
     *****************************************************/
    
    int       mat_size,i,j,pixStamp,istamp,k,i1,i2,j1,j2,ibg,jbg,ivecbg,jj;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    **matrix0,p0,q,fx,fy;
    double    **vec;
    
    int ideg1, ideg2, xstamp, ystamp;
    double a1, a2;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    pixStamp   = fwKSStamp * fwKSStamp;
    
    mat_size = ncomp1 * ncomp2 + nbg_vec + 1;
    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",mat_size,ncomp2,ncomp1,nbg_vec);
    
    for(i = 0; i <= mat_size; i++)
        for(j = 0; j <= mat_size; j++)
            matrix[i][j] = 0.0;
    
    for(i = 0; i < nS; i++)
        for(j = 0; j < ncomp2; j++)
            wxy[i][j] = 0.0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        /* skip over any bad stamps along the way */
        while(stamps[istamp].sscnt >= stamps[istamp].nss) { 
            ++istamp; 
            if (istamp >= nS) break;
        }
        if (istamp >= nS) break;
        
        vec    = stamps[istamp].vectors;
        xstamp = stamps[istamp].xss[stamps[istamp].sscnt];
        ystamp = stamps[istamp].yss[stamps[istamp].sscnt];
        
        /* build weight function *HERE*, implicitly giving bad stamps zero weight */
        /*   because we skip them over above... */
        k  = 0;
        
        fx = (xstamp - rPixX/2) / rPixX/2; 
        fy = (ystamp - rPixY/2) / rPixY/2; 
        
        for (ideg1 = 0, a1 = 1.0; ideg1 <= kerOrder; ideg1++, a1 *= fx)
            for (ideg2 = 0, a2 = 1.0; ideg2 <= kerOrder - ideg1; ideg2++, a2 *= fy)
                wxy[istamp][k++] = a1 * a2;
        
        matrix0 = stamps[istamp].mat;
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            
            for (j = 0; j <= i; j++) {
                j1 = j / ncomp2;
                j2 = j - j1 * ncomp2;
                
                /* spatially weighted W_m1 and W_m2 integrals */
                matrix[i+2][j+2] += wxy[istamp][i2] * wxy[istamp][j2] * matrix0[i1+2][j1+2];
            }
        }
        
        matrix[1][1] += matrix0[1][1];
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            matrix[i+2][1] += wxy[istamp][i2] * matrix0[i1+2][1];
        }
        
        for (ibg = 0; ibg < nbg_vec; ibg++) {
            i = ncomp + ibg + 1;
            ivecbg = ncomp1 + ibg + 1;
            for (i1 = 1; i1 < ncomp1 + 1; i1++) { 
                p0 = 0.0;
                
                /* integrate convolved images over all order backgrounds */
                for (k = 0; k < pixStamp; k++)
                    p0 += vec[i1][k] * vec[ivecbg][k];
                
                /* spatially weighted image * background terms */
                for (i2 = 0; i2 < ncomp2; i2++) {
                    jj = (i1 - 1) * ncomp2 + i2 + 1;
                    matrix[i+1][jj+1] += p0 * wxy[istamp][i2];
                }
            }
            
            p0 = 0.0;
            for (k = 0; k < pixStamp; k++)
                p0 += vec[0][k] * vec[ivecbg][k];
            matrix[i+1][1] += p0;
            
            /* background * background */
            for (jbg = 0;jbg <= ibg; jbg++) {
                for (k = 0, q = 0.0; k < pixStamp; k++)
                    q += vec[ivecbg][k] * vec[ncomp1+jbg+1][k];
                matrix[i+1][ncomp+jbg+2] += q;
            }
        }
    }
    
    /* fill lower half of matrix */
    for (i = 0; i < mat_size; i++) {
        for (j = 0; j <= i; j++) {
            matrix[j+1][i+1] = matrix[i+1][j+1];
            /* fprintf(stderr, "matrix[%i][%i]: %lf\n", i,j,matrix[i+1][j+1]); */
        }
    }
    
    return;
}

void build_matrix(stamp_struct *stamps, int nS, double **matrix) {
    /*****************************************************
     * Build overall matrix including spatial variations
     *****************************************************/
    
    int       mat_size,i,j,pixStamp,istamp,k,i1,i2,j1,j2,ibg,jbg,ivecbg,jj;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    **matrix0,p0,q;
    double    **vec;
    float     rPixX2, rPixY2;
    
    int ideg1, ideg2, xstamp, ystamp;
    double a1, a2, fx, fy;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    pixStamp = fwKSStamp * fwKSStamp;
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
    mat_size = ncomp1 * ncomp2 + nbg_vec + 1;
    if (verbose >= 2) fprintf(stderr, " Mat_size: %i ncomp2: %i ncomp1: %i nbg_vec: %i \n",mat_size,ncomp2,ncomp1,nbg_vec);
    
    for(i = 0; i <= mat_size; i++)
        for(j = 0; j <= mat_size; j++)
            matrix[i][j] = 0.0;
    
    for(i = 0; i < nS; i++)
        for(j = 0; j < ncomp2; j++)
            wxy[i][j] = 0.0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        /* skip over any bad stamps along the way */
        while(stamps[istamp].sscnt >= stamps[istamp].nss) {
            ++istamp; 
            if (istamp >= nS) break;
        }
        if (istamp >= nS) break;
        
        vec    = stamps[istamp].vectors;
        xstamp = stamps[istamp].xss[stamps[istamp].sscnt];
        ystamp = stamps[istamp].yss[stamps[istamp].sscnt];
        /* RANGE FROM -1 to 1 */
        fx     = (xstamp - rPixX2) / rPixX2; 
        fy     = (ystamp - rPixY2) / rPixY2; 
        
        /* build weight function *HERE* */
        k  = 0;
        a1 = 1.0;
        for (ideg1 = 0; ideg1 <= kerOrder; ideg1++) {
            a2 = 1.0;
            for (ideg2 = 0; ideg2 <= kerOrder - ideg1; ideg2++) {
                wxy[istamp][k++] = a1 * a2;
                a2 *= fy;
            }
            a1 *= fx;
        }
        
        matrix0 = stamps[istamp].mat;
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            
            for (j = 0; j <= i; j++) {
                j1 = j / ncomp2;
                j2 = j - j1 * ncomp2;
                
                /* spatially weighted W_m1 and W_m2 integrals */
                matrix[i+2][j+2] += wxy[istamp][i2] * wxy[istamp][j2] * matrix0[i1+2][j1+2];
            }
        }
        
        matrix[1][1] += matrix0[1][1];
        for (i = 0; i < ncomp; i++) {
            i1 = i / ncomp2;
            i2 = i - i1 * ncomp2;
            matrix[i+2][1] += wxy[istamp][i2] * matrix0[i1+2][1];
        }
        
        for (ibg = 0; ibg < nbg_vec; ibg++) {
            i = ncomp + ibg + 1;
            ivecbg = ncomp1 + ibg + 1;
            for (i1 = 1; i1 < ncomp1 + 1; i1++) { 
                p0 = 0.0;
                
                /* integrate convolved images over all order backgrounds */
                for (k = 0; k < pixStamp; k++)
                    p0 += vec[i1][k] * vec[ivecbg][k];
                
                /* spatially weighted image * background terms */
                for (i2 = 0; i2 < ncomp2; i2++) {
                    jj = (i1 - 1) * ncomp2 + i2 + 1;
                    matrix[i+1][jj+1] += p0 * wxy[istamp][i2];
                }
            }
            
            p0 = 0.0;
            for (k = 0; k < pixStamp; k++)
                p0 += vec[0][k] * vec[ivecbg][k];
            matrix[i+1][1] += p0;
            
            /* background * background */
            for (jbg = 0;jbg <= ibg; jbg++) {
                for (k = 0, q = 0.0; k < pixStamp; k++)
                    q += vec[ivecbg][k] * vec[ncomp1+jbg+1][k];
                matrix[i+1][ncomp+jbg+2] += q;
            }
        }
    }
    
    /* fill lower half of matrix */
    for (i = 0; i < mat_size; i++) {
        for (j = 0; j <= i; j++) {
            matrix[j+1][i+1] = matrix[i+1][j+1];
            /* fprintf(stderr, "matrix[%i][%i]: %lf\n", i,j,matrix[i+1][j+1]); */
        }
    }
    
    return;
}

void build_scprod(stamp_struct *stamps, int nS, float *image, double *kernelSol) {
    /*****************************************************
     * Build the right side of the complete least squares matrix 
     *****************************************************/
    
    int       istamp,xc,yc,xi,yi,i1,i2,k,ibg,i,ii;
    int       ncomp1, ncomp2, ncomp, nbg_vec;
    double    p0,q;
    double    **vec;
    
    ncomp1  = nCompKer - 1;
    ncomp2  = ((kerOrder + 1) * (kerOrder + 2)) / 2;
    ncomp   = ncomp1 * ncomp2;
    nbg_vec = ((bgOrder + 1) * (bgOrder + 2)) / 2;
    
    
    for (i = 0; i <= ncomp + nbg_vec + 1; i++)
        kernelSol[i]=0.0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        /* skip over any bad stamps along the way */
        while(stamps[istamp].sscnt >= stamps[istamp].nss) {
            ++istamp; 
            if (istamp >= nS) break;
        }
        if (istamp >= nS) break;
        
        vec= stamps[istamp].vectors;
        xi = stamps[istamp].xss[stamps[istamp].sscnt];
        yi = stamps[istamp].yss[stamps[istamp].sscnt];
        
        p0 = stamps[istamp].scprod[1];
        kernelSol[1] += p0;
        
        /* spatially weighted convolved image * ref image */
        for (i1 = 1; i1 < ncomp1 + 1; i1++) {
            p0 = stamps[istamp].scprod[i1+1];
            for (i2 = 0; i2 < ncomp2; i2++) {
                ii = (i1-1) * ncomp2 + i2 + 1;
                /* no need for weighting here */
                kernelSol[ii+1] += p0 * wxy[istamp][i2];
            }
        }
        
        /* spatially weighted bg model convolved with ref image */
        for (ibg = 0; ibg < nbg_vec; ibg++) {
            q = 0.0;
            for (xc = -hwKSStamp; xc <= hwKSStamp; xc++) {
                for (yc = -hwKSStamp; yc <= hwKSStamp;yc++) {
                    k  = xc + hwKSStamp + fwKSStamp * (yc + hwKSStamp);
                    q += vec[ncomp1+ibg+1][k] * image[xc+xi+rPixX*(yc+yi)];
                    
                }
            }
            kernelSol[ncomp+ibg+2] += q;
        }
    }
    return;
}

void getFinalStampSig(stamp_struct *stamp, float *imDiff, float *imNoise, double *sig) {
    int i, j, idx, nsig=0;
    int xRegion2, xRegion = stamp->xss[stamp->sscnt];
    int yRegion2, yRegion = stamp->yss[stamp->sscnt];
    float idat, indat;
    
    *sig = 0;
    
    for (j = 0; j < fwKSStamp; j++) {
        yRegion2 = yRegion - hwKSStamp + j;
        
        for (i = 0; i < fwKSStamp; i++) {
            xRegion2 = xRegion - hwKSStamp + i;
            
            idx   = xRegion2+rPixX*yRegion2;
            idat  = imDiff[idx];
            indat = 1. / imNoise[idx];
            
            /* this shouldn't be the case, but just in case... */
            if (mRData[idx] & FLAG_INPUT_ISBAD)
                continue;
            
            nsig++;
            *sig += idat * idat * indat * indat;
            
        }
    }
    if (nsig > 0) 
        *sig /= nsig;
    else
        *sig = -1;
    
    return;
}

void getStampSig(stamp_struct *stamp, double *kernelSol, float *imNoise, double *sig1, double *sig2, double *sig3) {
    int i, j, idx, sscnt, nsig, xRegion, yRegion, xRegion2, yRegion2;
    double cSum, cMean, cMedian, cMode, cLfwhm;
    double *im, tdat, idat, ndat, diff, bg;
    
    /* info */
    sscnt   = stamp->sscnt;
    xRegion = stamp->xss[sscnt];
    yRegion = stamp->yss[sscnt];
    
    /* the comparison image */
    im = stamp->krefArea;
    /* background from fit */
    bg = get_background(xRegion, yRegion, kernelSol);
    /* temp contains the convolved image from fit, fwKSStamp x fwKSStamp */
    make_model(stamp, kernelSol, temp); 
    
    /* get sigma of stamp diff */
    nsig = 0;
    *sig1 = 0;
    *sig2 = 0;
    *sig3 = 0;
    
    for (j = 0; j < fwKSStamp; j++) {
        yRegion2 = yRegion - hwKSStamp + j;
        
        for (i = 0; i < fwKSStamp; i++) {
            xRegion2 = xRegion - hwKSStamp + i;
            
            idx  = i+j*fwKSStamp;
            
            tdat = temp[idx];
            idat = im[idx];
            ndat = imNoise[xRegion2+rPixX*yRegion2];
            
            diff = tdat - idat + bg;
            
            if ((mRData[xRegion2+rPixX*yRegion2] & FLAG_INPUT_ISBAD) || (fabs(idat) <= ZEROVAL)) {
                continue;
            }
            else {
                temp[idx] = diff;
            }
            
            /* check for NaN */
            if ((tdat*0.0 != 0.0) || (idat*0.0 != 0.0)) {
                mRData[xRegion2+rPixX*yRegion2] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                continue;
            }
            
            nsig++;
            *sig1 += diff * diff / ndat;
            /*fprintf(stderr, "OK %d %d : %f %f %f\n", xRegion2, yRegion2, tdat, idat, ndat);*/
        }
    }
    if (nsig > 0) {
        *sig1 /= nsig;
        if (*sig1 >= MAXVAL)
            *sig1 = -1;
    }
    else
        *sig1 = -1;
    
    
    /* don't do think unless you need to! */
    if (strncmp(figMerit, "v", 1)!=0) {
        if (getStampStats3(temp, xRegion - hwKSStamp, yRegion - hwKSStamp, fwKSStamp, fwKSStamp,
                           &cSum, &cMean, &cMedian, &cMode, sig2, sig3, &cLfwhm, 0x0, 0xffff, 5)) {
            *sig2 = -1;
            *sig3 = -1;
        }
        else if (*sig2 < 0 || *sig2 >= MAXVAL)
            *sig2 = -1;
        else if (*sig3 < 0 || *sig3 >= MAXVAL)
            *sig3 = -1;
    }
    
    return;
}


char check_again(stamp_struct *stamps, double *kernelSol, float *imConv, float *imRef, float *imNoise, double *meansigSubstamps, double *scatterSubstamps, int *NskippedSubstamps) {
    /*****************************************************
     * Check for bad stamps after the global fit - iterate if necessary
     *****************************************************/
    
    int    istamp,nss,scnt;
    double sig,mean,stdev;
    char   check;
    double sig1, sig2, sig3;
    float  *ss;
    
    ss  = (float *)calloc(nS, sizeof(float));
    nss = 0;
    
    sig = 0;
    check = 0;
    mean = stdev = 0.0;
    *NskippedSubstamps=0;
    
    for (istamp = 0; istamp < nS; istamp++) {
        
        /* if was fit with a good legit substamp */
        if (stamps[istamp].sscnt < stamps[istamp].nss) {
            
            getStampSig(&stamps[istamp], kernelSol, imNoise, &sig1, &sig2, &sig3);
            
            if ((strncmp(figMerit, "v", 1)==0 && (sig1 == -1)) ||
                (strncmp(figMerit, "s", 1)==0 && (sig2 == -1)) ||
                (strncmp(figMerit, "h", 1)==0 && (sig3 == -1))) {
                
                /* something went wrong with this one... */
                if (verbose>=2) fprintf(stderr, "\n    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i ITERATE substamp (BAD)\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        sig, stamps[istamp].sscnt, stamps[istamp].nss);
                
                stamps[istamp].sscnt++;
                fillStamp(&stamps[istamp], imConv, imRef);
                if (verbose>=2) fprintf(stderr, "\n");
                
                check = 1;
                
            } else {
                if (strncmp(figMerit, "v", 1)==0)
                    sig = sig1;
                else if (strncmp(figMerit, "s", 1)==0)
                    sig = sig2;
                else if (strncmp(figMerit, "h", 1)==0)
                    sig = sig3;
                
                if (verbose>=2) fprintf(stderr, "    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i OK\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        sig, stamps[istamp].sscnt, stamps[istamp].nss);
                
                stamps[istamp].chi2 = sig;
                ss[nss++]           = sig;
                
            }
        } else {
            (*NskippedSubstamps)++;
            if (verbose>=2) fprintf(stderr, "    xs : %4i ys : %4i skipping... \n",stamps[istamp].x, stamps[istamp].y);
        }
    }
    
    sigma_clip(ss, nss, &mean, &stdev, 10);
    fprintf(stderr, "    Mean sig: %6.3f stdev: %6.3f\n", mean, stdev);
    fprintf(stderr, "    Iterating through stamps with sig > %.3f\n", mean + kerSigReject * stdev);
    
    /* save the mean and scatter so that it can be saved in the fits header */
    (*meansigSubstamps)=mean;
    (*scatterSubstamps)=stdev;
    
    scnt = 0;
    for (istamp = 0; istamp < nS; istamp++) {
        /* if currently represented by a good substamp */
        if (stamps[istamp].sscnt < stamps[istamp].nss) {
            
            /* no fabs() here, keep good stamps kerSigReject on the low side! */
            if ((stamps[istamp].chi2 - mean) > kerSigReject * stdev) { 
                if (verbose>=2) fprintf(stderr, "\n    # %d    xss: %4i yss: %4i sig: %6.3f sscnt: %2i nss: %2i ITERATE substamp (poor sig)\n",
                                        istamp,
                                        stamps[istamp].xss[stamps[istamp].sscnt],
                                        stamps[istamp].yss[stamps[istamp].sscnt],
                                        stamps[istamp].chi2,
                                        stamps[istamp].sscnt,stamps[istamp].nss);
                
                stamps[istamp].sscnt++;
                scnt += (!(fillStamp(&stamps[istamp], imConv, imRef)));
                if (verbose>=2) fprintf(stderr, "\n");
                
                check = 1;
            }
            else
                scnt += 1;
        }
    }
    
    fprintf(stderr, "    %d out of %d stamps remain\n", scnt, nS);
    
    free(ss);
    return check;
}

void spatial_convolve(float *image, float **variance, int xSize, int ySize, double *kernelSol, float *cRdata, int *cMask) {
    /*****************************************************
     * Take image and convolve it using the kernelSol every kernel width
     *****************************************************/
    
    int       i1,j1,i2,j2,nsteps_x,nsteps_y,i,j,i0,j0,ic,jc,ik,jk,nc,ni,mbit,dovar;
    double    q, qv, kk, aks, uks;
    float     *vData=NULL;
    
    if ((*variance) == NULL)
        dovar = 0;
    else
        dovar = 1;
    
    if (dovar) {
        if ( !(vData = (float *)calloc(xSize*ySize, sizeof(float)))) {
            return;
        }
    }
    
    nsteps_x = ceil((double)(xSize)/(double)kcStep);
    nsteps_y = ceil((double)(ySize)/(double)kcStep);
    
    for (j1 = 0; j1 < nsteps_y; j1++) {
        j0 = j1 * kcStep + hwKernel;
        
        for(i1 = 0; i1 < nsteps_x; i1++) {
            i0 = i1 * kcStep + hwKernel;
            
            make_kernel(i0 + hwKernel, j0 + hwKernel, kernelSol);
            
            for (j2 = 0; j2 < kcStep; j2++) {
                j = j0 + j2;
                if ( j >= ySize - hwKernel) break;
                
                for (i2 = 0; i2 < kcStep; i2++) {
                    i = i0 + i2;
                    if (i >= xSize - hwKernel) break;
                    
                    ni = i+xSize*j;
                    qv = q = aks = uks = 0.0;
                    mbit = 0x0;
                    for (jc = j - hwKernel; jc <= j + hwKernel; jc++) {
                        jk = j - jc + hwKernel;
                        
                        for (ic = i - hwKernel; ic <= i + hwKernel; ic++) {
                            ik = i - ic + hwKernel;
                            
                            nc = ic+xSize*jc;
                            kk = kernel[ik+jk*fwKernel];
                            
                            q     += image[nc] * kk;
                            if (dovar) {
                                if (convolveVariance)
                                    qv += (*variance)[nc] * kk;
                                else
                                    qv += (*variance)[nc] * kk * kk;			   
                            }
                            
                            mbit  |= cMask[nc];
                            aks   += fabs(kk);
                            if (!(cMask[nc] & FLAG_INPUT_ISBAD)) {
                                uks += fabs(kk);
                            }
                        }
                    }
                    
                    cRdata[ni]   = q;
                    if (dovar)
                        vData[ni] = qv;
                    
                    /* mask propagation changed in 5.1.9 */
                    /* mRData[ni]  |= mbit; */ 
                    /* mRData[ni]  |= FLAG_OK_CONV      * (mbit > 0);*/
                    mRData[ni]  |= cMask[ni];
                    mRData[ni]  |= FLAG_OUTPUT_ISBAD * ((cMask[ni] & FLAG_INPUT_ISBAD) > 0);
                    
                    if (mbit) {
                        if ((uks / aks) < kerFracMask) {
                            mRData[ni] |= (FLAG_OUTPUT_ISBAD | FLAG_BAD_CONV);
                        }
                        else {
                            mRData[ni] |= FLAG_OK_CONV;
                        }
                    }
                    
                }       
            }
        }
    }
    if (dovar) {
        free(*variance);
        *variance = vData;
    }
    return;
}

double make_kernel(int xi, int yi, double *kernelSol) {
    /*****************************************************
     * Create the appropriate kernel at xi, yi, return sum
     *****************************************************/
    
    int    i1,k,ix,iy,i;
    double ax,ay,sum_kernel;
    double xf, yf;
    
    k  = 2;
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
            /* fprintf(stderr, "bg: %d %d %d %d %f %f %f\n", xi, yi, i, j, ax, ay, kernelSol[ncompBG+k-1]); */
            ay *= yf;
        }
        ax *= xf;
    }
    return background;
}

void make_model(stamp_struct *stamp, double *kernelSol, float *csModel) {
    /*****************************************************
     * Create a model of the convolved image
     *****************************************************/
    
    int       i1,k,ix,iy,i,xi,yi;
    double    ax,ay,coeff;
    double    *vector;
    float     rPixX2, rPixY2;
    double    xf, yf;
    
    rPixX2   = 0.5 * rPixX;
    rPixY2   = 0.5 * rPixY;
    
    xi = stamp->xss[stamp->sscnt];
    yi = stamp->yss[stamp->sscnt];
    
    /* RANGE FROM -1 to 1 */
    xf = (xi - 0.5 * rPixX) / (0.5 * rPixX);
    yf = (yi - 0.5 * rPixY) / (0.5 * rPixY);
    
    for (i = 0; i < fwKSStamp * fwKSStamp; i++) csModel[i] = 0.0;
    
    vector = stamp->vectors[0];
    coeff  = kernelSol[1];
    for (i = 0; i < fwKSStamp * fwKSStamp; i++) csModel[i] += coeff * vector[i];
    
    k=2;
    for (i1 = 1; i1 < nCompKer; i1++) {
        vector = stamp->vectors[i1];
        coeff  = 0.0; 
        ax     = 1.0;
        for (ix = 0; ix <= kerOrder; ix++) {
            ay = 1.0;
            for (iy = 0; iy <= kerOrder - ix; iy++) {
                coeff += kernelSol[k++] * ax * ay;
                ay *= yf;
            }
            ax *= xf;
        }
        
        for (i = 0; i < fwKSStamp * fwKSStamp; i++) {
            csModel[i] += coeff*vector[i];
        }
    }
    return;
}

int ludcmp(double **a, int n, int *indx, double *d)
#define TINY 1.0e-20;
{
    int     i,imax=0,j,k;
    double  big,dum,sum,temp2;
    double  *vv,*lvector();
    void    lnrerror();
    
    vv=(double *)malloc((n+1)*sizeof(double));
    
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp2=fabs(a[i][j])) > big) big=temp2;
        if (big == 0.) {
            fprintf(stderr," Numerical Recipies run error....");
            fprintf(stderr,"Singular matrix in routine LUDCMP\n");
            fprintf(stderr,"Goodbye ! \n");
            
            /*
              for (i = 1; i <= n; i++) {
              for (j = 1; j <= n; j++) {
              fprintf(stderr, "%d %d %f\n", i, j, a[i][j]);
              }
              }
            */
            return (1);
        }
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free(vv);
    return 0;
}

void lubksb(double **a, int n, int *indx, double  *b)
{
    int i,ii=0,ip,j;
    double  sum;
    
    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}
