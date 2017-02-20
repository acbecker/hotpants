hotpants
========

Import of v5.1.11 of High Order Transform of Psf ANd Template Subtraction code (hotpants).

Note on usage: Your mileage will vary based on the configuration of the software.  The most important tuning parameter is the size of the gaussians that you use.  A good rule of thumb is, asssuming you have measured the widths of the Psfs in the science and template image:

 * Sigma_image < Sigma_template : This requires deconvolution (sharpening) of the template.  This will lead to false positives, in practice.  Consider convolving the science image instead (-c i).  OR, since you really don't want to mess with the science pixels unnecessarily, consider convolving the science image with its Psf *before* matching the template to it.  This process is typically done after image subtraction for optimal point source filtering; in this case, the image should not be convolved with anything before detection, or just convolved with a delta function.  I.e.

   * Difference Image: D = I - T x K
   * Detect on difference image: D' = D x PSF = I x PSF - T x K x PSF
   * Instead, prefilter with Psf: I' = I x PSF
   *                              D' = I' - T x K'
   * Ideally K' = K x PSF
   * This effectively makes the image you match T to (I') have a larger PSF by sqrt(2) compared to I, avoiding deconvolution in many cases.

   

 * Sigma_image > Sigma_template : This leads to smoothing of the template.  Assume that both Psfs are Gaussian, in which case the Gaussian that matches the two has Sigma_match = sqrt(Sigma_image**2 - Sigma_template**2).  It is recommended that this be the central Gaussian in your kernel basis, with the smallest one being 0.5 * Sigma_match and the largest being 2.0 * Sigma_match.  Set these using the -ng flag.  E.g. -ng 3 6 0.5*Sigma_match 4 Sigma_match 2 2.0*Sigma_match.

######### All command line options
```
Version 5.1.11
Required options:
   [-inim fitsfile]  : comparison image to be differenced
   [-tmplim fitsfile]: template image
   [-outim fitsfile] : output difference image

Additional options:
   [-tu tuthresh]    : upper valid data count, template (25000)
   [-tuk tucthresh]  : upper valid data count for kernel, template (tuthresh)
   [-tl tlthresh]    : lower valid data count, template (0)
   [-tg tgain]       : gain in template (1)
   [-tr trdnoise]    : e- readnoise in template (0)
   [-tp tpedestal]   : ADU pedestal in template (0)
   [-tni fitsfile]   : input template noise array (undef)
   [-tmi fitsfile]   : input template mask image (undef)
   [-iu iuthresh]    : upper valid data count, image (25000)
   [-iuk iucthresh]  : upper valid data count for kernel, image (iuthresh)
   [-il ilthresh]    : lower valid data count, image (0)
   [-ig igain]       : gain in image (1)
   [-ir irdnoise]    : e- readnoise in image (0)
   [-ip ipedestal]   : ADU pedestal in image (0)
   [-ini fitsfile]   : input image noise array (undef)
   [-imi fitsfile]   : input image mask image (undef)

   [-ki fitsfile]    : use kernel table in image header (undef)
   [-r rkernel]      : convolution kernel half width (10)
   [-kcs step]       : size of step for spatial convolution (2 * rkernel + 1)
   [-ft fitthresh]   : RMS threshold for good centroid in kernel fit (20.0)
   [-sft scale]      : scale fitthresh by this fraction if... (0.5)
   [-nft fraction]   : this fraction of stamps are not filled (0.1)
   [-mins spread]    : Fraction of kernel half width to spread input mask (1.0)
   [-mous spread]    : Ditto output mask, negative = no diffim masking (1.0)
   [-omi  fitsfile]  : Output bad pixel mask (undef)
   [-gd xmin xmax ymin ymax]
                     : only use subsection of full image (full image)

   [-nrx xregion]    : number of image regions in x dimension (1)
   [-nry yregion]    : number of image regions in y dimension (1)
   -- OR --
   [-rf regionfile]  : ascii file with image regions 'xmin:xmax,ymin:ymax'
   -- OR --
   [-rkw keyword num]: header 'keyword[0->(num-1)]' indicates valid regions

   [-nsx xstamp]     : number of each region's stamps in x dimension (10)
   [-nsy ystamp]     : number of each region's stamps in y dimension (10)
   -- OR --
   [-ssf stampfile]  : ascii file indicating substamp centers 'x y'
   -- OR --
   [-cmp cmpfile]    : .cmp file indicating substamp centers 'x y'

   [-afssc find]     : autofind stamp centers so #=-nss when -ssf,-cmp (1)
   [-nss substamps]  : number of centroids to use for each stamp (3)
   [-rss radius]     : half width substamp to extract around each centroid (15)

   [-savexy file]    : save positions of stamps for convolution kernel (undef)
   [-c  toconvolve]  : force convolution on (t)emplate or (i)mage (undef)
   [-n  normalize]   : normalize to (t)emplate, (i)mage, or (u)nconvolved (t)
   [-fom figmerit]   : (v)ariance, (s)igma or (h)istogram convolution merit (v)
   [-sconv]          : all regions convolved in same direction (0)
   [-ko kernelorder] : spatial order of kernel variation within region (2)
   [-bgo bgorder]    : spatial order of background variation within region (1)
   [-ssig statsig]   : threshold for sigma clipping statistics  (3.0)
   [-ks badkernelsig]: high sigma rejection for bad stamps in kernel fit (2.0)
   [-kfm kerfracmask]: fraction of abs(kernel) sum for ok pixel (0.990)
   [-okn]            : rescale noise for 'ok' pixels (0)
   [-fi fill]        : value for invalid (bad) pixels (1.0e-30)
   [-fin fill]       : noise image only fillvalue (0.0e+00)
   [-convvar]        : convolve variance not noise (0)

   [-oni fitsfile]   : output noise image (undef)
   [-ond fitsfile]   : output noise scaled difference image (undef)
   [-nim]            : add noise image as layer to sub image (0)
   [-ndm]            : add noise-scaled sub image as layer to sub image (0)

   [-oci fitsfile]   : output convolved image (undef)
   [-cim]            : add convolved image as layer to sub image (0)

   [-allm]           : output all possible image layers

   [-nc]             : do not clobber output image (0)
   [-hki]            : print extensive kernel info to output image header (0)

   [-oki fitsfile]   : new fitsfile with kernel info (under)

   [-sht]            : output images 16 bitpix int, vs -32 bitpix float (0)
   [-obs bscale]     : if -sht, output image BSCALE, overrides -inim (1.0)
   [-obz bzero]      : if -sht, output image BZERO , overrides -inim (0.0)
   [-nsht]           : output noise image 16 bitpix int, vs -32 bitpix float (0)
   [-nbs bscale]     : noise image only BSCALE, overrides -obs (1.0)
   [-nbz bzero]      : noise image only BZERO,  overrides -obz (0.0)

   [-ng  ngauss degree0 sigma0 .. degreeN sigmaN]
                     : ngauss = number of gaussians which compose kernel (3)
                     : degree = degree of polynomial associated with gaussian #
                                (6 4 2)
                     : sigma  = width of gaussian #
                                (0.70 1.50 3.00)
                     : N = 0 .. ngauss - 1

                     : (3 6 0.70 4 1.50 2 3.00
   [-pca nk k0.fits ... n(k-1).fits]
                     : nk      = number of input basis functions
                     : k?.fits = name of fitsfile holding basis function
                     : Since this uses input basis functions, it will fix :
                     :    hwKernel 
                     :    
   [-v] verbosity    : level of verbosity, 0-2 (1)
 NOTE: Fits header params will be added to the difference image
       COMMAND             (what was called on the command line)
       NREGION             (number of regions in image)
       PHOTNORM            (to which system the difference image is normalized)
       TARGET              (image which was differenced)
       TEMPLATE            (template for the difference imaging)
       DIFFIM              (output difference image)
       MASKVAL             (value for masked pixels)
       REGION??            (IRAF-format limits for each region in the image)
       CONVOL??            (which image was convolved for each region)
       KSUM??              (sum of the convolution kernel for each region)
```
