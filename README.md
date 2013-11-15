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
