; testing output from the kcorrect (v4_1_4) web site examples

kcorrect, [1., 4.78, 10.96, 14.45, 19.05],  $
  [1100., 28., 7.7, 4.4, 2.5], $
  0.03, kcorrect, band_shift=0.1, chi2=chi2, $
  rmatrix=rmatrix, zvals=zvals, /silent
print, kcorrect

kcorrect, [1., 4.73, 11.26, 14.25, 18.85],  $
  [1100., 28., 7.7, 4.4, 2.5], $
  0.03, kcorrect, band_shift=0.1, chi2=chi2, $
  rmatrix=rmatrix, zvals=zvals, /silent, vname='default'
print, kcorrect

kcorrect, [1., 4.78, 10.96, 14.45, 19.05],  $
          [1100., 28., 7.7, 4.4, 2.5], $
          0.03, kcorrect, band_shift=0.1, chi2=chi2, $
          vmatrix=vmatrix, lambda=lambda, coeffs=coeffs, /silent
plot, lambda, vmatrix#coeffs, xra=[2000., 12000.]

k_reconstruct_maggies, coeffs, 0.03, maggies, vmatrix=vmatrix, $
    lambda=lambda, filterlist=['sdss_u0.par', 'sdss_g0.par', $
    'sdss_r0.par', 'sdss_i0.par', 'sdss_z0.par'], /silent
print, maggies

k_reconstruct_maggies, coeffs, 0.03, maggies, vmatrix=vmatrix, $
    lambda=lambda, filterlist=['bessell_B.par', 'bessell_V.par'], /silent
vega2ab=k_vega2ab(filterlist=['bessell_B.par', 'bessell_V.par'], /kurucz)
bessellmags=-2.5*alog10(maggies)-vega2ab
BminusV=bessellmags[0]-bessellmags[1]

print, BminusV

kphotoz, [1., 4.78, 10.96, 14.45, 19.05],  $
         [1100., 28., 7.7, 4.4, 2.5], $
         photoz, rmatrix=rmatrix, zvals=zvals
print, photoz

exit
