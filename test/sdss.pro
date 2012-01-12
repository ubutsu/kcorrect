openr, lun, 'sdss.dat', /GET_LUN
header = STRARR(1)
READF, lun, header

rows = 10
data = FLTARR(17, rows)  
READF, lun, data

redshift = data(0,*)
lup = data(2:6,*)
lup_sig = data(7:11,*)
red = data(12:16,*)


lupcorr = lup - red

mag = fltarr(5,10)
mag_ivar = fltarr(5,10)
k_sdssfix,lupcorr,lup_sig,mag,mag_ivar

nmag = mag * 1e9
nmag_ivar = mag_ivar / 1e9 / 1e9

kcorr = fltarr(5,10)
kcorr = sdss_kcorrect(redshift,nmgy=nmag,ivar=nmag_ivar,absmag=absmag)

print, '### k-corrections: SDSS ugriz ###'
print, kcorr
print, ''

print, '### absolute magnitudes: SDSS ugriz ###'
print, absmag 
print, ''

kcorr = sdss_kcorrect(redshift,nmgy=nmag,ivar=nmag_ivar,absmag=absmag, $
                      coeffs=coeffs, band_shift=0.1, $
                      omaggies=om,oivar=oi)
print, '### k-corrections: SDSS ugriz band-shifted to z=0.1 ###'
print, kcorr
print, ''

print, '### absolute magnitudes: SDSS ugriz band-shifted to z=0.1 ###'
print, absmag
print, ''

; need to get vmatrix, lambda, and coeffs
kcorrect, om,oi, redshift, kcdm, vmatrix=vmatrix, lambda=lambda, $
  coeffs=coeffs, rmatrix=rmatrix, zvals=zvals

ngals = n_elements(redshift)

redshift0 = replicate(0.,n_elements(redshift))
distmod = lf_distmod(redshift)
distmod0 = lf_distmod(redshift0)

fl_sdss = 'sdss_'+['u','g','r','i','z']+'0.par'
fl = 'bessell_'+['U','B','V','R','I']+'.par'

k_reconstruct_maggies, coeffs, redshift0, rm0, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl
k_reconstruct_maggies, coeffs, redshift, rm1, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_sdss

kcorr = 2.5*alog10(rm0/rm1)

abmag = fltarr(n_elements(fl),ngals)

for i = 0L, n_elements(fl)-1L do $
    abmag[i,*] = -2.5*alog10(mag[i,*]) - distmod - kcorr[i,*]

print, '### absolute magnitudes: Bessell UBVRI ###'
print, abmag
print, ''

band_shift = 0.1
bs = replicate(band_shift,ngals)
k_reconstruct_maggies, coeffs, bs, rm0, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl ;, band_shift=band_shift
rm0 = rm0 / (1.+band_shift)
kcorr = 2.5*alog10(rm0/rm1)
for i = 0L, n_elements(fl)-1L do $
    abmag[i,*] = -2.5*alog10(mag[i,*]) - distmod - kcorr[i,*]

print, '### absolute magnitudes: Bessell UBVRI band-shifted to z=0.1 ###'
print, abmag
print, ''

fl_deep = 'deep_'+['B','R','I']+'.par'
k_reconstruct_maggies, coeffs, redshift, maggies, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_deep
appmag = -2.5*alog10(maggies)
print, '### apparent magnitudes: DEEP BRI ###'
print, appmag
print, ''

vega2ab = k_vega2ab(filterlist=fl_deep,/kurucz)
vegamag = fltarr(n_elements(fl_deep), ngals)
for i = 0L, n_elements(fl_deep)-1L do $
  vegamag[i,*] = appmag[i,*] - vega2ab[i]
print, '### apparent magnitudes: DEEP BRI on Vega scale ###'
print, vegamag
print, ''

redshift1 = replicate(1.0, n_elements(redshift))
distmod1 = lf_distmod(redshift1)
k_reconstruct_maggies, coeffs, redshift1, maggies, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_deep
for i = 0L, n_elements(fl_deep)-1L do $
  appmag[i,*] = -2.5*alog10(maggies[i,*]) - distmod + distmod1
print, '### apparent magnitudes: DEEP BRI if at z = 1 ###'
print, appmag
print, ''

exit
