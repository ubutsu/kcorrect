openr, lun, 'deep.dat', /GET_LUN
header = STRARR(1)
READF, lun, header

rows = 10
data = FLTARR(5, rows)  
READF, lun, data

redshift = data(0,*)
mag = data(2:4,*)

; prepare zcat like structure
zcatel = create_struct('magb',0.,'magr',0.,'magi',0.)
zcat = replicate(zcatel,10)
zcat.magb = reform(mag[0,*])
zcat.magr = reform(mag[1,*])
zcat.magi = reform(mag[2,*])

kcorr_nub = fltarr(3,10)
kcorr_nub = deep_kcorrect(redshift,zcat=zcat,absmag=absmag_nub, $
                          omaggies=om,oivar=oi)


fl_deep = 'deep_'+['B','R','I']+'.par'
ngals = n_elements(redshift)

kcorrect, om,oi, redshift, kcorr_deep, vmatrix=vmatrix, lambda=lambda, $
  coeffs=coeffs, rmatrix=rmatrix, zvals=zvals, filterlist=fl_deep, $
  absmag = absmag_deep

redshift0 = replicate(0.,n_elements(redshift))
distmod = lf_distmod(redshift)
distmod0 = lf_distmod(redshift0)

print, '### k-corrections: DEEP BRI ###'
print, kcorr_deep
print, ''

print, '### absolute magnitudes: DEEP BRI ###'
print, absmag_deep
print, ''

band_shift = 1.0

kcorrect, om,oi, redshift, kcorr_deep, $
  filterlist=fl_deep, $
  absmag = absmag_deep, band_shift = band_shift

print, '### k-corrections: DEEP BRI band-shifted to z=1 ###'
print, kcorr_deep
print, ''

print, '### absolute magnitudes: DEEP BRI band-shifted to z=1 ###'
print, absmag_deep
print, ''

print, '### k-corrections: NUB (IDL default) ###'
print, kcorr_nub
print, ''

print, '### absolute magnitudes: NUB (IDL default) ###'
print, absmag_nub
print, ''


fl_bessell = 'bessell_'+['U','B','V','R','I']+'.par'
redshift0 = replicate(0.,n_elements(redshift))
distmod = lf_distmod(redshift)
distmod0 = lf_distmod(redshift0)

appmag = fltarr(n_elements(fl_bessell), ngals)
k_reconstruct_maggies, coeffs, redshift0, rm0, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_bessell
k_reconstruct_maggies, coeffs, redshift, rm1, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_bessell

kcorr = 2.5*alog10(rm0/rm1)

appmag = -2.5*alog10(rm0)

vega2ab = k_vega2ab(filterlist=fl_bessell,/kurucz)
vegamag = fltarr(n_elements(fl_bessell), ngals)
for i = 0L, n_elements(fl_bessell)-1L do $
  vegamag[i,*] = appmag[i,*] -vega2ab[i] - distmod; - kcorr[i]
print, '### absolute magnitudes: Bessell UBVRI on Vega scale ###'
print, vegamag
print, ''

fl_sdss = 'sdss_'+['u','g','r','i','z']+'0.par'

redshift01 = replicate(0.1,ngals)
distmod01 = lf_distmod(redshift01)
k_reconstruct_maggies, coeffs, redshift01, rm, $
  vmatrix=vmatrix, lambda=lambda, filterlist=fl_sdss
appmag = -2.5*alog10(rm)

for i = 0L, n_elements(fl_bessell)-1L do $
  appmag[i,*] = appmag[i,*] - distmod + distmod01

print, '### apparent magnitude SDSS ugriz if at z = 0.1 ###'
print, appmag
print, ''










exit



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
