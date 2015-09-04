; this set of routines integrates cross sections over a Maxwell
; velocity distribution and computes the rate coefficients

; Paul Barklem 
; written 2005
; general rewrite March 2006

function rkernel, Ecm

; Ecm (centre of mass frame) in cgs

common temp, T
common xsects, EcmQ, Q
common masses, mtarg, mpert
common debug, ip, ie
common interp, ilog
common threshold, thresh_cgs

kB = 1.38066d-16                              ; Boltzmann in ergs/K
bohrsq = 2.800028d-17                         ; Bohr radius squared in cm^2
Ry = 2.17991d-11                              ; Rydberg constant in ergs
mu = mtarg*mpert/(mtarg+mpert) * 1.66057d-24  ; reduced mass in g

; interpolate cross section 
; interpolation is linear (in normal or log space)

if ilog gt 0 then begin
  cross = 10.^interpol(alog10(Q), alog10(EcmQ), alog10(Ecm))
endif else begin
  cross = interpol(Q, EcmQ, Ecm)
endelse
ind = where(cross lt 0., nind)
if nind gt 0 then cross[ind] = 0.

ind = where(Ecm lt thresh_cgs, nind)
if nind gt 0 then cross[ind] = 0.

Emax = EcmQ[n_elements(EcmQ)-1]
Emin = EcmQ[0]
ind2 = where((Ecm gt Emax) or (Ecm lt Emin ), nind2)

if nind2 gt 0 and ie eq 0 then begin
  cross[ind2] = 0.d0
endif

if nind2 gt 0 and ip gt 0 and ie eq 1 then begin
  Extmin = min(Ecm[ind2])*13.6058/2.17991d-11
  Extmax = max(Ecm[ind2])*13.6058/2.17991d-11
  Exmin = Emin*13.6058/2.17991d-11
  Exmax = Emax*13.6058/2.17991d-11
  print, 'WARNING (RKERNEL) - extrapolating cross sections  '
  print, 'Min and Max interrogated: ',string(Extmax, '(e8.2)'), ' - ',string(Extmin, '(e8.2)'), ' eV'
  print, 'Input data spans :        ',string(Exmin, '(e8.2)'), ' - ', string(Exmax, '(e8.2)')
endif

const = sqrt(8.d0/!pi/mu) * (kB*T)^(-1.5d0)

ind = where(Ecm/kB/T lt 500.d0, nind)         ; numerical limit for exp(-arg)
kernel = 0.d0 * Ecm
if nind gt 0 then $
  kernel[ind] = cross[ind] * exp(-1.d0*Ecm[ind]/kB/T) * Ecm[ind] * const 
return, kernel
end

pro integ_rate, E, Q, mtarg, mpert, T, C, cofm=cofm, printdebug=printdebug, $
                ilin=ilin, ilog=ilog, thresh=thresh, noextrap=noextrap, intmeth=intmeth

; main routine
;
;  inputs:
;    E = array of energies       [eV]   (default: laboratory frame)
;    Q = array of cross sections [cm^2] 
;    mtarg = mass of target      [amu; 1 amu = 1.66057E-27 kg]
;    mpert = mass of pertuber    [amu]
;    T = temperatures            [K]
;
;  optional inputs:
;    /cofm = E is in centre of mass frame, instead of lab frame
;    /printdebug = turns on warning printouts
;    /ilin = use linear interpolation, (default is linear, but this kept)
;    /ilog = use logarithmic interpolation
;            NB: be VERY careful with log interp if there are zero xsections
;    /thresh = input threshold in eV (default, assume zero or lowest E point)
;              should be in centre-of-mass system
;    /noextrap = no extrapolation of cross sections
;    /intmeth = 0  : adaptive grid
;             = 1  : no adaptive grid, use input grid for integration with 5 point Newton-Coates
;                    not usually good since needs smooth function to work well
;             = 2  : no adaptive grid, use input grid for integration with trapezoidal method
; 
;  output
;    C = rate coefficients at T  [cm^3 /s]   

common temp, Tp
common xsects, Ecm, Qp
common masses, mtargp, mpertp
common debug, ip, ie
common interp, iilog
common threshold, thresh_cgs

ip = 0
if keyword_set(printdebug) then ip = 1

ie = 1
if keyword_set(noextrap) then ie = 0

ia = 0
if keyword_set(intmeth) then ia = intmeth

if keyword_set(thresh) then begin
   thresh_ev = thresh
endif else begin
   thresh_ev = E[0]
endelse

iilog = 0
if keyword_set(ilog) then iilog = 1 

Ew = E
Qw = Q

if iilog gt 0 then begin
  ind = where(Ew lt 1d-42, nind)  ; fudge zeros since log interpolation
  if nind gt 0 then Ew[ind] = 1d-42
  ind = where(Qw lt 1d-42, nind)  ; fudge zeros since log interpolation
  if nind gt 0 then Qw[ind] = 1d-42
endif

if keyword_set(cofm) then begin
  Ecm = Ew * 2.17991d-11/13.6058                        ; -> ergs
  thresh_cgs = thresh_ev * 2.17991d-11/13.6058
endif else begin
  Ecm = mtarg/(mpert+mtarg) * Ew * 2.17991d-11/13.6058  ; -> CofM frame in ergs
  thresh_cgs = mtarg/(mpert+mtarg) * thresh_ev * 2.17991d-11/13.6058
endelse

Qp = Qw                             ; pass to common
mtargp = mtarg
mpertp = mpert

Emin = 0.d0
if iilog gt 0 then Emin = 1.d-5 * 2.17991d-11/13.6058  ; 1e-5 eV is low limit
                                                       ; note, can't use 0 since log interp
Emx  = 1.d2 * 2.17991d-11/13.6058	; 100 eV	
				    
C = T * 0.d0                        ; allow arrays, but have to do 1-by-1

for i = 0, n_elements(T)-1 do begin
  Tp = T[i]
  Emax = max([1.38d-16*Tp * 1d2, Emx])
  if ie eq 0 then Emax = max(Ecm)   ; added PB 20090525
  ; new options 20090525 to not use adaptive grid - good for sharp resonances
  case ia of
    ; adaptive grid
    0: igral = qpint1d('rkernel', Emin, Emax, epsrel=1d-8, status=st)
    ; tabulated grid with 5 point Newton-Coates
    1: begin
         ind = where(Ecm ge 0.)
         igral = int_tabulated(Ecm[ind], rkernel(Ecm[ind]), /double)
       end
    ; trapezoidal
    2: begin
         ind = where(Ecm ge 0.)
         Ecmt = Ecm(ind)
         nEcmt = n_elements(Ecmt)
         igral = total(( rkernel(Ecmt(0l:nEcmt-2l)) + rkernel(Ecmt(1l:nEcmt-1l)) ) /2. * ( Ecmt(1l:nEcmt-1l) - Ecmt(0l:nEcmt-2l) ))
       end
  endcase
  C[i] = igral 
endfor

end
