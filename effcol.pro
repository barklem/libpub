pro effcol, T, C, g, dE, colstr, down=down
;
; computes effective collision strengths for an array of rate coefficients 
;
;  T = array of temperatures
;  C = array of rate coefficients in cgs units, cm^3/s
;  g = the statistical weight of the initial state
;  dE = energy of transition, should be positive, in eV
;
;  /down -> is rate coefficients are for downward process
; see Hummer et al 1993

k = 8.6173324d-5  ; Boltzmann ev/K
Ac = 8.629e-6 ; cm^3 s^-1 K^1/2
if keyword_set(down) then begin
   colstr = C/Ac *g *sqrt(T) 
endif else begin
   ;colstr = C/Ac *g *sqrt(T) * exp(dE/(k*T)) 
   colstr = exp(alog(C/Ac *g *sqrt(T))+dE/(k*T))  ; limits numerical error at small T
end
end