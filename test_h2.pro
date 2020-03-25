PRO test_h2

RESTORE, 'catalogue_of_molecules.sav'

t = FINDGEN(100)*50 + 1000.
name = 'H2'
id = WHERE(cat.name EQ name)
d0 = cat[id].d0
mab = cat[id].mab
kk_erg = 1.380D-16
ev2erg  = 6.24150934D11        ; 1 erg = 6.241509D11 eV 
kk_ev  = 8.61734D-5 
hh      = 6.62606957D-27       ; erg s                    (from NIST)
theta = 5040.D0/t
nnuc = 2
d0 = 4.4780068
m_h     = 1.67333D-24

pf_sauval = diatomics_sauval_and_tatum(t, name, type = 'pf', d0 = d0)

kp_sauval = diatomics_sauval_and_tatum(t, name, type = 'kp', d0 = d0)
kp_irwin = diatomics_irwin(t, name, type = 'kp', d0 = d0, mab = mab)
kp_tsuji = diatomics_tsuji(t, name, type = 'kp')
kp_rossi = molecules_rossi_and_maciel(t, name, type = 'kp')

; They all match:
; plot, alog10(t), alog10(kp_sauval)
; oplot, alog10(t), alog10(kp_rossi), col = getcolor('salmon')
; oplot, alog10(t), alog10(kp_tsuji), col = getcolor('gold')  
; oplot, alog10(t), alog10(kp_irwin), col = getcolor('turquoise')

; They are all in CGS

eint_sauval = diatomics_sauval_and_tatum(t, name, type = 'eint')
eint_irwin = diatomics_irwin(t, name, type = 'eint')

; Compute first PF from Tsuji
@pf_path 
RESTORE, path+'pf_irwin.sav'

; Let's first check it we get the PF right.
pf_tsuji = (kk_erg * t) * (2.D0 * !dpi * mab *  m_h * kk_erg * t / hh^2)^1.5 * (PARTITION_IRWIN(t, 1)).u1^2 / kp_tsuji * EXP(-d0/kk_ev/t)

pf_sauval2 = (kk_erg * t) * (2.D0 * !dpi * mab * m_h * kk_erg * t / hh^2)^1.5 * (PARTITION_IRWIN(t, 1)).u1^2 / kp_sauval * EXP(-d0/kk_ev/t)

PLOT, ALOG10(pf_sauval)
OPLOT, ALOG10(pf_sauval2), col = GETCOLOR('turquoise')
OPLOT, ALOG10(pf_tsuji), col = GETCOLOR('gold')

; Let's get, for comparison, Eint from Sauval using PF and numerical derivative (and neglecting the 
; derivative of H partition function:

; This is nearly identical to the values of Eint provided by diatomics_sauval_and_tatum
eint2_sauval = - kk_erg*t* deriv(alog10(theta), alog10(pf_sauval))

; This one identical as well 
eint3_sauval = kk_erg*t*(2.5 + deriv(alog10(theta), alog10(kp_sauval))) + d0/ev2erg

eint1_tsuji = kk_erg*t*(2.5 + deriv(alog10(theta), alog10(kp_tsuji))) + d0/ev2erg


; And now compute Eint from explicit kp by Tsuji
  
index = WHERE(tsuji_eqc.name EQ 'H2')
coeffs = tsuji_eqc[index].coeffs
 
dlogkp = 0.D0
FOR ii = 1, N_ELEMENTS(coeffs)-1 DO dlogkp += coeffs[ii]*ii*theta^(ii-1)
; dkp = 10.D0^dkp

eint2_tsuji = kk_erg*t*(2.5 + theta*ALOG(10.)*dlogkp) + d0/ev2erg


eint3_tsuji = kk_erg*t*(2.5 + 2*deriv(t, ALOG((PARTITION_IRWIN(t, 1)).u1))  + (theta *ALOG(10.)) * dlogkp) + d0/ev2erg


; They all match as well (Tsuji somewhat different from Sauva, as expected, Irwin is ringing)      
plot, alog10(t), alog10(eint2_sauval)
oplot, alog10(t), alog10(eint_sauval), col = getcolor('gold')
oplot, alog10(t), alog10(eint3_sauval), col = getcolor('salmon')                    
oplot, alog10(t), alog10(eint1_tsuji), col = getcolor('magenta') 
oplot, alog10(t), alog10(eint2_tsuji), col = getcolor('chocolate'), psym = -1
oplot, alog10(t), alog10(eint3_tsuji), col = getcolor('turquoise')
oplot, alog10(t), alog10(eint_irwin), col = getcolor('red')







END
