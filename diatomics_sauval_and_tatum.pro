;+
; NAME:
;
;   DIATOMICS_SAUVAL_AND_TATUM
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial of the partition 
;   function and the chemical equilibrium constant for a given temperature  
;   and diatomic specie. The coefficients of the polynomial are from Sauval
;   and Tatum (1984). It also computed the energy of the rotational-vibrational 
;   bound states based on these polynomials. For more information, check
;   the header of diatomics.pro.
;
; AUTHOR:
;
;       Nikola Vitas
;       Instituto de Astrofisica de Canarias (IAC)
;       C/ Via Lactea, s/n
;       E38205 - La Laguna (Tenerife), Espana
;       Email: n.vitas@iac.es, nikola.vitas@gmail.com
;       Homepage: nikolavitas.blogspot.com
;
; CATEGORY:
;
;   Atomic data.
;
; CALLING SEQUENCE:
;
;   f = DIATOMICS_SAUVAL_AND_TATUM(t, name, type = type, d0 = d0, mab = mab)
;
; INPUTS:
;
;   t =      Scalar or array, float. Temperature (K)
;
;   name =   Scalar string. Name = chemical formula of the molecule
;
; OUTPUTS:
; 
;   f  =     Array, float. The output function (the partition function,
;            the chemical equilibrium constant or the internal energy).
;
; KEYWORDS:
;
;   type = String. Specifies the output function ('pf', the partition function,
;          default; 'kp', the chemical equilibrium constant for pressure;
;          or 'eint', the internal energy).
;
;   d0   = Scalar, float. The dissociation constant in eV.
;
;   mab  = Scalar, float. The reduced mass, mab = ma*mb/(ma+mb) in a.m.u.
;   
; COMMENT:
;
; DEPENDENCIES:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;   Written by: Nikola Vitas (February 2014)
;   Updated to match the new molecules routine (Apr 2019)
;
;-
;================================================================================
; DIATOMICS_SAUVAL_AND_TATUM, IDL routine by Nikola Vitas is licensed under a  
; Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in sauvalandtatum1984_*.sav are not integral part of this routine. If  
; you use the data in a publication, please acknowledge it by citing the proper  
; reference (1984ApJS...56..193S).
;
; This software is provided by NV ''as is'' and any express or implied warranties, 
; including, but not limited to, the implied warranties of merchantability and 
; fitness for a particular purpose are disclaimed. In no event shall NV be liable 
; for any direct, indirect, incidental, special, exemplary, or consequential 
; damages (including, but not limited to, procurement of substitute goods or 
; services; loss of use, data, or profits; loss of use, data, or profits; or 
; business interruption) however caused and on any theory of liability, whether 
; in contract, strict liability, or tort (including negligence or otherwise) 
; arising in any way out of the use of this software, even if advised of the 
; possibility of such damage.
;================================================================================

FUNCTION diatomics_sauval_and_tatum, t, name, type = type, d0 = d0, mab = mab

IF NOT KEYWORD_SET(type) THEN type = 'pf'

kkev    = 8.6173324D-5         ; eV K^-1                  (from NIST)
kkerg   = 1.3806488D-16        ; erg K^-1                 (from NIST)
theta = ALOG10(EXP(1.0))/(kkev * t)

@molecules_path

IF type EQ 'pf' OR type EQ 'eint' THEN BEGIN

  RESTORE, path + 'catalogue_of_molecules.sav'
  data = st_pf
  index = WHERE(name EQ data.name)

  u = t*0.D0
  coeffs = data[index].coeffs
  nc = N_ELEMENTS(coeffs)

  logu = 0.D0
  FOR ii = 0, nc-1 DO $
    logu = logu + coeffs[ii]*(ALOG10(theta))^ii
  u = 10.D0^logu

  result = u

ENDIF

IF type EQ 'kp' THEN BEGIN

  RESTORE, path + 'catalogue_of_molecules.sav'
  data = st_eqc
  index = WHERE(name EQ data.name)

  kp = t*0.D0
  coeffs = data[index].coeffs
  nc = N_ELEMENTS(coeffs)

  FOR ii = 0, nc-1 DO $
    kp = kp + coeffs[ii]*(ALOG10(theta))^ii
  kp = kp - d0*theta + 1.D0

  kp = 10.D0^kp

; kkerg   = 1.3806488D-16        ; erg K^-1                 (from NIST)
; m_h     = 1.67333D-24          ; g per H atom
; amu     = 1.660538921D-24      ; g                        (from NIST)
; ev2erg  = 6.24150934D11        ; 1 erg = 6.241509D11 eV   (from NIST)
; hh      = 6.62606957D-27       ; erg s                    (from NIST)
; hm_ion1 = 0.754209             ; eV                       (Pekeris, Phys.Rev. 1958, 1962)
;    Alternative (gives identical result)
;    kp_st2 = (2.*!pi*mab*m_h*kkerg*t/hh^2)^1.5 * EXP(-d0/(kkev*t)) * uh.u1^2 / ust   
;    kp_st2 = kp_st2 * kkerg * t

  result = kp

ENDIF

IF type EQ 'eint' THEN BEGIN
  ; eint = kT d ln U / d ln T = kT d log U / d log T
  ; d log U / d log T = d log U / d log theta * d log theta / d log T =
  ;   =  (d log U / d log theta) * (T/theta) d theta / d T = 
  ;   =  - (d log U / d log theta)

  eint = 0.D0
  FOR ii = 1, nc-1 DO $
    eint = eint - coeffs[ii]*ii*(ALOG10(theta))^(ii-1)
  
  eint = eint * kkerg * t

  result = eint
ENDIF

RETURN, result
END
