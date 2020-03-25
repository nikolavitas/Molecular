;+
; NAME:
;
;   DIATOMICS_IRWIN
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial of the partition 
;   function for a given temperature for H2 and CO. The coefficients of 
;   the polynomial are from Irwin (1981). It also computes the energy of 
;   the rotational-vibrational bound states and the chemical equilibrium
;   conctant based on these polynomials. For more information, check the
;   header of diatomics.pro.
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
;   f = DIATOMICS_IRWIN(t, name, type = type, d0 = d0, mab = mab)
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
;
;-
;================================================================================
; DIATOMICS_IRWIN, IDL routine by Nikola Vitas is licensed under a Creative   
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data used in the computation are from Irwin (1981, 1981ApJS...45..621I). If  
; you use the data in a publication, please acknowledge it by citing the proper  
; reference.
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

FUNCTION diatomics_irwin, t, name, type = type, d0 = d0, mab = mab, fit = fit

;-------------------------------------------------------------------------------
; Constants
;-------------------------------------------------------------------------------
kkerg   = 1.3806488D-16        ; erg K^-1                 (from NIST)
kkev    = 8.6173324D-5         ; eV K^-1                  (from NIST)
m_h     = 1.67333D-24          ; g per H atom
amu     = 1.660538921D-24      ; g                        (from NIST)
ev2erg  = 6.24150934D11        ; 1 erg = 6.241509D11 eV   (from NIST)
hh      = 6.62606957D-27       ; erg s                    (from NIST)

IF NOT KEYWORD_SET(type) THEN type = 'pf'
IF NOT KEYWORD_SET(fit) THEN fit = 'i81'

theta = ALOG10(EXP(1.0))/(kkev * t)

IF name EQ 'H2' THEN BEGIN

  IF fit EQ 'i81' THEN BEGIN
    ; Coefficients A from Irwin (1981, first column of Table 8, p.356)
    coeffs = [1.67298118410D4, -1.49945289142D4, 5.74838863349D3, -1.22210505066D3, $
              1.55637569965D2, -1.18744926193D1, 5.02617615447D-1, -9.10563051348D-3]  
    lnu = 0.D0
    FOR ii = 0, 7 DO $
      lnu = lnu + coeffs[ii]*ALOG(t)^ii 
    u = EXP(lnu)
  ENDIF

  IF fit EQ 'st' THEN BEGIN
    ; Coefficients B from Irwin (1981, second column of Table 8, p.356)
    coeffs = [1.69179D0, -1.7227D0, 7.98033D-1, -1.57089D-1, $
             -5.35313D-1, 1.75818D0, -2.63895D0, 1.35708D0]
    logu = 0.D0
    FOR ii = 0, 7 DO $
      logu = logu + coeffs[ii]*(ALOG10(5040./t))^ii
    u = 10.^logu
  ENDIF
ENDIF
 
IF name EQ 'CO' THEN BEGIN

  IF fit EQ 'i81' THEN BEGIN
    ; Coefficients A from Irwin (1981, third column of Table 8, p.356)
    coeffs = [-5.05610415417D+4, -5.19600025580D+4, 2.33277267148D+4, -5.97599449706D+3, $
               9.55509531681D+2, -9.76517012179D+1, 6.22988547018D+0, -2.26856284960D-1, $
               3.61025385248D-3]
    lnu = 0.D0
    FOR ii = 0, 7 DO $
      lnu = lnu + coeffs[ii]*ALOG(t)^ii 
    u = EXP(lnu)
  ENDIF

  IF fit EQ 'st' THEN BEGIN
    ; Coefficients B from Irwin (1981, fourth column of Table 8, p.356)
    coeffs = [3.615300D+0, -1.773848D+0,  3.516181D-1, 8.620792D-2, 2.911791D-1, $
             -1.141469D+0,  2.513133D+0, -2.886502D+0, 1.238932D0]
    logu = 0.D0
    FOR ii = 0, 7 DO $
      logu = logu + coeffs[ii]*(ALOG10(5040./t))^ii
    u = 10.^logu
  ENDIF
ENDIF

IF type EQ 'pf' THEN result = u ELSE BEGIN
 
; Two other quantities we have to derive.

  IF type EQ 'eint' THEN BEGIN
    eint = 0.D0

    ; This is the derivative of u
    FOR ii = 1, N_ELEMENTS(coeffs)-1 DO $
      eint = eint + ii*coeffs[ii]*ALOG(t)^(ii-1) 

    result = eint*kkerg*t
  ENDIF

  IF type EQ 'kp' THEN BEGIN

    IF name EQ 'H2' THEN BEGIN
      ua = (PF(t, 1, data = 'g')).u1 
      uu = ua^2
    ENDIF

    IF name EQ 'CO' THEN BEGIN
      ua = (PF(t, 6, data = 'g')).u1
      ub = (PF(t, 8, data = 'g')).u1 
      uu = ua*ub
    ENDIF 

    kp = (2.*!pi*mab*m_h*kkerg*t/hh^2)^1.5 * EXP(-d0/(kkev*t)) * uu / u
    kp = kp * kkerg * t

    result = kp
  ENDIF

ENDELSE

RETURN, result
END
