;+
; NAME:
;
;   DIATOMICS_BOHN
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial of the partition 
;   function and for a given temperature and for H2 and CO. The coefficients
;   of the polynomial are from Bohn and Wolf (1984). It also computed the energy
;   of the rotational-vibrational bound states and the chemical equilibrium 
;   constant based on these polynomials. For more information, check the header
;   of diatomics.pro.
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
;   f = DIATOMICS_BOHN(t, name, type = type, d0 = d0, mab = mab)
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
; DIATOMICS_BOHN, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in used in this routine comes from Bohn and Wolf (1984, 
; 1984A&A...130..202B). If you use the data in a  publication, please acknowledge 
; it by citing that references.

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
FUNCTION diatomics_bohn, t, name, type = type, d0 = d0, mab = mab

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

theta = ALOG10(EXP(1.0))/(kkev * t)

; This coefficients are from Table 1. (Bohn and Wolf, 1984, p.204)
 IF name EQ 'H2' THEN BEGIN

  a0 = [9.5547D-1, 1.0395D-4, -7.1985D-8, 1.1544D-11, -6.4035D-16]  
  b0 = [5.7971D-3, 4.6761D-8, 1.0886D-11]

ENDIF
 
IF name EQ 'CO' THEN BEGIN

  a0 = [1.1377D-1, -1.7072D-4, -1.2576D-8, 6.5767D-12, -5.2274D-16]  
  b0 = [3.5987D-1, 2.3872D-6, 5.1280D-11]

ENDIF

u = t * (b0[0] + b0[1]*t + b0[2]*t^2)
u = u / (a0[0] + a0[1]*t + a0[2]*t^2 + a0[3]*t^3 + a0[4]*t^4)


IF type EQ 'pf' THEN result = u ELSE BEGIN
 
; Two other quantities we have to derive.

  IF type EQ 'eint' THEN BEGIN
    eint = 0.D0
    print, 'eint'
    print, a0
    print, b0
    ; This is the derivative of u 
    f = t * (b0[0] + b0[1]*t + b0[2]*t^2)
    g = a0[0] + a0[1]*t + a0[2]*t^2 + a0[3]*t^3 + a0[4]*t^4
    fprim = b0[0] + 2*b0[1]*t + 3*b0[2]*t^2
    gprim = a0[1] + 2*a0[2]*t + 3*a0[3]*t^2 + 4*a0[4]*t^3
    eint = fprim/g - gprim*f/g^2
    eint = eint * kkerg * t^2 / u
    result = eint
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
