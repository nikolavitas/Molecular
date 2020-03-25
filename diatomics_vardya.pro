;+
; NAME:
;
;   DIATOMICS_VARDYA
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial of the chemical   
;   equilibrium constant and the internal energy of the rotational-vibrational 
;   bound states of H2 and H2+ using the polynomial fits of Vardya (1961, 1965). 
;   For more information, check the header of diatomics.pro.
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
;   f = DIATOMICS_VARDYA(t, name, type = type, d0 = d0, mab = mab)
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
; DIATOMICS_VARDYA, IDL routine by Nikola Vitas is licensed under a Creative 
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in used in this routine comes from two papers of Vardya 
; (1961ApJ...133..107V and 1965MNRAS.129..205V).  If you use the data in a 
; publication, please acknowledge it by citing the corresponding references.

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

FUNCTION diatomics_vardya, t, name, type = type, d0 = d0, mab = mab

kkerg   = 1.3806488D-16        ; erg K^-1                 (from NIST)
kkev    = 8.6173324D-5         ; eV K^-1                  (from NIST)

IF NOT KEYWORD_SET(type) THEN type = 'kp'

theta = ALOG10(EXP(1.0))/(kkev * t)

IF type EQ 'kp' THEN BEGIN

  IF name EQ 'H2' THEN BEGIN
    kp = 12.533505D0 - theta*4.9251644D0 + theta^2*5.6191273D-2 - theta^3*3.2687661D-3
    kp = 10.^kp
  ENDIF

  IF name EQ 'H2+' THEN BEGIN
    kp = 11.206998D0 - theta*2.7942767D0 + theta^2*7.9196803D-2 - theta^3*2.4790744D-2
    kp = 10.^kp
  ENDIF

  result = kp

ENDIF

IF type EQ 'eint' THEN BEGIN

  IF name EQ 'H2' THEN BEGIN
    eint = (2.6757D0 - 1.4772D0*theta + 0.60602D0*theta^2 - 0.12427D0*theta^3 + 0.0097503D0*theta^4) * (kkerg * t)
  ENDIF

  IF name EQ 'H2+' THEN BEGIN
   eint = (2.9216D0 - 2.0036D0*theta + 1.7231D0*theta^2 - 0.82685D0*theta^3 + 0.15253D0*theta^4) * (kkerg * t)
  ENDIF

  result = eint

ENDIF

IF type EQ 'pf' THEN BEGIN
  PRINT, 'Not available.'
  result = -1
ENDIF

RETURN, result
END
