;+
; NAME:
;
;   DIATOMICS_TSUJI
;
; PURPOSE:
;
;   This routine evaluates the interpolation polynomial of the chemical   
;   equilibrium constant for a given temperature and diatomic specie. The 
;   coefficients of the polynomial are from Tsuji (1973). For more 
;   information, check the header of diatomics.pro.
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
;   f = DIATOMICS_TSUJI(t, name, type = type, d0 = d0, mab = mab)
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
; DIATOMICS_TSUJI, IDL routine by Nikola Vitas is licensed under a Creative  
; Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
;
; The data in tsuji1973_eqc.sav are not integral part of this routine. If  
; you use the data in a publication, please acknowledge it by citing the proper  
; reference (1973A&A....23..411T).
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
FUNCTION diatomics_tsuji, t, name, type = type, d0 = d0, mab = mab

IF NOT KEYWORD_SET(type) THEN type = 'kp'

kkev    = 8.6173324D-5         ; eV K^-1                  (from NIST)
theta = ALOG10(EXP(1.0))/(kkev * t)

@molecules_path

IF type EQ 'kp' THEN BEGIN

  RESTORE, path + 'catalogue_of_molecules.sav'
  data = tsuji_eqc
  index = WHERE(name EQ data.name)

  kp = t*0.
  coeffs = data[index].coeffs
  nc = N_ELEMENTS(coeffs)

  FOR ii = 0, nc-1 DO kp += coeffs[ii]*theta^ii
  kp = 10.D0^kp

  result = kp

ENDIF

IF type EQ 'pf' OR type EQ 'eint' THEN BEGIN
  PRINT, 'Not available.'
  result = -1
ENDIF

RETURN, result
END
