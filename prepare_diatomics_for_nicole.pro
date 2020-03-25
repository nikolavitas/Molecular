PRO prepare_diatomics_for_nicole

@diatomics_path
filename = path + 'catalogue_of_diatomics.sav'
RESTORE, filename

nmol = N_ELEMENTS(catalogue.name)

IF !D.NAME EQ 'WIN' THEN newline = STRING([13B, 10B]) ELSE newline = STRING(10B)
  

OPENW, 1, 'diatomics_for_nicole.f90.txt'

; Molecular name, string
line = '  Data Molecule_name/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += "'" + catalogue[im].name +"'"
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 8 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Molecular identifier
line = '  Data Molecule_id/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += STRCOMPRESS(STRING(im), /rem)
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 14 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Constituent 1, string
line = '  Data Molecule_char1/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += "'" + catalogue[im].constituents[0] +"'"
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 8 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Constituent 2, string
line = '  Data Molecule_char2/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += "'" + catalogue[im].constituents[1] +"'"
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 8 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Molecular code: z1z1.ion
line = '  Data Molecule_code/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += "'" + catalogue[im].code + "'"
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 6 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

el = load_list_of_elements(indgen(92)+1)

; Constituent 1, atomic number
line = '  Data Molecule_z1/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  z1 = WHERE(el EQ catalogue[im].constituents[0]) + 1
  line += STRCOMPRESS(STRING(z1), /rem)
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 12 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Constituent 2, atomic number
line = '  Data Molecule_z2/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  z2 = WHERE(el EQ catalogue[im].constituents[1]) + 1
  line += STRCOMPRESS(STRING(z2), /rem)
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 14 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Charge (0 or 1)
line = '  Data Molecule_charge/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += STRCOMPRESS(STRING(catalogue[im].charge), /rem)
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 14 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Number of unique elements (1 or 2)
line = '  Data Molecule_nuc/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += STRCOMPRESS(STRING(catalogue[im].nuc), /rem)
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 14 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Dissociation constant
line = '  Data Molecule_d0/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += STRING(catalogue[im].d0, format = '(F13.8)')
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 4 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

; Reduced mass
line = '  Data Molecule_mab/'
linecount = 0
FOR im = 0, nmol-1 DO BEGIN
  line += STRING(catalogue[im].mab, format = '(F13.8)')
  IF im LT nmol-1 THEN line += ", "
  IF linecount EQ 4 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

RESTORE, path + 'sauvalandtatum1984_pf.sav'

; Number in Sauval's list is the same as i

; Reduced mass
line = '  Data Molecule_ST1984_pf/ &' + newline
linecount = 0
FOR im = 0, N_ELEMENTS(data.name)-1 DO BEGIN
  FOR ic = 0, 4 DO BEGIN
    line += STRING(data[im].coeffs[ic], format = '(F13.8)')
    IF ic LT 4 THEN line += ", "
  ENDFOR
  IF im LT N_ELEMENTS(data.name)-1 THEN line += ", "
  IF linecount EQ 0 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line

RESTORE, path + 'sauvalandtatum1984_eqc.sav'

; Reduced mass
line = '  Data Molecule_ST1984_eqc/ &' + newline
linecount = 0
FOR im = 0, N_ELEMENTS(data.name)-1 DO BEGIN
  FOR ic = 0, 5 DO BEGIN
    line += STRING(data[im].coeffs[ic], format = '(F13.8)')
    IF ic LT 5 THEN line += ", "
  ENDFOR
  IF im LT N_ELEMENTS(data.name)-1 THEN line += ", "
  IF linecount EQ 0 THEN BEGIN
    line += ' &' + newline + '     '
    linecount = 0
  ENDIF ELSE linecount +=1
ENDFOR
line += '/'
PRINTF, 1, line


CLOSE, 1
stop


END