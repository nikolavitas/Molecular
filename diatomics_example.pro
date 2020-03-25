pro diatomics_example

!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')
WINDOW, xsi = 500, ysi = 400
t = 2000. + DINDGEN(6)*2000.

WINDOW, xsi = 500, ysi = 400
PLOT,  t,diatomics(t, 'H2', data = 's'), lines = 0, thick = 2, title = "Partition function for neutral H2", $
       pos = [0.12, 0.14, 0.94, 0.92], xtit = 'T [K]', ytit = 'U', charsi = 1.4
oplot, t, diatomics(t, 'H2', data = 'i'), col = getcolor('Steel Blue') , thick = 2
oplot, t, diatomics(t, 'H2', data = 'b'), col = getcolor('Orange Red') , thick = 2
OPLOT, [2400, 3400], [200, 200], thick = 2
OPLOT, [2400, 3400], [185, 185], color = getcolor('Steel Blue'), thick = 2
OPLOT, [2400, 3400], [170, 170], color = getcolor('Orange Red'), thick = 2
XYOUTS, 3450, 198, 'Sauval and Tatum', chars = 1.4
XYOUTS, 3450, 183, 'Irwin', chars = 1.4
XYOUTS, 3450, 168, 'Bohn', chars = 1.4
WRITE_PNG, 'example_diatomic_pf_h2.png', tvrd(/true)


!p.background = GETCOLOR('white')
!p.color = GETCOLOR('black')
WINDOW, xsi = 500, ysi = 400
WINDOW, xsi = 500, ysi = 400
PLOT,  t, alog10(diatomics(t, 'H2', data = 's', type = 'kp')), lines = 0, thick = 2, title = "Chemical equilibrium const. for neutral H2", $
       pos = [0.12, 0.14, 0.94, 0.92], xtit = 'T [K]', ytit = 'log Kp', charsi = 1.4
oplot, t, alog10(diatomics(t, 'H2', data = 'i', type = 'kp')), col = getcolor('Steel Blue') , thick = 2
oplot, t, alog10(diatomics(t, 'H2', data = 'b', type = 'kp')), col = getcolor('Orange Red') , thick = 2
oplot, t, alog10(diatomics(t, 'H2', data = 't', type = 'kp')), col = getcolor('Forest Green') , thick = 2
oplot, t, alog10(diatomics(t, 'H2', data = 'v', type = 'kp')), col = getcolor('Chocolate') , thick = 2
;PLOT,  t, diatomics(t, 'H2', data = 's', type = 'kp'), lines = 0, thick = 2, title = "Chemical equilibrium const. for neutral H2", $
;       pos = [0.12, 0.14, 0.94, 0.92], xtit = 'T [K]', ytit = 'log Kp', charsi = 1.4
;oplot, t, diatomics(t, 'H2', data = 'i'), col = getcolor('Steel Blue') , thick = 2
;oplot, t, diatomics(t, 'H2', data = 'b'), col = getcolor('Orange Red') , thick = 2
;oplot, t, diatomics(t, 'H2', data = 't'), col = getcolor('Forest Green') , thick = 2
;oplot, t, diatomics(t, 'H2', data = 'v'), col = getcolor('Chocolate') , thick = 2
OPLOT, [6400, 7400], [6, 6], thick = 2
OPLOT, [6400, 7400], [5.2, 5.2], color = getcolor('Steel Blue'), thick = 2
OPLOT, [6400, 7400], [4.4, 4.4], color = getcolor('Orange Red'), thick = 2
OPLOT, [6400, 7400], [3.6, 3.6], color = getcolor('Forest Green'), thick = 2
OPLOT, [6400, 7400], [2.8, 2.8], color = getcolor('Chocolate'), thick = 2
XYOUTS, 7450, 5.8, 'Sauval and Tatum', chars = 1.4
XYOUTS, 7450, 5.0, 'Irwin', chars = 1.4
XYOUTS, 7450, 4.2, 'Bohn', chars = 1.4
XYOUTS, 7450, 3.4, 'Tsuji', chars = 1.4
XYOUTS, 7450, 2.6, 'Vardya', chars = 1.4
WRITE_PNG, 'example_diatomic_kp_h2.png', tvrd(/true)



OPLOT, [2400, 3400], [200, 200], thick = 2
OPLOT, [2400, 3400], [185, 185], color = getcolor('Steel Blue'), thick = 2
OPLOT, [2400, 3400], [170, 170], color = getcolor('Orange Red'), thick = 2
XYOUTS, 3450, 198, 'Sauval and Tatum', chars = 1.4
XYOUTS, 3450, 183, 'Irwin', chars = 1.4
XYOUTS, 3450, 168, 'Bohn', chars = 1.4
WRITE_PNG, 'example_diatomic_pf_h2.png', tvrd(/true)


end
