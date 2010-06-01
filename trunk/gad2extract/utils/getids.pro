filevec=['4096_R300__0_dump_370', $
       '4096_R300__1_dump_370',$ 
       '4096_R300__2_dump_370']
openw, 1, "id.idx"
printf, 1, '# number of halos, np ID1 ID2... etc'
printf, 1, 3
window, xsize=300, ysize=900
!P.multi=[0, 1, 3]
for i=0, 2 do begin
file=filevec[i]
readnew, file, posin, 'POS ', parttype=4, /quiet
readnew, file, idin, 'ID  ', parttype=4, /quiet
readnew, file+'_rho_4', rhoin, 'H4  '
WCOM=[0.0,0.0,0.0]
R=50.0
GetIDInReg, posin, WCOM,R, idinreg, SQ=1
get_w_com, posin(*,idinreg), rhoin(idinreg), WCom

plot, posin(0,*), posin(1,*), psym=3, /xs, /ys, /iso
plots, WCOM(0), WCOM(1), psym=1, color=fsc_color('blue'), symsize=10

R=4.0
GetIDInReg, posin, WCOM,R, idinreg, SQ=1
np=n_elements(idinreg)
print,'Halo ',i, '  Np=', np
printf, 1, np
printf, 1, [idin(idinreg)]
endfor 


close, 1
end
