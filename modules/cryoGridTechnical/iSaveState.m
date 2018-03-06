function iSaveState( fname, T, wc, t, SEB, PARA, GRID)

FINAL.T=T;
FINAL.wc=wc;
FINAL.t=t;
FINAL.SEB=SEB;
FINAL.PARA=PARA;
FINAL.GRID=GRID;

save( fname, 'FINAL' );
end