function SEB = initializeSEB()

size_L_star_smoothing=1;

SEB.Qnet=0; 
SEB.Qh=0;
SEB.Qe=0; 
SEB.Qg=0;
SEB.Sout=0; 
SEB.Lout=0;
SEB.newSnow=0;
<<<<<<< HEAD
%SEB.sublim=0;
%SEB.meltwater=0;
=======
>>>>>>> origin/xice_mpi_polygon_TC
SEB.L_star=-100000+zeros(1,size_L_star_smoothing);
SEB.u_star=10;


SEB.Qsurf = 0;  % for EB checks
