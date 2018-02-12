function water_albedo = waterAlbedo(sun_elevation,windspeed)
%water albedo parameterization after Wayne and Burt (1954)
%assuming a fixed atmospheric turbitiy T'=3.
%This is a very simple parameterization and might be replaced by a more
%sophisticated model in the future.

u_a=windspeed;

%translate wind speeds [m/s] into wave slopes following suggestions of
%Wayne and Burt (1954)
h2 =      (0.0<=u_a & u_a<0.2) .* 100;
h2 = h2 + (0.2<=u_a & u_a<5.5) .* 30;
h2 = h2 + (5.5<=u_a & u_a<13.9).* 20;
h2 = h2 + (13.9<=u_a).*10;

%tabulation of wave slopes after Wayne and Burt (1954) 
H2 = [100; 30; 20; 10]; 
[~, i_h]=min(abs(H2-h2));

p=sun_elevation;
%tabulation of sun elevation angles after Wayne and Burt (1954)
P=[90 50 30 10]; 
[~, i_p]=min(abs(P-p));

water_albedo_tab = [0.044 0.053 0.089 0.168;
                    0.045 0.050 0.086 0.202;
                    0.045 0.050 0.084 0.219;
                    0.046 0.050 0.080 0.281];
                
water_albedo = water_albedo_tab(i_h,i_p);
                
                
                
