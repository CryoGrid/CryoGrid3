clear
clc
close all

if ispc
    slash='\';
else
    slash='/';
end

addpath(['..' slash 'cryoGridSoil'])

geomSetup(2).descr=[];

%% Set up of the 26/04/2019

% General description
geomSetup(1).descr='round palsa';
geomSetup(1).numlabs=7;

% Topological relashioships
rad=2; % radius of central palsa
ring=0.4; % width of the surounding rings
geomSetup(1).weight = [ pi*(8^2-(rad+5*ring)^2)   pi*((rad+5*ring)^2-(rad+4*ring)^2)   pi*((rad+4*ring)^2-(rad+3*ring)^2)   pi*((rad+3*ring)^2-(rad+2*ring)^2)   pi*((rad+2*ring)^2-(rad+1*ring)^2)   pi*((rad+1*ring)^2-rad^2)   pi*rad^2;]; % Has to be an integer
area_tot = sum(geomSetup(1).weight);
geomSetup(1).area = geomSetup(1).weight ./ sum(geomSetup(1).weight) .* area_tot ; % in m^2
geomSetup(1).weight=round(geomSetup(1).weight);
dist=[2.20 0.20 0.20 0.20 0.20 2.20];
A=zeros(geomSetup(1).numlabs,geomSetup(1).numlabs);
idx = sub2ind(size(A),[1:geomSetup(1).numlabs-1 2:geomSetup(1).numlabs],[2:geomSetup(1).numlabs 1:geomSetup(1).numlabs-1]);
A(idx) = [dist dist];
geomSetup(1).distanceBetweenPoints= A; %   %in m. Put 0 for all non-connected ensemble members
A = double( geomSetup(1).distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)
geomSetup(1).A=A;

% Topographical relationships
geomSetup(1).initial_altitude=linspace(300,301.5,7);

% Heat exchange
B=A;
perimeters=[2*pi*(rad+5*ring) 2*pi*(rad+4*ring) 2*pi*(rad+3*ring) 2*pi*(rad+2*ring) 2*pi*(rad+1*ring) 2*pi*rad];
B(idx)= [perimeters perimeters];
geomSetup(1).thermal_contact_length = B;
B=A;
thdist=[0.20 0.20 0.20 0.20 0.20 0.20];
B(idx)= [thdist thdist];
geomSetup(1).thermalDistance = B;

% Water exchange
geomSetup(1).boundaryCondition={'DarcyReservoir','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC'};
geomSetup(1).Darcy_elevation=[300 nan nan nan nan nan nan];
geomSetup(1).Darcy_fluxFactor=[1000*5*1e-5/1.20 nan nan nan nan nan nan];

% Snow exchange
geomSetup(1).immobile_snow_height = [ 0.05 0.05 0.05 0.05 0.05 0.05 0.05];

% Thermal Init
thermalInit(7).ActiveLayer=NaN;
ActiveLayer={NaN;0.900000000000000;0.900000000000000;0.850000000000000;0.800000000000000;0.750000000000000;0.700000000000000};
[thermalInit.ActiveLayer]=ActiveLayer{:};
thermalInit(1).layer_properties=[0 0.800000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.500000000000000 0.800000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(2).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.900000000000000 0.821276595744681 0.0446808510638298 0.134042553191489 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(3).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.900000000000000 0.838461538461539 0.0403846153846154 0.121153846153846 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(4).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.850000000000000 0.851724137931035 0.0370689655172414 0.111206896551724 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(5).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.800000000000000 0.862500000000000 0.0343750000000000 0.103125000000000 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(6).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.750000000000000 0.871428571428571 0.0321428571428572 0.0964285714285714 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(7).layer_properties=[0 0.550000000000000 0.0500000000000000 0.150000000000000 1 0.800000000000000;0.700000000000000 0.878947368421053 0.0302631578947368 0.0907894736842105 1 0.800000000000000;3 0.500000000000000 0.500000000000000 0 2 0.500000000000000;10 0.0300000000000000 0.970000000000000 0 1 0.0300000000000000];
thermalInit(1).Tinitial=[-5 10;0 5;0.100000000000000 2;0.500000000000000 1.50000000000000;1 2;2 2;10 2;30 2;500 4;1000 10];
thermalInit(2).Tinitial=[-5 10;0 5;0.100000000000000 2;0.900000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
thermalInit(3).Tinitial=[-5 10;0 5;0.100000000000000 2;0.900000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
thermalInit(4).Tinitial=[-5 10;0 5;0.100000000000000 2;0.850000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
thermalInit(5).Tinitial=[-5 10;0 5;0.100000000000000 2;0.800000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
thermalInit(6).Tinitial=[-5 10;0 5;0.100000000000000 2;0.750000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
thermalInit(7).Tinitial=[-5 10;0 5;0.100000000000000 2;0.700000000000000 0;1 -1;2 -1;10 1;30 2;500 4;1000 10];
geomSetup(1).thermalInit=thermalInit;
% clear related var
clear A area_tot B dist idx perimeters rad ring thdist

%% Set up of the 02/05/2019

% General description
geomSetup(2).descr='round palsa';
geomSetup(2).numlabs=14;

% Topological relashioships
radius=[0.6:0.6:6 7:1:9 12];
area=pi*radius.^2;
area= [area(1) area(2:end)-area(1:end-1)];
weight=round(area);
geomSetup(2).area = fliplr(area); % fliplr because so far we have been putting the mire on the left and the palsa on the right
geomSetup(2).weight= fliplr(weight);
dist=[2.10 0.6*ones(1,12)];
A=zeros(geomSetup(2).numlabs,geomSetup(2).numlabs);
idx = sub2ind(size(A),[1:geomSetup(2).numlabs-1 2:geomSetup(2).numlabs],[2:geomSetup(2).numlabs 1:geomSetup(2).numlabs-1]);
A(idx) = [dist dist];
geomSetup(2).distanceBetweenPoints= A; %   %in m. Put 0 for all non-connected ensemble members
A = double( geomSetup(2).distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)
geomSetup(2).A=A;

% Topographical relationships
geomSetup(2).initial_altitude=[300 300:0.25:302 302*ones(1,4)];

% Heat exchange
B=A;
perimeters=fliplr(2*pi*radius(1:end-1));
B(idx)= [perimeters perimeters];
geomSetup(2).thermal_contact_length = B;
B=A;
thdist=0.60*ones(1,13);
B(idx)= [thdist thdist];
geomSetup(2).thermalDistance = B;

% Water exchange
geomSetup(2).boundaryCondition={'DarcyReservoir','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC', 'NoBC'};
geomSetup(2).Darcy_elevation=[300 nan(1,13)];
geomSetup(2).Darcy_fluxFactor=[1000*5*1e-5/1.20 nan(1,13)];

% Snow exchange
geomSetup(2).immobile_snow_height = 0.05.*ones(1,14);

% Thermal init
% Create variable
thermalInit(14).ActiveLayer=[];
thermalInit(14).layer_properties=[];
thermalInit(14).Tinitial=[];

% Active layers
ActiveLayer=[NaN NaN linspace(0.9,0.7,12)]; % Input active layers
ispf=~isnan(ActiveLayer); % Define who is permafrost

% Stratigraphy
Strati_mire =[    0.0     0.80    0.05    0.15    1   0.80    ;...
                  0.5     0.80    0.05    0.15    1   0.80    ;...
                  3.0     0.50    0.50    0.00    2   0.50    ;...
                 10.0     0.03    0.97    0.00    1   0.03   ];
Strati_palsa_initial=Strati_mire; % Start from the mire strati
Strati_palsa_initial(1,2)=0.55; % Set the upper layer at FC
palsaHeight=geomSetup(2).initial_altitude - min(geomSetup(2).initial_altitude);

for i=1:14;
    if ispf(i)==0;
        thermalInit(i).layer_properties=Strati_mire;
    else
        thermalInit(i).layer_properties=stratiXice(Strati_palsa_initial, palsaHeight(i), ActiveLayer(i)); % Modify the ice, mineral and organic content
    end
end

% Initial T profiles
T_z    =[-5 0 0.1 0.5 1 2 10 30 500 1000]';
T_mire =[10 5 2 1.5 2 2 2 2 4 10]';
T_palsa=[10 5 2 0 -1 -1 1 2 4 10]';
ActiveLayer(isnan(ActiveLayer))=0.5;
assert(min(ActiveLayer)>0.1 && max(ActiveLayer)<1,'Check initial active layers and initial T profile')
for i=1:14;
    T_z(4)=ActiveLayer(i);
    if ispf(i)==0;
        thermalInit(i).Tinitial=[T_z T_mire];
    else
        thermalInit(i).Tinitial=[T_z T_palsa];
    end
end

ActiveLayer=num2cell(ActiveLayer);
[thermalInit.ActiveLayer]=ActiveLayer{:};

geomSetup(2).thermalInit=thermalInit;

% clear related var
clear area radius weight A B dist idx perimeters thdist ActiveLayer i ispf palsaHeight Strati_mire Strati_palsa_initial T_mire T_palsa T_z thermalInit

%% Set up of the 02/05/2019 - Higher Palsa

% General description
geomSetup(3).descr='round palsa';
geomSetup(3).numlabs=14;

% Topological relashioships
radius=[0.6:0.6:6 7:1:9 12];
area=pi*radius.^2;
area= [area(1) area(2:end)-area(1:end-1)];
weight=round(area);
geomSetup(3).area = fliplr(area); % fliplr because so far we have been putting the mire on the left and the palsa on the right
geomSetup(3).weight= fliplr(weight);
dist=[2.10 0.6*ones(1,12)];
A=zeros(geomSetup(3).numlabs,geomSetup(3).numlabs);
idx = sub2ind(size(A),[1:geomSetup(3).numlabs-1 2:geomSetup(3).numlabs],[2:geomSetup(3).numlabs 1:geomSetup(3).numlabs-1]);
A(idx) = [dist dist];
geomSetup(3).distanceBetweenPoints= A; %   %in m. Put 0 for all non-connected ensemble members
A = double( geomSetup(3).distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)
geomSetup(3).A=A;

% Topographical relationships
geomSetup(3).initial_altitude=[300 linspace(300,303,9) 303*ones(1,4)];

% Heat exchange
B=A;
perimeters=fliplr(2*pi*radius(1:end-1));
B(idx)= [perimeters perimeters];
geomSetup(3).thermal_contact_length = B;
B=A;
thdist=0.60*ones(1,13);
B(idx)= [thdist thdist];
geomSetup(3).thermalDistance = B;

% Water exchange
geomSetup(3).boundaryCondition={'DarcyReservoir','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC','NoBC', 'NoBC'};
geomSetup(3).Darcy_elevation=[300 nan(1,13)];
geomSetup(3).Darcy_fluxFactor=[1000*5*1e-5/1.20 nan(1,13)];

% Snow exchange
geomSetup(3).immobile_snow_height = 0.05.*ones(1,14);

% Thermal init
% Create variable
thermalInit(14).ActiveLayer=[];
thermalInit(14).layer_properties=[];
thermalInit(14).Tinitial=[];

% Active layers
ActiveLayer=[NaN NaN linspace(0.9,0.7,12)]; % Input active layers
ispf=~isnan(ActiveLayer); % Define who is permafrost

% Stratigraphy
Strati_mire =[    0.0     0.80    0.05    0.15    1   0.80    ;...
                  0.5     0.80    0.05    0.15    1   0.80    ;...
                  3.0     0.50    0.50    0.00    2   0.50    ;...
                 10.0     0.03    0.97    0.00    1   0.03   ];
Strati_palsa_initial=Strati_mire; % Start from the mire strati
Strati_palsa_initial(1,2)=0.55; % Set the upper layer at FC
palsaHeight=geomSetup(3).initial_altitude - min(geomSetup(3).initial_altitude);

for i=1:14;
    if ispf(i)==0;
        thermalInit(i).layer_properties=Strati_mire;
    else
        thermalInit(i).layer_properties=stratiXice(Strati_palsa_initial, palsaHeight(i), ActiveLayer(i)); % Modify the ice, mineral and organic content
    end
end

% Initial T profiles
T_z    =[-5 0 0.1 0.5 1 2 10 30 500 1000]';
T_mire =[10 5 2 1.5 2 2 2 2 4 10]';
T_palsa=[10 5 2 0 -1 -1 1 2 4 10]';
ActiveLayer(isnan(ActiveLayer))=0.5;
assert(min(ActiveLayer)>0.1 && max(ActiveLayer)<1,'Check initial active layers and initial T profile')
for i=1:14;
    T_z(4)=ActiveLayer(i);
    if ispf(i)==0;
        thermalInit(i).Tinitial=[T_z T_mire];
    else
        thermalInit(i).Tinitial=[T_z T_palsa];
    end
end

ActiveLayer=num2cell(ActiveLayer);
[thermalInit.ActiveLayer]=ActiveLayer{:};

geomSetup(3).thermalInit=thermalInit;

% clear related var
clear area radius weight A B dist idx perimeters thdist ActiveLayer i ispf palsaHeight Strati_mire Strati_palsa_initial T_mire T_palsa T_z thermalInit

%% Save
clear slash
geomSetup2=geomSetup;
save('geomSetup2.mat','geomSetup2')