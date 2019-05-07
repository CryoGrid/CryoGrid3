clear
clc
close all

load geomSetup2.mat 
addpath ..\cryoGridExcessIceInfiltration
oldGRID=load('..\..\runs\190502_14w50y_roundPalsa3m\190502_14w50y_roundPalsa3m_realization4_settings.mat');
oldGRID=oldGRID.GRID;
GRID.general.cT_grid=oldGRID.general.cT_grid;
GRID.soil.cT_domain=oldGRID.soil.cT_domain;
GRID.general.K_delta=oldGRID.general.K_delta;
clear oldGRID

% Stuff necessary to get where we want
PARA.soil.layer_properties=geomSetup2(3).thermalInit(4).layer_properties;
PARA.location.initial_altitude=geomSetup2(3).initial_altitude(4);
GRID=createStratigraphy(PARA,GRID);

% Do the calcuation
deltaSoil=GRID.general.K_delta(GRID.soil.cT_domain);
indi=find(GRID.soil.cT_mineral~=0.05);
indi=indi(1);
indf=find(GRID.soil.cT_mineral==0.5);
indf=indf(1)-1;
deltaSoil=deltaSoil(indi:indf);
heightXice_before = sum(deltaSoil);

[ ~ , finalThick ] = excessGroundIceCheckThickness(GRID, 0.02);
heightXice_after=sum(finalThick(indi:indf));

deltaShrink=heightXice_before-heightXice_after
finalElev= PARA.location.initial_altitude - deltaShrink