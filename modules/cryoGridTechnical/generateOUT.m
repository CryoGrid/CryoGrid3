function OUT = generateOUT()

    OUT.snow.outSnow_i=[];
    OUT.snow.outSnow_a=[];
    OUT.snow.outSnow_w=[];

    OUT.cryoGrid3=[];
    OUT.water=[];
    OUT.liquidWater=[]; % for distinction between total and liquid water content
    OUT.timestamp=[];
    OUT.TIMESTEP=[];

    OUT.SEB.Lsta=[];
    OUT.SEB.QE=[];
    OUT.SEB.QH=[];
    OUT.SEB.QNET=[];
    OUT.SEB.QG=[];
    OUT.SEB.Tsurf=[];
    OUT.SEB.albedo_stored=[];
    OUT.SEB.Qsurf=[];

    OUT.soil.topPosition=[];
    OUT.soil.lakeFloor=[];
    OUT.soil.soil=cell(0);

    OUT.snow.topPosition=[];
    OUT.snow.botPosition=[];

end