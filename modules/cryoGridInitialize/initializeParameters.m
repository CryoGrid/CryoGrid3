function PARA = initializeParameters(PARA, FORCING)

PARA.surf.albedo=PARA.soil.albedo;
PARA.surf.epsilon=PARA.soil.epsilon;  
PARA.surf.z0=PARA.soil.z0; 
PARA.surf.rs=PARA.soil.rs;  
PARA.snow.albedo=PARA.snow.max_albedo; % sets the initial albedo

if isempty(PARA.technical.starttime)
    PARA.technical.starttime=FORCING.data.t_span(1,1); 
end

if isempty(PARA.technical.endtime)
    PARA.technical.endtime=FORCING.data.t_span(end-1,1); %be on the save side end-1 
end
