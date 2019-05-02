% Work on input of thermal initialisation

% Create variable
thermalInit(7).Tinitial=[];

% T initial
load('190426_7w100y_roundPalsa15m_realization1_settings.mat')
thermalInit(1).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization2_settings.mat')
thermalInit(2).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization3_settings.mat')
thermalInit(3).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization4_settings.mat')
thermalInit(4).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization5_settings.mat')
thermalInit(5).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization6_settings.mat')
thermalInit(6).Tinitial=PARA.Tinitial;
load('190426_7w100y_roundPalsa15m_realization7_settings.mat')
thermalInit(7).Tinitial=PARA.Tinitial;

% Active Layer + reorder
thermalInit(7).ActiveLayer=[];
thermalInit(7).layer_properties=[];
thermalInit = orderfields(thermalInit, [2,1,3]);
thermalInit = orderfields(thermalInit, [1,3,2]);
ActiveLayer=[NaN 0.9 0.9 0.85 0.8 0.75 0.7];
ActiveLayer=num2cell(ActiveLayer);
[thermalInit.ActiveLayer]=ActiveLayer{:};

% Layer properties
load('190426_7w100y_roundPalsa15m_realization1_settings.mat')
thermalInit(1).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization2_settings.mat')
thermalInit(2).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization3_settings.mat')
thermalInit(3).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization4_settings.mat')
thermalInit(4).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization5_settings.mat')
thermalInit(5).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization6_settings.mat')
thermalInit(6).layer_properties=PARA.soil.layer_properties;
load('190426_7w100y_roundPalsa15m_realization7_settings.mat')
thermalInit(7).layer_properties=PARA.soil.layer_properties;