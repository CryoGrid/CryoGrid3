# CryoGrid3
***Depreciation note:*** *We encourage potentially interested users of CryoGrid to use the newer [community version](https://github.com/CryoGrid/CryoGrid) of the model. It has a much more flexible and modular structure and comprises (almost) all functionalities of the older CryoGrid3 model.*

## Background

**CryoGrid 3** is a land-surface scheme dedicated to modeling of ground temperatures in permafrost environments. Its excess ice module **Xice** is capable of simulating ground subsidence and thermokarst lake formation due to melting of excess ground ice. 

## References

The basic version of **CryoGrid 3** **Xice** is described in the following article:

Westermann, S., Langer, M., Boike, J., Heikenfeld, M., Peter, M., Etzelmüller, B., & Krinner, G. (2016). Simulating the thermal regime and thaw processes of ice-rich permafrost ground with the land-surface model CryoGrid 3. *Geosci. Model Dev.*, 9(2), 523–546. [https://doi.org/10.5194/gmd-9-523-2016](https://doi.org/10.5194/gmd-9-523-2016)

This version has been extended by the lake mode **FLake** and has been used for simulations described in the following article:

Langer, M., Westermann, S., Boike, J., Kirillin, G., Grosse, G., Peng, S., & Krinner, G. (2016). Rapid degradation of permafrost underneath waterbodies in tundra landscapes—Toward a representation of thermokarst in land surface models. *Journal of Geophysical Research: Earth Surface*, *121*(12), 2446–2470. [https://doi.org/10.1002/2016JF003956](https://doi.org/10.1002/2016JF003956)

Version [`v1.0.0`](https://github.com/CryoGrid/CryoGrid3/releases/tag/v1.0.0) has been extended by a hydrology scheme for unfrozen ground conditions, and schemes for the lateral transport of heat, water, and snow between adjacent parts of the simulated environment. It has been set up to study the degradation of ice-wedge polygons as described in the following article:

Nitzbon, J., Langer, M., Westermann, S., Martin, L., Aas, K. S., & Boike, J. (2019). Pathways of ice-wedge degradation in polygonal tundra under different hydrological conditions. *The Cryosphere*, 13(4), 1089–1123. [https://doi.org/10.5194/tc-13-1089-2019](https://doi.org/10.5194/tc-13-1089-2019)

Very similar model versions have been used to study the degradation of peat plateaus and palsas as described in the following articles:

Martin, L. C. P., Nitzbon, J., Aas, K. S., Etzelmüller, B., Kristiansen, H., & Westermann, S. (2019). Stability Conditions of Peat Plateaus and Palsas in Northern Norway. *Journal of Geophysical Research: Earth Surface*, *124*, 705–719. [https://doi.org/10.1029/2018JF004945](https://doi.org/10.1029/2018JF004945)

Martin, L. C. P., Nitzbon, J., Scheer, J., Aas, K. S., Eiken, T., Langer, M., Filhol, S., Etzelmüller, B., & Westermann, S. (2020). Thermal erosion patterns of permafrost peat plateaus in northern Norway. *The Cryosphere Discussions*, 1–33. [https://doi.org/10.5194/tc-2020-338](https://doi.org/10.5194/tc-2020-338)

Version [`v1.1.0`](https://github.com/CryoGrid/CryoGrid3/releases/tag/v1.1.0) has been extended by a scheme for the lateral transport of sediment within the simulated environment. It has been used for the simulations related to a research article on the response of ice-rich permafrost to a warming climate:

Nitzbon, J., Westermann, S., Langer, M., Martin, L. C. P., Strauss, J., Laboor, S., & Boike, J. (2020). Fast response of cold ice-rich permafrost in northeast Siberia to a warming climate. *Nature Communications*, 11, 2201. [https://doi.org/10.1038/s41467-020-15725-8](https://doi.org/10.1038/s41467-020-15725-8)

Version [`v1.2.0`](https://github.com/CryoGrid/CryoGrid3/releases/tag/v1.2.0) has been extended by a multi-scale tiling scheme. It has been used for simulations presented in a research article published in *The Cryosphere*: 

Nitzbon, J., Langer, M., Martin, L. C. P., Westermann, S., Schneider von Deimling, T., & Boike, J. (2020). Effects of multi-scale heterogeneity on the simulated evolution of ice-rich permafrost lowlands under a warming climate. *The Cryosphere*, 15(3), 1399–1422. [https://doi.org/10.5194/tc-15-1399-2021](https://doi.org/10.5194/tc-15-1399-2021)