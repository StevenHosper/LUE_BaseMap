# LUE Hydrologic BaseMAP
## Background information
The LUE Framework is being developed by the [computational geography team of the University Utrecht](https://github.com/computationalgeography).
This model is a first attempt at creating a basic hydrologic model for the Netherlands using the LUE Framework.
In this project Nelen & Schuurmans and the University Utrecht work together in order to create a state-of-the-art hydrologic model.

## Research goals
Currently many of the NHI uses a lumped model for many of the watersheds in the Netherlands. 
In order to improve this and make it easier to adapt when a situation changes it seems preferable to create a fully distributed model. 
With recent technological improvements, like the LUE modelling framework, the possibility of a fully distributed model for the Netherlands has increased.
The long-term goal is to accurately represent the Netherlands at high-resolution with simulations based on atmospheric inputs.
The initial study will investigate the capabilities of LUE framework to build a hydrologic model that can accurately represent the Hupselse Beek region in the East of the Netherlands. 
This region is only 25 km2, and ideally the entire Netherlands would be modelled eventually, but hardware is limited. 
During development 16 GB of RAM and an 11th Gen Intel(R) Core(TM) i5-1135G7 @ 2.40GHz processor are used.
Developing a smaller array will improve the speed of model development and allow for quick adaptations.

## Model concept
The main goal was to create a surface water routing model. 
During initial creating of the model it became aparent that groundwater would have to be included.
Model configuration is currently determined within an ini file, of which a template can be found within the config folder.
Data for atmospheric fluxes is extracted from csv files in mm/hour and converted to a m/s rate.

Initial conditions are set if files are suplied for the initial conditions of one of the following: initial groundwater storage, initial interception storage, initial surface water height. The default value is zero.

Currently the model includes the following fluxes:
*Vertical*
- Precipitation
- Evaporation
- Infiltration
- Interception

*Lateral*
- Discharge
- Groundwater flow
- Seepage

Currently new flux values are set every 60 seconds, at this same time interval outputs are reported.
