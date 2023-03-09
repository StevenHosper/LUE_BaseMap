# LUE Hydrologic BaseMAP
## Background information
The LUE Framework is being developed by the [computational geography team of the University Utrecht](https://github.com/computationalgeography).
This model is a first attempt at creating a basic hydrologic model for the Netherlands using the LUE Framework.
In this project Nelen & Schuurmans and the University Utrecht work together in order to create a state-of-the-art hydrologic model.

## Research goals
Currently many of the NHI uses a lumped model for many of the watersheds in the Netherlands. 
In order to improve this and make it easier to adapt when a situation changes it seems preferable to create a fully distributed model. 
With recent technological improvements, like LUE, the possibility of a fully distributed model for the Netherlands increased.
The goal is to accurately represent large parts of the Netherlands at high-resolution with simulations based on atmospheric inputs.
Large, in ideal situations, would mean the entire Netherlands, however, the hardware during this project is limited.
Therefore initially we will work with a array of 1000x1000 and 5000x5000 size.
This will improve the speed of model development and allow for quick adaptations.

## Model concept
Currently the model is mainly planned to function for surface water routing.
This will supply data required by the ModFlow-MetaSwap groundwater model that is currently being used for modelling purposes in general.
Ideally during this project the coupling will have already occured to make the results of the model more accurate.
However, as of yet it is not clear how fast development will go, and if LUE is even capable of the attaining accurate results for high-resolution large-scale modelling.
For now the model is planned to incorporate atleast the following surface processes:
- Precipitation
- Evaporation
- Infiltration
- Abstraction
- Injection

If groundwater processes are added, percolation and seepage will be the first to be added.
A good connection to the groundwater is preferred as many areas in the Netherlands are heavily influenced by groundwaterflows. 
However, initially to deminish model complexities, groundwater flow will be ignored.
