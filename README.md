## Biotic velocity (i.e., migration rate)
### Data: 
#### Biotic velocities (in meters/year) for 15 North American tree taxa (plus an 'other' group that combines remaining tree taxa that weren't modeled individually). 
Biotic velocities were calculated using the bioticVelocity function in the enmSdm R package (https://github.com/adamlilith/enmSdm), which requires rasterstacks of habitat suitability over time. 
Suitability surfaces were modeled using estimation and prediction steps (functions used: pg_stlm_latent_overdispersed and predict_pg_stlm_latent_overdispersed, respectively) in the pgR R package (https://github.com/jtipton25/pgR).

### R: 
R code used to plot biotic velocities.

### Figures: 
Plots of biotic velocities. 
