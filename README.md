# rough_surface_contact

This repository contains the source code to reproduce the contact evolution results in the paper "Elastic contact of rough surfaces with fractal and Hurst effects." 

The file "rough_surface_gen.R" generates rough surfaces with either Cauchy or Dagum autocorrelation functions. 

The file "dcfft_g.m" solves the elastic contact problem for a given rough surface and outputs the contact evolution data. The dataset "data_emwes2048tm.mat" file corresponds to the influence coefficient data, which is necessary to run "dcfft_g.m" file. 

The files "kappap_g.m" and "slope_g.m" plots $\kappa$, $\kappa^{\prime}$, and $\kappa^{\prime\prime}$ as a function of contact area evolution. 

The file "Supplementary_Material.pdf" contains the contact area results corresponding to both the Cauchy and Dagum random surfaces
without contact area correction near the perimeter of the contact zones.  
