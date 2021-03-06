# JAMES2022-Optimal

[![DOI](https://zenodo.org/badge/497347062.svg)](https://zenodo.org/badge/latestdoi/497347062)

Matlab scripts and functions relevant to the paper:
  
Snyder, C. and G. J. Hakim, 2022: An optimal linear transformation for data assimilation, *Journal of Advances in Modeling Earth Systems*, accepted. 

tmp_2scalesRevisted.m pertains to the illustrative example of section 2.  It generates realizations of the state, observations, and state ensemble, then performs updates of various forms.

tmp_2scalesPlot.m plots results from tmp_2scaleRevisited.m.

plot_realizations.m plot realizations of specified random processes, similar to Fig. 2.

test_updateWithOptlTransf.m pertains to the one-dimensional example of section 4. It generates realizations of the state, observations, and state ensemble, then performs updates of various forms.

wrapper_bestFixedLocalization.m sets parameters and calls test_updateWithOptlTransf.m, allowing accumulation of results as ensemble size and localization length vary.

fig3_optimalPaper.m reads mmse_6April_Ne1024.mat, which contains results from wrapper_updateAsFnAnything.m, and plots Fig. 3.
