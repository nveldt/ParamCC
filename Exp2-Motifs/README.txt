README

All code for running experiments and making plots for motif clustering experiments in paper:

“Parameterized Correlation Clustering \\in Hypergraphs and Bipartite Graphs”
Nate Veldt, Anthony Wirth, and David F. Gleich

—————————————————

OUTSIDE CODE: 
Running experiments requires one to first download/install/compile Metis, Graclus, and hMetis, and add local path to the software in the MATLAB files. 


save_hypergraphs_for_hmetis.jl: 
Run this to save the hypergraphs in the correct format required by hMetis


—————————————————

Florida Bay Experiments

Louvain-based algorithms are run in Florida-Bay-Exps.jl

Recursive Spectral is run in Run_RecSpectral_Flbay.jl

All other algorithms are run using a MATLAB interface in Florida_day_run_algorithms.m

Plots are generated automatically in Florida-Bay-Exps.jl

—————————————————

Email Experiments

Run_RecSpectral_Email.jl for recursive spectral
Email_run_algorithms.m for Metis, Graclus, and hMetis
Email-Louvain-based.jl for Louvain-based algorithms (including for optimizing the HyperLam framework)

Plots generated in Email-all-plots.jl