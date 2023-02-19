This file is created to briefly explain the use of R scripts in the repository. 

HT.R 
The R script includes the functions (and its helpers) to perform the rank-based hypothesis testing procedure (Section 2.3 in the shared Overleaf). 

get_data.R 
The R script includes the helper functions that are used to generate testing functions, given copula samples. The output is a data matrix with each column being a functional data.  

tfplot.r 
This is R script origianlly created by Huang Huang. The function provides visualizationsfor the rank-based testing results. For example, Figure 2 in the shared Overleaf is created by tfplot() in DA_nutritions.R.  

Table3.R
It was used to perform simulations for the setup in Table 3 of the shared Overleaf. Note that it includes some helper functions that can be simply called from HT.R or get_data.R.  

Table4_radial.R 
It was used to perform simulations of the setup in Table 4 (radial symmetry) of the shared Overleaf. Note that it includes some helper functions that can be simply called from HT.R or get_data.R 

Table4_joint.R 
It was used to perform simulations of the setpu in Table 4 (joint symmetry) of the shared Overleaf. Note thatit includes some helper functions that can be simply called from HT.R or get_data.R 

DA_nutrition.R
It was used to generate Figure 2 in the shared Overleaf (DA on the nutrition dataset). 

DA_wind_1.R
It was used to perform DA over Wind_Dataset1.rds. To generate plots as in Figure 2, one needs to make some modifications. 

DA_wind_2.R
It was used to perform DA over Wind_Dataset2.rds. To generate plots as in Figure 2, one needs to make some modifications. 


copula_simulation.R
The R script includes the functions to sample from a copula of desired symmetric structures. The input supports empirical copulas. 

simulate_plots.R
The R script was used to generate the scatterplots in pages 13-14 of the shared Overlea.The purpose of the figure is to visually demonstrate that our chosen mixtures of empirical copulas succeed in yielding desired symmetry strutures.   

copula_structure_visualization.R 
The R script has nothing to do with the shared Overleaf. I simply used it to familiarize the fucntional boxplot under the context of copulas. 
