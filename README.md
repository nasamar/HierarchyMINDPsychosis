# Hierarchy and MIND in schizophrenia spectrum disorders (SSD)

This repository contains code and data created in support of the project García-San-Martín, N.; Bethlehem, R. AI; Sebenius, I. et al. *"Long-term morphometric similarity gradients relate to cortical hierarchy and psychiatric symptoms in schizophrenia"* medRxiv (2026). All code was written in Matlab, R, and Python. Folders, files, and first steps are described below.


## **Data**

The `Data` folder contains all the data required for running the analyses. Here are the files that need to be downloaded and stored in a specific location. The remaining files will be automatically generated:

- `volumes` folder contains the regional volume and velocity peaks in `Table_2_2.csv` obtained from https://doi.org/10.1038/s41586-022-04554-y. 
  
-	The code used to compute MIND networks is available at https://github.com/isebenius/MIND and corresponds to [MIND_long_01_MIND.py](Code/MIND_long_01_MIND.py).

- The `sensorimotor-association_axis_ranking_DK.csv` file was derived from https://doi.org/10.1016/j.neuron.2021.06.016.

-	The `all_microsc_DesikanKilliany68.csv` file is available at https://github.com/netneurolab/netneurotools.



## **Code**

The `Code` folder contains all the code required for running the analyses and generate data and figures. All scripts are designed to be sequentially executed. Don't forget to change the location variable regularly. 

-	[MIND_long_01_MIND.py](Code/MIND_01_MIND.py) – computes MIND networks from FreeSurfer directory (by default stored in the surf/ folder). It returns a .csv file for each individual.
  
-	[MIND_long_02_COMBATLS.R](Code/MIND_long_02_COMBATLS.R) - applies the COMBATLS harmonization method to MIND networks.
  
-	[MIND_long_03_degree_and_edges_PAFIP_COMBATLS.m](Code/MIND_long_03_degree_and_edges_PAFIP_COMBATLS.m) – calculates the MIND degrees for each HC and SSD individual.
    
-	[MIND_long_04_gradients.m](Code/MIND_long_04_gradients.m) – calculates the MIND gradients for each HC and SSD individual, and generates brain maps of cortical and subcrotical degrees, as well as cortical gradients.

```matlab
gm_ref = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');

```

- [MIND_long_05_longitudinal_dx_MIND.m](Code/MIND_long_05_longitudinal_dx_MIND.m) - performs multiple regression and linear mixed modeling (LMM) on MIND networks.

- [MIND_long_06_longitudinal_BPRS_MIND.m](Code/MIND_long_05_longitudinal_BPRS_MIND.m) - performs generalized linear modeling (GLM) and generalized LMM (GLMM) on MIND networks and psychiatric symptoms. 

-	[MIND_long_07_hierarchy_SCZ_brain_maps.m](Code/MIND_long_07_hierarchy_SCZ_brain_maps.m) – generates the brain maps of cortical hierarchy and SCZ epicenters, and their correlation (co-localization) with MIND associations.

-	[MIND_long_08_cortical_MIND_association_maps.m](Code/MIND_long_08_cortical_MIND_association_maps.m) – generates the regional brain maps of cortical MIND associations.

-	[MIND_long_09_subcortical_MIND_association_maps.m](Code/MIND_long_09_subcortical_MIND_association_maps.m) – generates the regional brain maps of subcortical MIND associations.

  
### **Function calls**

This section contains the functions that are essential for running the scripts but must not be executed.

-	[mix_dx.m](Code/mix_dx.m) – creates randomized groups by mixing patients with different diagnoses or group membership. It is called by [MIND_02_degree_and_edges_PAFIP.m](Code/MIND_02_degree_and_edges_PAFIP.m) and [MIND_05_maturational_features.m](Code/MIND_05_maturational_features.m).

-	[regional_brainmap_representation.R](Code/regional_brainmap_representation.R) - generates regional brain maps from .csv files.

-	[regional_brainmap_representation_borders.R](Code/regional_brainmap_representation_borders.R) - generates from .csv files regional brain maps highlighting the significant regions.


## **License**

This project is licensed under the terms of the [GNU General Public License v3.0 license](LICENSE).


## **Cite**
If you use this software, please cite the following paper and software:

- García-San-Martín, N., Bethlehem, R.A., Sebenius, I. et al. Long-term morphometric similarity gradients relate to cortical hierarchy and psychiatric symptoms in schizophrenia. medRxiv (2026). https://doi.org/
  
- Natalia-García-San-Martín. NeuroimagingBrainNetworks/HierarchyLongitudinalGradientsPsychosis: v1.0.0-alpha. Zenodo https://doi.org/ (2026).
