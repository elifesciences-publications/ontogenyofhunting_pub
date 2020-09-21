## Ontogeny of Hunting in Larval Zebrafish  
Konstantinos Lagogiannis 2019 (c o s t a s l a g [at] gmail.com)
Example analysis code used to replicate statistical inference on zebrafish behaviour shown on manuscript: https://doi.org/10.7554/eLife.55119

*Note*: This code repository may change in response to feedback.
The following files contain example code in [R](https://www.r-project.org/ "R, *free* software environment for statistical computing and graphics"), along with the data replicating the analysis of zebrafish behaviour in the article *"Learning steers the ontogeny of an efficient hunting sequence in zebrafish larva."* , by K. Lagogiannis, G. Diana and M. P Meyer, of King's College London.

- stat_3DLarvaGroupBehaviour.r
  * Hierarchical model of revealing the typical hunting behaviour of each rearing group. 
- stat_FishDev-Length.r
  * Gaussian models used for comparing larval lengths between rearing groups.
- stat_ClusterCaptureBouts.r
  * The mixture of 2 Gaussians used to classify capture bouts based on speed and distance from prey.
- stat_HuntEfficiency.r
  * Model file showing how hunt rates and capture efficiency data was modeled jointly, and how the distribution of consumption of each group was estimated (Fig. 3)
 
The above files make use of the data in the dat subfolder. The data is preprocessed and was extracted from high-speed video recording that were analysed using a tracker I developed, which is not included here. The outcome of hunting, which is used to estimate the capture efficiency of each larva, was obtained by manually scoring process that was blind to rearing group.

