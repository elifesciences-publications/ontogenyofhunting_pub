## Ontogeny of Hunting in Larval Zebrafish  
Konstantinos Lagogiannis 2019 (c o s t a s l a g [at] gmail.com)
Example analysis code used to replicate statistical inference on zebrafish behaviour shown on manuscript : 
[BioArxiv Pre-Print Available](https://www.biorxiv.org/content/10.1101/2019.12.19.883157v2 "BioArxiv Pre-Print available")

*Note*: This code repository may change in response to feedback.
The following files contain example code and data replicating the analysis of behaviour in the article *"Learning steers the ontogeny of an efficient hunting sequence in zebrafish larva."* , by K. Lagogiannis, G. Diana and M. P Meyer, of King's College London.

- stat_3DLarvaGroupBehaviour.r
  * Hierarchical model of revealing the typical hunting behaviour of each rearing group. 
- stat_FishDev-Length.r
  * Gaussian models used for comparing larval lengths between rearing groups.
- stat_ClusterCaptureBouts.r
  * The mixture of 2 Gaussians used to classify capture bouts based on speed and distance from prey.
- stat_HuntEfficiency.r
  * Model file showing how hunt rates and capture efficiency data was modeled jointly, and how the distribution of consumption of each group was estimated (Fig. 3)
 
The above files make use of the data in the dat subfolder. The data is preprocessed and was extracted from high-speed video recording that were analysed using a tracker I developed, which is not included here. The outcome of hunting, which is used to estimate the capture efficiency of each larva, was obtained by manually scoring process that was blind to rearing group.

