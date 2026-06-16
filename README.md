## Ancestral State Reconstruction Using Dollo Parsimony
Orthogroup (OG) presence/absence was reconstructed across the phylogeny using a custom Dollo parsimony implementation in R, applied to a binary species × OG character matrix (OrthoFinder table OrthoGroups.GeneCount.tsv) and a rooted phylogeny. Under Dollo parsimony, each character (OG) can be gained only once but lost independently multiple times. For each OG, the algorithm first identified tips scored as present, then traversed the tree upward to infer ancestral states. Direct parent nodes of present-state tips were marked present, followed by all internal nodes whose descendant sets included at least one present tip. The MRCA of all present-state nodes was then computed and added to the present list if applicable, along with any intermediate nodes connecting it to already-identified present nodes. Finally, the root was assigned the present state only if both of its immediate descendant subtrees contained at least one present tip, enforcing the single-gain assumption. The procedure was applied across all OGs, producing a node × OG matrix of inferred ancestral states. Each node of the tree was then classified by comparing inferred states to its immediate ancestral nodes: labeled as gain (absent → present), loss (present → absent), inherited, or nonexistent. To evaluate which topology required significantly fewer OG losses, pairwise Fisher's exact tests were applied to a 2×2 contingency table of total retained versus lost OGs under each tree hypothesis (using pairwise_fisher_test from the rstatix R package), with p-values adjusted for multiple comparisons. Effect sizes were quantified as the difference in loss-to-retention ratios between each pair of topologies. The topology requiring significantly fewer losses was interpreted as the better-supported hypothesis under Dollo parsimony.

## Analysis in:
Kulkarni et al. 2026. No silver bullet: Patterns of macrosynteny recapitulate systemic conflicts in the higher-level relationships of the arachinds. In prep.

## Files:
**Dollo.R** 
- Ancestral state inference of orthogroups under Dollo parsimony

**Arachnids.Rmd** 
- Tabulates ASR output from Dollo.R and compares loss-to-retention dynamics among competing topologies.
