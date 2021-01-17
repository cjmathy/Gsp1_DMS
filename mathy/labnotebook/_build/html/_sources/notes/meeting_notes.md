# Meeting notes for DMS project

## 2020-05-04 David
- is codon signal consistent across experiments for beneficial mutants (i.e. likely that they really are  beneficial
- what if instead the Z-score is defined based on the STOP-like => that would likely get rid of "beneficial", which would mean GAPs are STOP-like and GEFs are WT-like (rahter than beneficial)

- can we use the fitness data prospectively predict the class of mutants (GAP-perturbed or GEF-perturbed)
- would likely throw out the mutants that are predicted to destabilize first
- then maybe look at ddGs of mutations in interface positions and the active site to see if that makes sense with the fitness
- ultimately, Dan's big thing is trying **to make the fitness scores predictive?**
  - can the fitness scores be argued to predict the GAP and GEF mutations in a meaningful way?

## 2020-07-21 Tanja
- Tanja says to run cartesian_ddg without proline flag

## 2020-08-27 with Tina
- about CONSURF analysis
- More negative score is more conserved
- CONSURF bin: higher is more conserved
- First distinction: effector lobe, allosteric lobe, tail
- Second distinction:
    - Structural core: always buried
    - Interface core: in at least one interface core
    - Surface: always exposed
    - mixed: sometimes interface sometimes surface

## 2020-08-28 with Tanja

- definitely the nucleotide changing conformations leading to the three different clusters of total_score
- this is OK because the ddg correlation with total score is consistent, i.e. the slope between total_score mutant and ddg is the same, even though the intercept is different. So it just reflects an energy difference between simulations that should be accounted for when subtracting MUT-WT
- plot fitness vs. ddg, showing foldx and rosetta in different colors. Also make a plot for the subset of points for which foldx and rosetta disagree

## 2020-10-30 with Tanja
- Tanja will email Dan to check in
- Need dominant negatives kinetics data to know how paper order will go
- agreed that being careful about aggregation is important
- suggested inserting standards in between runs on analytical when running a lot of samples
- re-run ddG with another structure or two. 3gj0 for Ran GDP-state and maybe another GTP-bound Ran structure

## 2020-11-23 thesis committee meeting
- from Martin: could we further understand which core mutations are permissive or tolerant by looking at conservation in yeast? what do the permissive core mutations have to do with the fact that yeast may be growing at a different temperature than the one in the assay?
