# 2020-07-09 Gsp1 DMS Plan PROGRESS

## ROSETTA DDG
1. Compile benchmark datasets
    - Download Kellogg and Protherm datasets from the kortemme lab website
      - git cloned Kortemme-Lab/ddg, the datasets are in
        - KELLOGG: in guybrush:ddg_benchmark/ddg/input/csv/kellogg.csv
        - PROTHERM: in guybrush:ddg_benchmark/ddg/input/csv/curatedprotherm.csv
    - Note the best score function as those results should be used in analysis
      - talaris 2013 (scaled version, but nonscaled is good too) for the Kellogg dataset
      - talaris 2013 for the Protherm dataset
    - Get Brandon's dataset. We already have Hahnboem's dataset (he emailed it)
      - supp table 1 in the preprint but not posted on bioRxiv, so downlaoded from slack where brandon shared it:
      - now in ddg_benchmark/datasets/balanced_benchmark.csv
      - Hahnboem's was saved in ddg_benchmark/datasets/ref2015.Lizddg.txt
    - Write code that will read in each dataset and have the list of mutants to include so that we can compute correlations based on dataset subsets and compare our results with published results (that use those subsets)
    - Get mutfiles from Shane and Kyle's protocol capture
2. Run cartesian ddG on all structures in benchmark
    - Rosetta version: use master after Frank and Brandon's recent commits fixing the relax repacking behavior are merged
    - Use the flags from Park et al 2016 JCTC for the preminimization and the ddG run, except -beta:cart should NOT be included because it is standard now (according to Hahnboem)
    - Pre-minimize all structures that appear in any dataset
    - Run ddG with all structures
3. Analyze `cartesian_ddg` benchmark
    - Compute correlation with experimental values for the full dataset, and for all of the subsets necessary to compare performance with published runs
    - Can also do the classifier task (destabilizing if ddg < -1, stabilizing if ddg > 1)
    - If performance is good (>0.6 correlation, and hopefully at least as good as Park et al 2016), then use the protocol on Gsp1
    - If performance is not good, ask Hahnboem again for the code version used in Park et al 2016, and troubleshoot on that exact replicate run
4. Run `cartesian_ddg` on Gsp1 (3m1i, residues 10-183)
    - preminimize, check to make sure side chains are not repacked (or if they are, that behavior is expected)
    - make sure WT and MUT structures are written out for each run, so we can compare WT noise from run and look for sources of noise in MUT calculations
    - run full saturation mutagenesis
5. Simple Ala scan (`mutate_residue` + minimize) 
    - Write the method (RosettaScripts or pyrosetta)
    - Run this method on benchmark set (preminimized with `cartesian_ddg`)
    - Compute correlation with experimental data, compare to `cartesian_ddg` results - can we use the simple Ala scan method for our purposes?
    - Run method on Gsp1
6. Analysis:
    - Compare Gsp1 results from 3 methods: `cartesian_ddg`, simple Ala scan, foldX
    - Correlate each of three methods with experimental EMPIRIC data
    - Can compute SASA vs. fitness and SASA vs. ddG metrics to see if the ddG matches SASA, which we expect it too
    - May try excluding residues from analysis if in van der waals contact with nucleotide (likely to have deleterious effects)
    - A good quality check is all heavy atom RMSD for the three lowest energy structures of each mutant (as ddG is calculated as mean for these 3 structures)
    - Predictability of experimental data: 
7. Extra ddG things we can do:
    - compare yeast-to-human mutations against the whole distribution to inform if we can use human structures for ddG
    - can also just exclude residues from analysis if in van der waals contact with non-conserved residues
    - if we use human structures: run `cartesian_dd`, simple Ala scan, and/or `flex_ddg` on:
      - open and closed monomer structures
      - GAP and GEF structures
      - All complexes
    - Can run ddG methods on ubiquitin, and see how predictive it is for the experimental data, as an extra confirmation of the method being appropriate for Gsp1 analysis

## SASA COMPUTATION
- Recompute SASA (burial) for the structures in PyRosetta or RosettaScripts
- Recompute interface SASA for the interface residues
  - "NONE" needs to be split up into "protein core" and "surface (non-interface)". the latter may be relatively small n
  - definitely combine rim/support, consider combining with core
- Send David dataframe of interface positions and burial
- Make a scatterplot of distance from nucleotide vs. fitness and compare to SASA (burial) vs. fitness

## COMPARISON HISTOGRAMS
1. Make histograms with new binnings so we can look at them all together
2. Make stacked barplots, where "in van der waals contact with nucleotide" is shaded
3. Comparisons:
    - Functional loops are enriched compared to other loops
      - Do we need helices vs. loops vs. sheets plot, perhaps distinguishing non-functional vs. functional loops, and core vs. surface
    - C-terminal linker with C-terminal helix
      - There is a deleterious proline in the helix, is it an n1 proline? Highlight this.
    - Protein core vs. surface exposed (non-interface) vs. surface exposed (interface)
      - We expect interface core to be less important for Gsp1 than one would think for protein interaction interfaces in general
    - Effector lobe vs. allosteric lobe
      - note that the new binning loses many beneficial mutants
      - is effector lobe being worse only due t functional loops? Compare helices vs. sheets vs. func loops vs. non-func loops between lobes
      - what can we really say about this? talk to tanja
4. Come up with a single metric to summarize interface conclusions across many partners. Also could use this metric for the other comparison categories.
    - One metric for interfaces could be number drop-outs in core of interface divided by number of drop-outs in interface (so including rim and support)
5. Make a GEF interface plot:
    - Stacked barplot showing GEF interface residues with those residues shared with GAP shaded
    - other categories are all other surface residues in an interface, all other surface residues not in an interface, and protein core residues
    - can consider doing this for other partners if it makes sense

## GSP1 INTRODUCTION WRITING
- What Gsp1 does
- Difficulty of studying Gsp1 as an essential gene
- Studying point mutations to understand function, cite our paper

## ALLOSTERIC LOBE SECTION
If we include this section:
- For the structure figure, show side chains of beneficial mutations, looking like the contiguous blobs from Rama Ranganathan's lab's papers
- Look at conservation of allosteric lobe ineractions
  - align GEFs of all GTPases, perhaps Carla Mattos has a paper on this
  - Look at structural variation of the lobe across all GTPases

## GTPASE CONSERVATION
- Look into Alfonso Valencia's papers
- How conserved is S. cerevisiae Gsp1 compared to other yeast homologs?
- How conserved are other S.c. GTPases compared to their other yeast homologs?
- How conserved is S.c. Gsp1 compared to other S.c. GTPases?

## RAS DMS COMPARISON
- Run structural alignment, to find corresponding residues between HRas and Gsp1
  - which positions behave most similarly or most differently between the two experiments?
- NOTE: Tina says the Kuriyan lab is currently re-doing these experiments, as the deletion of the last half of alpha helix 5 (6?), as was done in Pradeep's paper, seems to have destabilized the allosteric lobe mutants that were screened as deleterious (suggesting they may be false positives)

## DOMINANT NEGATIVE ANALYSIS
- Can look into gene dosage datasets from Weissman and Kirschner labs
- Is Gsp1 conserved because mutants are dominant negative with WT? Look at Ran T24N and Ran Q71L cell biology papers