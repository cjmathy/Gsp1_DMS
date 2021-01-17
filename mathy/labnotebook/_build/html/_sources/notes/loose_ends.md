# Loose ends - some general things to think about for DMS 

## Questions that our manuscript needs to address:
- how do we explain neutral mutations that in functional loops?
- how do we explain interface mutations that are neutral?
- how do we explain core mutations that are neutral

## ddG scatterplots
- neutral mutatioons with 10 REU score and above, look at lennard jones repulsion and dunbrack, does that explain high scores?
- need to examine a couple of regions of outliers for ddGs
  - negative ddG but negative fitness is perhaps enriched in interface positions, i.e. bad for fitness but  
- are negative ddGs but negative fitnesses in interface?
- what are the residues with positive ddG but neutral fitness
- remove glycines and prolines?
- two questions: how many core mutations bad in ddG and fitness? contacting nucleotide and bad in fitness?

## ddG rigor
- re-run with another structure to see if finding is robust?

## per position correlations of ddG and fitness
- look at number of points in correlation, may be spurious
- maybe see some interesting backbone conformational encoding (look at one-body terms)

## One possible grided summary figure:
- groups of mutatons by bin along the x-axis
- functional group in y-axis
- each grid cell is a plot comparing the two groups of the y-axis label for that subset of mutants
- one half is "traditional" y-axis groups => core vs surface, for example
- other half is "expert" or "informed" groups that better explain => which interface it's in, allosteric vs effector lobe, functional loops vs. other

## Points to think about for manuscript
- surprising that there are beneficial mutations
  - because partner interfaces are sites of regulation and the mutations mimic that?
  - perhaps losing binding to a subset of partners is good (in this assays) => David is looking for gene sets that show fitness increase over expectation when deleted in yeast
- all other DMS does biochem in cell (is this still true)
- This paper is just fitness, and we are considering mutations in background of all interactions => contextualize results with all interactions/functions

## EMAP average fitness (not urgent)
- need to check code:
  - what is being measured
  - are batch effects being handled
- is this worth it? 


  
