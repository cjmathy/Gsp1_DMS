# 2020-06-09 Mutation info for David

## 31P NMR results
- Ran/Gsp1 bound to GTP shows 2 peaks for its gamma phosphate, corresponding to two discrete states in slow exchange (on NMR time scale). Upper limit at 5C is 200 per second, lower limit is 1e-3.
- State 2 corresponds to the hydrolytically compatible state. We think this for two reasons:
  1. Intrinsic and GAP-mediated hydrolysis are both generally faster for mutants that are more in state 2
  2. RanBP1 (yeast Yrb1) binding stabilizes state 2, catalyzes GAP-mediated hydrolysis, and is present in a solved ternary complex with the GAP
- WT is 87% state 2
- We consider four sites as allosteric sites, since they deviate significantly from WT in the 31P NMR:
  - T34:
    - T34E/Q = 0% state 2, T34A = 12%, T34L = 25%, and T34G = 64%
    - T34 is in the binding interface with RanBP1/Yrb1, so we can see how either a binding event or a mutation at this site can affect the active site
  - Q147: 
    - Q147E is 53.9% in state 2, so also shifted to be less active. This was a nice result because it is far from the T34 site, so we found two different allosteric sites that favor state 1 more, strengthening the idea that the active site is coupled to interaction interaces
    - This is in the karyopherin interfaces
  - Y157, H141 (and K132):
    - Y157A, H141R, and K132H are all only observed in state 2.
    - Y157 and H141 are in the karyopherin interface
    - K132H is in the karyopherin and GAP interface
    - we didn't refer to K132H as explicitly allosteric in the paper because it is in the GAP interface, meaning it's GAP-mediated hydrolysis (which would be predicted to increase) is decreased by the interface perturbation. I think it is, however, fair to call this interface coupling as well

## GEF and GAP results
- GAP:
  - The GAP activity of the mutants generally correlates with the NMR data, so decreased activity mutants are the T34's, D79'S, Q147E.
  - Also decreased: K132H (interface mutant), R78K (part of active site loop - switch II), and R112S (likely due to being destabilized - unable to get NMR data because of this)
  - Some R108's are slower in GAP, some are faster
- GEF:
  - All mutants except R108A decrease the GEF activity
  - especially bad ones are K101R and the R108's
  - R78K and R112S, are also bad in GEF

## Stability
- WT apparent Tm is 76. Unstable ones seem to be R112's (R112S - 71), H141's (63-66), and Y157A (63)