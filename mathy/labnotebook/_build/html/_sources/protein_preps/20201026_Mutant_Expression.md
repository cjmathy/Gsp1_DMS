# 2020-10-26 Gsp1 Mutant Expression

Goal: express pET28a+ mutants in BL21 cells. Tina previously prepped the plasmids and the competent BL21 cells

## 2020-10-26 Transformation

- Transforming BL21 cells with the following plasmids (transformed 8, want to make sure to have 6 to grow this week):
    - pE44_G35A
    - pE45_F54A
    - pE47_F54W
    - pE50_F163L
    - pE51_H32I
    - pE52_H32L
    - pE53_H32D
    - pE54_H32E
- Heat shocked, recovered 40 min at 37C, spun down 2m at 4000 rpm, plated on LB+Kan, overnight at 37C

## 2020-10-27 Inoculation

- All plates are covered in colonies
- Inoculated in ZY media + 50 mL NPS + 20 mL 5052 buffer + 1 mL 1M MgSO4 + 1 mL trace metals + 1 mL Kan. Started at 13h30
    - pE44_G35A
    - pE45_F54A
    - pE47_F54W
    - pE50_F163L
    - pE51_H32I
    - pE52_H32L
- Also inoculated one colony per plate into overnight cultures (5 mL LB+Kan) to do a colony PCR to confirm I transformed the correct plasmids

## 2020-10-28 Colony PCR

Running Colony PCR on the overnights from the transformed BL21s, to confirm I transformed the correct strains

- Diluted 1 uL of overnights into 25 uL of ddH2O
- Incubated at 95C for 5 min
- Spun down at max, 2 min
- Expected size: 725, using pET28A_6xHis-Gsp_WT as target

### Primers:

oligo id | name             | sequence              | dir | len | %GC | Tm |
---------|------------------|-----------------------|-----|-----|-----|----|
oCJM0046 | pET_Gsp1_CDS_fwd | CATCATCATCATCACAGCAGC | F   | 21  | 48  | 55 |
oCJM0047 | pET_Gsp1_CDS_rev | TGTCCACCAGTCATGCTAG   | R   | 19  | 53  | 55 |

### PCR recipe:
Reagent          | Volume (uL) 
---------------- | ------
2x Q5 Master Mix | 10
template         | 1
10 uM primers    | 1
H2O              | 8
**TOTAL**        | **20**

### PCR Protocol:
- 98C, 3 min
- Repeat 35x:
  - 98C, 30 sec
  - 52C, 30 sec
  - 72C, 1 min/kb rule of thumba
- 72C, 2 min

## 2020-10-29 Colony PCR continued
- Samples were at 4C overnight in the PCR machine.
- Re-ran PCR with 1 uL of PCR product as template
- Ran 1% DNA gel, 150 V for 20 min, with GeneRuler 100 bp DNA Ladder
  - Left to Right: PCR 1 (44, 45, 47, 50, 51, 52), PCR 2 (44, 45, 47, 50, 51, 52), Ladder
- Results at '../images/20201029_Gsp1_mutant_colony_PCRs_labeled.png'

<img src="../../images/20201029_Gsp1_mutant_colony_PCRs_labeled.png" class="centered" width="500px">

- First PCR worked fine. Second PCR didn't work, probably because I didn't resuspend the first PCR before taking template. Or because I took too much template?
- Diluted first PCR samples 1:2 in ddH2O to increase volume and sent to Quintara for sequencing to confirm mutants. Used oCJM0046 as a primer for sequencing (sent 100 uL of 10 uM).

## 2020-10-30 Sequencing results
- got back results, all are correct