# 2020-06-08 Tina on EMAP SASA

- I now rewrote the script to output a more general file called Data/Gsp1_interfaces_SASA_and_conservation.txt
- The script that makes that file is in `Scripts/complex_structure_analyses/calculate_interfaces_and_interface_conservation_for_all_chains.R.`
- this script only works because Gsp1 is always chain A and it needs a different index file than before
- which is in the folder with the script
- then there is a downstream script that takes that file and calculates stuff
- I needed that to make the EDF1 for reviewer 1
- cos now I needed both sides of the interface
- the downstream script is interfaces_and_complexes_analyses.R
- but I can't remember if that one also makes EDF1 last panel
- let me check
- if yes, we can kind of abandon the interface_residues_and_mutations.txt file
- or we can just parse the Data/Gsp1_interfaces_SASA_and_conservation.txt to have only the Gsp1 interfaces and name that interface_residues_and_mutations.txt