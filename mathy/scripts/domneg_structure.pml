# last updated 2020-09-29

run ~/gdrive/gsp1_dms/mathy/scripts/set_ucsf_colors_in_pymol.pml

bg_color white

delete everything

load "~/gdrive/gsp1_dms/Data/pdbs_raw/1k5d.pdb"

create Ran, 1k5d and chain A and resi 1-213
create GNP, 1k5d and chain A and resi 1250
create Mg, 1k5d and chain A and resi 1251
create RanBP1, 1k5d and chain B and resi 22-167
create RanGAP, 1k5d and chain C and resi 2-386
delete 1k5d

color ucsf_navy2, Ran
color ucsf_yellow1, RanBP1
color ucsf_orange1, RanGAP
util.cbaw not polymer

select dominant_negatives, Ran and resi 17-20+22-24+26+30+33+41+42+44+48+49+52+57+65-70+73+75+76+79+91+94+96+99+123+125+130+150+151+153-155+157+158+161+181+181+184+185+186+190+200+204+209

util.cbac Ran and dominant_negatives
show sticks, Ran and dominant_negatives

set_view (\
    -0.174172312,    0.599935174,   -0.780852139,\
    -0.081910811,    0.781400144,    0.618628025,\
     0.981306493,    0.171710983,   -0.086957209,\
    -0.000241444,   -0.000040218, -249.231109619,\
    -0.728888988,   -3.126214027,  -27.023597717,\
  -212.465332031,  710.923828125,  -20.000000000 )
  
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3

png ~/gdrive/gsp1_dms/mathy/figures/domneg_structure_1k5d.png, width=20cm, dpi=300, ray=1
