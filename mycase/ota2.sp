subckt ota2  GND IBIAS VCM VDD VIM VIP VOM VOP

M0 (net43 net43 GND GND) nch_lvt l=3u w=2u multi=4 nf=1
M1 (net46 net43 GND GND) nch_lvt l=3u w=2u multi=4 nf=1

M2 (out net46 GND GND) nch_lvt l=2u w=15u multi=1 nf=1
C0 (out net46) cfmom nr=40 l=3u w=8u s=100n stm=1 spm=3 multi=1 ftip=140.0n

M8 (net46 vinn net45 VDD) pch l=3u w=2u multi=4 nf=1
M9 (net43 vinp net45 VDD) pch l=3u w=2u multi=4 nf=1

M10 (net45 vb VDD VDD) pch l=6u w=2u multi=1 nf=1
M11 (out vb VDD VDD) pch l=6u w=2u multi=1 nf=1

ends ota2
