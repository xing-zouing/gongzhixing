subckt   OTA2_my    GND  VDD vb OUT VINN VINP

M12 (net17 vb GND GND) nch_lvt l=4u w=3u multi=1 nf=1

M0 (net22 VINP net19 GND) nch_lvt l=1.5u w=2u multi=4 nf=10
M7 (net24 VINN net19 GND) nch_lvt l=1.5u w=2u multi=4 nf=10
M6 (net21 net24 VDD VDD) pch l=10u w=3u multi=1 nf=1
M8 (OUT net22 VDD VDD) pch l=10u w=3u multi=1 nf=1
M10 (net21 net21 net17 GND) nch_lvt l=2u w=4u multi=1 nf=6
M11 (OUT net21 net17 GND) nch_lvt l=2u w=4u multi=1 nf=6

M9 (net19 vb GND GND) nch_lvt l=4u w=3u multi=1 nf=1

M3 (net24 net24 VDD VDD) pch l=5u w=2.5u multi=1 nf=1
M4 (net22 net22 VDD VDD) pch l=5u w=2.5u multi=1 nf=1

M2 (net22 net24 VDD VDD) pch l=5u w=2.5u multi=1 nf=1
M5 (net24 net22 VDD VDD) pch l=5u w=2.5u multi=1 nf=1

C1 (net24 net21) cfmom nr=6 l=1.5u w=1u s=70n stm=1 spm=3 multi=1 ftip=140.0n dmflag=0
C2 (net22 OUT) cfmom nr=6 l=1.5u w=1u s=70n stm=1 spm=3 multi=1  ftip=140.0n  dmflag=0


ends  OTA2_my
