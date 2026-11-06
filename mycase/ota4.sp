subckt   OTA3_my    VSS  VDD VB VOUT VN VP

M1 (VOUT VB VDD VDD) pch l=1.5u w=4.5u multi=1 nf=8
M2 (VB VB VDD VDD) pch l=1.5u w=4.5u multi=1 nf=8

M3 (net42 VN net45 net45) pch l=2u w=1.5u multi=4 nf=4
M4 (net46 VP net45 net45) pch l=2u w=1.5u multi=4 nf=4

M0 (net45 VB VDD VDD) pch  l=2u w=4.5u multi=2 nf=8

C1 (net46 VOUT) cfmom_2t nr=20 l=2.5u w=4u s=100n stm=1 spm=3 multi=1  ftip=140.0n

M5 (net46 net42 VSS VSS) nch_lvt l=2u w=1.5u multi=4 nf=4
M6 (net42 net42 VSS VSS) nch_lvt l=2u w=1.5u multi=4 nf=4

M7 (VOUT net46 VSS VSS) nch_lvt l=3u w=4.5u multi=1 nf=8

ends  OTA3_my
