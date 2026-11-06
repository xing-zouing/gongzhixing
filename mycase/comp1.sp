subckt   comp1   GND  VDD CK  von vop

M6 (net22 CK VDD VDD) pch l=2.14u w=4.28u multi=1 nf=4
M5 (von CK VDD VDD) pch l=2.14u w=4.28u multi=1 nf=4
M4 (net19 CK VDD VDD) pch l=2.14u w=4.28u multi=1 nf=4
M3 (vop CK VDD VDD) pch l=2.14u w=4.28u multi=1 nf=4
M1 (von vop VDD VDD) pch l=2.33u w=4.7u multi=1 nf=5
M2 (vop von VDD VDD) pch l=2.33u w=3.55u multi=1 nf=5
M10 (net5 CK GND GND) nch_lvt l=1.5u w=1.5u multi=1 nf=1
M9 (net22 vip net5 GND) nch_lvt l=1.5u w=1u multi=4 nf=1
M8 (net19 vin net5 GND) nch_lvt l=1.5u w=1u multi=4 nf=1
M7 (von vop net22 GND) nch_lvt l=1.5u w=1u multi=4 nf=1
M0 (vop von net19 GND) nch_lvt l=1.5u w=1u multi=4 nf=1

ends  comp1
