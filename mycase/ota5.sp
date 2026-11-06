subckt   ota5    VSS  VDD VREF OUTP OUTM   VBIAS_P  INM INP

M16 (vcmfb OUTM net87 VDD) pch l=700n w=800n multi=2 nf=1
M15 (net86 VREF net85 VDD) pch l=700n w=800n multi=2 nf=1
M14 (net86 VREF net87 VDD) pch l=700n w=800n multi=2 nf=1
M13 (vcmfb OUTP net85 VDD) pch l=700n w=800n multi=2 nf=1

M12 (net87 VBIAS_P VDD VDD) pch l=1u w=1.2u multi=1 nf=1

M10 (net81 VBIAS_P VDD VDD) pch l=700n w=800n multi=2 nf=1
M9 (net83 VBIAS_P VDD VDD) pch l=700n w=800n multi=2 nf=1

M8 (OUTM net83 VDD VDD) pch l=1u w=1.2u multi=1 nf=1
M7 (OUTP net81 VDD VDD) pch l=1u w=1.2u multi=1 nf=1
M11 (net85 VBIAS_P VDD VDD) pch l=1u w=1u multim=1 nf=1
M18 (net86 net86 VSS VSS) nch_lvt l=1u w=1u multi=1 nf=1
M17 (vcmfb vcmfb VSS VSS) nch_lvt l=1u w=1u multi=1 nf=1
M6 (OUTP vbias_n VSS VSS) nch_lvt l=1u w=1u multi=1 nf=1
M5 (OUTM vbias_n VSS VSS) nch_lvt l=1u w=1u multi=1 nf=1
M4 (net82 vcmfb VSS VSS) nch_lvt l=1u w=1u multi=1 nf=1

M3 (net83 INM net82 VSS) nch_lvt l=1u w=1.2u multi=4 nf=4
M2 (net81 INP net82 VSS) nch_lvt l=1u w=1.2u multi=4 nf=4

M1 (net82 vbias_n VSS VSS) nch_lvt l=1u w=1.2u multi=1 nf=1

ends  ota5
