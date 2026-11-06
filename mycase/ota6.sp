///ota6(ota4)
subckt   ota6    VSS  VDD VREF OUTP OUTM   VBIAS_P  INM INP

M22 (vbias_n VBIAS_P VDD VDD) pch l=320n w=320n multi=1 nf=1
M21 (VBIAS_P VBIAS_P VDD VDD) pch l=320.0n w=320.0n multi=1 nf=1
M19 (net103 VSS net88 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M0 (net104 VSS net86 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M16 (vcmfb OUTM net93 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M15 (net92 VREF net91 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M14 (net92 VREF net93 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M13 (vcmfb OUTP net91 VDD) pch l=200.0n w=200.0n multi=2 nf=1
M12 (net93 VBIAS_P VDD VDD) pch l=320n w=320n multi=1 nf=1
M10 (net86 VBIAS_P VDD VDD) pch l=150.0n w=150.0n multi=2 nf=1
M9 (net88 VBIAS_P VDD VDD) pch l=150.0n w=150.0n multi=2 nf=1
M8 (OUTM net88 VDD VDD) pch l=320n w=320n multi=1 nf=1
M7 (OUTP net86 VDD VDD) pch l=320n w=320n multi=1 nf=1
M11 (net91 VBIAS_P VDD VDD) pch l=320n w=320n multi=1 nf=1
M20 (vbias_n vbias_n VSS VSS) nch_lvt l=320n w=320n multi=1 nf=1
M18 (net92 net92 VSS VSS) nch_lvt l=320n w=320n multi=1 nf=1
M17 (vcmfb vcmfb VSS VSS) nch_lvt l=320n w=320n multi=1 nf=1
M6 (OUTP vbias_n VSS VSS) nch_lvt l=320n w=320n multi=1 nf=1
M5 (OUTM vbias_n VSS VSS) nch_lvt l=320n w=320n multi=1 nf=1
M4 (net87 vcmfb VSS VSS) nch_lvt l=320.0n w=500n multi=1 nf=1
M3 (net88 INM net87 VSS) nch_lvt l=100.0n w=100.0n multi=4 nf=4
M2 (net86 INP net87 VSS) nch_lvt l=100.0n w=100.0n multi=4 nf=4
M1 (net87 vbias_n VSS VSS) nch_lvt l=150.0n w=500n multi=1 nf=1
C1 (OUTP net104) cfmom  nr=32 nh=32 l=1.2u w=800n s=70n stm=1 spm=3 multi=1
C0 (OUTM net103) cfmom  nr=32 nh=32 l=1.2u w=800n s=70n stm=1 spm=3 multi=1
ends  ota6
