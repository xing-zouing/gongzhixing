subckt   Telescopic_Three_stage_flow INM INP OUTM OUTP VBN1 VDD VREF VSS

M60 (net057 net056 net058 VSS) nch l=1u w=2u multi=4 nf=4

M105 (net058 VBN VSS VSS) nch l=1u w=2u multi=1 nf=8

M93 (VCMFB3 VREF net058 VSS) nch l=1u w=2u multi=4 nf=4

M110 (VOM2 VCMFB2 VSS VSS) nch l=1u w=2u multi=1 nf=8

M75 (VOP2 VCMFB2 VSS VSS) nch l=1u w=2u multi=1 nf=8

M97 (VBP VBN VSS VSS) nch l=1u w=2u multi=1 nf=4

M24 (VBN1 VBN1 VBN VSS) nch l=1u w=2u multi=1 nf=10

M96 (VBP1 VBN VSS VSS) nch l=1u w=2u multi=1 nf=4

M21 (VBN VBN VSS VSS) nch l=1u w=2u multi=1 nf=4

M113 (OUTP VOM2 VSS VSS) nch l=1u w=2u multi=1 nf=16

M13 (OUTM VOP2 VSS VSS) nch l=1u w=2u multi=1 nf=16

M102 (VOM2 VBN VSS VSS) nch l=1u w=2u multi=1 nf=4

M18 (VOP2 VBN VSS VSS) nch l=1u w=2u multi=1 nf=4

M106 (VOP1 VBN1 net2 VSS) nch l=1u w=2u multi=1 nf=10


M4 (net1 VBN VSS VSS) nch l=1u w=2u  multi=1 nf=28

M3 (VOM1 VBN1 net4 VSS) nch l=1u w=2u multi=1 nf=10

M63 (net057 net057 VDD VDD) pch l=1u w=2u multi=1 nf=4

M112 (OUTP VCMFB3 VDD VDD) pch l=1u w=2u multi=1 nf=40

M73 (OUTM VCMFB3 VDD VDD) pch l=1u w=2u multi=1 nf=40

M114 (VCMFB3 net057 VDD VDD) pch l=1u w=2u multi=1 nf=4

M26 (VBP1 VBP1 VDD VDD) pch l=1u w=2u multi=1 nf=10

M90 (VBP VBP VDD VDD) pch l=1u w=2u multi=1 nf=4

M111 (OUTP VBP VDD VDD) pch l=1u w=2u multi=1 nf=40

M14 (OUTM VBP VDD VDD) pch l=1u w=2u multi=1 nf=40

M109 (VOM2 VOP1 VDD VDD) pch l=1u w=2u multi=1 nf=80

M9 (VOP2 VOM1 VDD VDD) pch l=1u w=2u multi=1 nf=80

M7 (net07 VCMFB1 VDD VDD) pch l=1u w=2u multi=1 nf=4


M6 (VOM1 VBP1 net07 VDD) pch l=1u w=2u multi=1 nf=10


M107 (VOP1 VBP1 net08 VDD) pch l=1u w=2u multi=1 nf=10

M108 (net08 VCMFB1 VDD VDD) pch l=1u w=2u multi=1 nf=4


C7 (VOP1 OUTM ) cfmom nr=50 l=6u w=2u s=70n stm=1 spm=6 multi=1

C6 (VOM1 OUTP ) cfmom nr=50 l=6u w=2u s=70n stm=1 spm=6 multi=1

C5 (VOM2 OUTM ) cfmom nr=88 l=8u w=2u s=70n stm=1 spm=6 multi=1

C1 (VOM2 OUTP ) cfmom nr=88 l=8u w=2u  s=8u stm=1 spm=6 multi=1

R8 (VCMFB2 VOM2  ) rppolywo_m l=10u w=6u multi=1

R4 (VOP1 VCMFB1 ) rppolywo_m l=12u w=4u multi=1

R9 (OUTP net056  ) rppolywo_m l=14u w=6u multi=1

R6q (VCMFB1 VOM1  ) rppolywo_m l=12u w=4u multi=1

R10 (net056 OUTM  ) rppolywo_m l=14u w=6u multi=1

R7 (VOP2 VCMFB2 ) rppolywo_m l=10u w=6u multi=1

M12 (net4 INP net1 VSS) nch l=2u w=1u multi=1 nf=10

M2 (net2 INM net1 VSS) nch l=2u w=1u multi=1 nf=10

ends Telescopic_Three_stage_flow

