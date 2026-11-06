subckt   OPA_my    GND  VDD VB VOUT VIN VIP

M7 (net103 net96 net98 GND) nch_lvt w=2u l=5u nf=1 multi=1
M6 (net103 net98 net116 GND) nch_lvt w=2u l=5u nf=1 multi=1
M15 (net98 net98 net117 GND) nch_lvt w=2u l=5u nf=1 multi=1
M8 (net96 net98 GND GND) nch_lvt w=2u l=5u nf=1 multi=1
R0 (GND net116) rppolywo l=5.5u w=20u series=6 segspace=250n
C0 (net111 VOUT) cfmom nr=60 l=8u w=5u s=70n stm=1 spm=6 multi=1 ftip=140.0n dmflag=0
M9 (net97 net97 VDD VDD) pch w=3.5u l=5u  nf=1 multi=1
M10 (net108 net103 VDD VDD) pch w=4u l=2u nf=1 multi=1
M11 (VOUT net103 VDD VDD) pch w=3.5u l=5u  nf=1 multi=1
M12 (net103 net103 VDD VDD) pch w=4u l=2u nf=1 multi=1
M13 (net98 net103 VDD VDD) pch w=1.2u l=1.2u nf=1 multi=4
M14 (net96 net98 net97 VDD) pch w=1.2u l=1.2u nf=1 multi=4
M5 (net117 GND GND GND) nch_lvt w=3.5u l=5u nf=1   multi=1
M2 (net105 net105 GND GND) nch_lvt w=3.5u l=5u nf=1 multi=1
M1 (net105 VIN net108 net108) pch w=3.5u l=5u nf=1 multi=1
M0 (net111 VIP net108 net108) pch w=3.5u l=5u nf=1 multi=1
M3 (VOUT net111 GND GND) nch_lvt w=4u l=1u nf=1 multi=1
M4 (net111 net105 GND GND) nch_lvt w=4u l=1u nf=1  multi=1
ends  OPA_my
