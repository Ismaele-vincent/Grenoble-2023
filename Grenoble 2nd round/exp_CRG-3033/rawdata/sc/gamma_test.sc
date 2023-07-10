START	standard
F	D:\data\Cycle 193\exp_CRG-3033\rawdata\sc\gamma_test.dat
R	room to enter remarks
CALL	setA1Zon.sc
I	AO_DC2Y	Amps	-0.600000	0
C	time	3
I	AO_DC2Y	Amps	-0.400000	0
C	time	3
I	AO_DC2Y	Amps	-0.200000	0
C	time	3
I	AO_DC2Y	Amps	0.000000	0
C	time	3
I	AO_DC2Y	Amps	0.200000	0
C	time	3
I	AO_DC2Y	Amps	0.400000	0
C	time	3
I	AO_DC2Y	Amps	0.600000	0
C	time	3
CALL	setA1Zoff.sc
E
CALL	setDC2ppi2.sc
END
