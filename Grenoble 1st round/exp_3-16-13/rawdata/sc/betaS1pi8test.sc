START	standard		
CALL	setS2Zoff.sc
M	S18.auxiliary.LinearBlue	mm	-0.010000	FALSE
F	D:\data\Cycle 192\exp_3-16-13\rawdata\sc\betaS1Pi8_ON.dat
F	D:\data\Cycle 192\exp_3-16-13\rawdata\sc\betaS1Pi8_OFF.dat
R	room to enter remarks
M	S18.auxiliary.LinearBlue	mm	0.000000	FALSE
I1	AO_S1Z	Amps	0.521000
C	time	2.000000
I0	AO_S1Z	Amps	0.000000
C	time	2.000000
M	S18.auxiliary.LinearBlue	mm	5.000000	FALSE
I0	AO_S1Z	Amps	0.000000
C	time	2.000000
I1	AO_S1Z	Amps	0.521000
C	time	2.000000
E
END
