START	standard	
CALL	RobotZ.sc
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
;I	PXIDAQ2/ao3	Amps	0.00000	0
;SC	6229/ao1	Amps	0.0
;C	time	120.000000
SC	6229/ao1	Amps	-0.5

I	PXIDAQ2/ao3	Amps	-1.17500	0
C	time	120.000000

I	PXIDAQ2/ao3	Amps	1.17500	0
C	time	120.000000

I	PXIDAQ2/ao3	Amps	0.00000	0
SC	6229/ao1	Amps	0.0
E
END
