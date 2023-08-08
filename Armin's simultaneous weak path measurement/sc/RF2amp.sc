START	standard	
CALL	RobotZ.sc
G	GPIB1::8	fn2=SIN	*am2=0.01	fr2=60000	ph2=0	of2=0
SC	PXIDAQ2/ao1	Amps	0.0
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
C	time	10.000000
G	GPIB1::8	start2
SC	PXIDAQ2/ao1	Amps	1.70
G	GPIB1::8	*am2=0.01
C	time	10.000000
G	GPIB1::8	*am2=0.100
C	time	10.000000
G	GPIB1::8	*am2=0.200
C	time	10.000000
G	GPIB1::8	*am2=0.300
C	time	10.000000
G	GPIB1::8	*am2=0.400
C	time	10.000000
G	GPIB1::8	*am2=0.500
C	time	10.000000
G	GPIB1::8	*am2=0.600
C	time	10.000000
G	GPIB1::8	*am2=0.700
C	time	10.000000
G	GPIB1::8	*am2=0.800
C	time	10.000000
G	GPIB1::8	*am2=0.900
C	time	10.000000
G	GPIB1::8	*am2=1.000
C	time	10.000000
G	GPIB1::8	*am2=1.100
C	time	10.000000
G	GPIB1::8	*am2=1.200
C	time	10.000000
G	GPIB1::8	*am2=1.300
C	time	10.000000
G	GPIB1::8	stop2
SC	PXIDAQ2/ao1	Amps	0
E
END
