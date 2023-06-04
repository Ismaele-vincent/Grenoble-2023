START	standard
G	GPIB1::8	fn1=SIN	*am1=0.01	fr1=62500	ph1=0	of1=0
G	GPIB1::8	start1
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
SC	AO_A1Z	Amps	0.0
G	GPIB1::8	*am1=0.00
G	GPIB1::8	stop1
C	time	10.000000
SC	AO_A1Z	Amps	1.50
G	GPIB1::8	start1
;G	GPIB1::8	*am1=0.100
;C	time	10.000000
G	GPIB1::8	*am1=0.200
C	time	10.000000
;G	GPIB1::8	*am1=0.300
;C	time	10.000000
G	GPIB1::8	*am1=0.400
C	time	10.000000
;G	GPIB1::8	*am1=0.500
;C	time	10.000000
G	GPIB1::8	*am1=0.600
C	time	10.000000
G	GPIB1::8	*am1=0.700
C	time	10.000000
G	GPIB1::8	*am1=0.750
C	time	10.000000
G	GPIB1::8	*am1=0.800
C	time	10.000000
G	GPIB1::8	*am1=0.850
C	time	10.000000
G	GPIB1::8	*am1=0.900
C	time	10.000000
G	GPIB1::8	*am1=1.000
C	time	10.000000
G	GPIB1::8	*am1=1.100
C	time	10.000000
G	GPIB1::8	*am1=1.200
C	time	10.000000
G	GPIB1::8	*am1=1.300
C	time	10.000000

G	GPIB1::8	stop1
SC	AO_A1Z	Amps	0
E
END
