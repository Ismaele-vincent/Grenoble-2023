START	standard
G	GPIB1::8	fn1=SIN	am1=0.75	fr1=62500	ph1=0	of1=0.00
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
G	GPIB1::8	stop1
I	AO_A1Z	Amps	0.000000	0
C	time	10.000000
G	GPIB1::8	start1
;I	AO_A1Z	Amps	0.500000	0
;C	time	10.000000
I	AO_A1Z	Amps	0.700000	0
C	time	10.000000
I	AO_A1Z	Amps	0.800000	0
C	time	10.000000
I	AO_A1Z	Amps	0.900000	0
C	time	10.000000
I	AO_A1Z	Amps	1.00000	0
C	time	10.000000
I	AO_A1Z	Amps	1.10000	0
C	time	10.000000
I	AO_A1Z	Amps	1.20000	0
C	time	10.000000
I	AO_A1Z	Amps	1.30000	0
C	time	10.000000
I	AO_A1Z	Amps	1.40000	0
C	time	10.000000
I	AO_A1Z	Amps	1.50000	0
C	time	10.000000
I	AO_A1Z	Amps	1.60000	0
C	time	10.000000
I	AO_A1Z	Amps	1.70000	0
C	time	10.000000
I	AO_A1Z	Amps	1.80000	0
C	time	10.000000
I	AO_A1Z	Amps	1.90000	0
C	time	10.000000
I	AO_A1Z	Amps	2.00000	0
C	time	10.000000
I	AO_A1Z	Amps	2.10000	0
C	time	10.000000
I	AO_A1Z	Amps	2.20000	0
C	time	10.000000
I	AO_A1Z	Amps	2.30000	0
C	time	10.000000

G	GPIB1::8	stop1
I	AO_A1Z	Amps	0.00000	0
E
END
