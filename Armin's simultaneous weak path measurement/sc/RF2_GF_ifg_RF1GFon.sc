START	standard		

CALL	nullsetzer.sc
CALL	setRF2_pi_contin.sc
SC	6229/ao0	Amps	1.80

F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
I	PXIDAQ2/ao1	Amps	0.200000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	0.700000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	1.000000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	1.1000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	1.2000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	1.50000
CALL	ifg_30s.sc
I	PXIDAQ2/ao1	Amps	2.000000
CALL	ifg_30s.sc

G	GPIB1::8	stop1
I	PXIDAQ2/ao1	Amps	0.00000
E
END
