START	standard		

CALL	setRF2.sc


F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
I	PXIDAQ2/ao1	Amps	1.000000
CALL	ifg_10s.sc
I	PXIDAQ2/ao1	Amps	1.1000
CALL	ifg_10s.sc
I	PXIDAQ2/ao1	Amps	1.2000
CALL	ifg_10s.sc
I	PXIDAQ2/ao1	Amps	1.30000
CALL	ifg_10s.sc

G	GPIB1::8	stop1
I	PXIDAQ2/ao1	Amps	0.00000
E
END
