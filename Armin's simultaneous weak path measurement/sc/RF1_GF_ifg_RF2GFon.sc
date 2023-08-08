START	standard		

CALL	setRF1_pi_contin_bothGF.sc
I	PXIDAQ2/ao1	Amps	1.180000

F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks

SC	6229/ao0	Amps	1.65000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.72500
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.80000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.90000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.05000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.20000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.5000
CALL	ifg_30s.sc

G	GPIB1::8	stop1
SC	6229/ao0	Amps	0.00000
E
END
