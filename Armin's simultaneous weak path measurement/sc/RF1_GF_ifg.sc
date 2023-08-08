START	standard		

CALL	setRF1_cont.sc


F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\*.dat
R	room to enter remarks
SC	6229/ao0	Amps	1.30000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.40000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.50000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.60000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.70000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.75000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.80000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	1.90000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.0000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.10000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.20000
CALL	ifg_30s.sc
SC	6229/ao0	Amps	2.30000
CALL	ifg_30s.sc

G	GPIB1::8	stop1
SC	6229/ao0	Amps	0.00000
E
END
