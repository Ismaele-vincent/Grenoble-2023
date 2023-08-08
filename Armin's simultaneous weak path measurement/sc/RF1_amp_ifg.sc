START	standard	


SC	6229/ao0	Amps	1.73
G	GPIB1::8	fn1=SIN	*am1=0.01	fr1=62500	ph1=0	of1=0
G	GPIB1::8	start1


G	GPIB1::8	*am1=0.5
CALL	ifg_30s.sc
G	GPIB1::8	*am1=0.70
CALL	ifg_30s.sc
G	GPIB1::8	*am1=0.75
CALL	ifg_30s.sc
G	GPIB1::8	*am1=0.80
CALL	ifg_30s.sc
G	GPIB1::8	*am1=1.0
CALL	ifg_30s.sc

E
END
