START	standard	


I	PXIDAQ2/ao1	Amps	1.0
G	GPIB1::8	fn2=SIN	*am2=0.01	fr2=55500	ph2=0	of2=0
G	GPIB1::8	start2

G	GPIB1::8	*am2=0.6
CALL	ifg_10s.sc
G	GPIB1::8	*am2=0.7
CALL	ifg_10s.sc
G	GPIB1::8	*am2=0.8
CALL	ifg_10s.sc
G	GPIB1::8	*am2=0.9
CALL	ifg_10s.sc
G	GPIB1::8	*am2=1.0
CALL	ifg_10s.sc
G	GPIB1::8	*am2=1.1
CALL	ifg_10s.sc
G	GPIB1::8	*am2=1.2
CALL	ifg_10s.sc
G	GPIB1::8	*am2=1.3
CALL	ifg_10s.sc



E
END
