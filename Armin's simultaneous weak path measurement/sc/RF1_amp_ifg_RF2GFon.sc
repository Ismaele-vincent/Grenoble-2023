START	standard	

CALL	setRF1_pi_contin_bothGF.sc
;SC	6229/ao0	Amps	1.73
;G	GPIB1::8	fn1=SIN	*am1=0.01	fr1=62500	ph1=0	of1=0
;G	GPIB1::8	start1

;I	PXIDAQ2/ao1	Amps	1.180000

G	GPIB1::8	*am1=0.50
CALL	ifg_30s.sc
G	GPIB1::8	*am1=0.60
CALL	ifg_30s.sc
G	GPIB1::8	*am1=0.65
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
