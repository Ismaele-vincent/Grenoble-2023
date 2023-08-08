START	standard	
G	GPIB1::8	fn2=SIN	*am2=0.5	fr2=55500	ph2=0	of2=0	mo2=BURST
G	GPIB1::8	start2
I	PXIDAQ2/ao1	Amps	1.1
END
