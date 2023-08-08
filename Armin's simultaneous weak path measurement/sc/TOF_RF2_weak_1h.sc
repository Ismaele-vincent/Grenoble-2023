START	standard
;R	RF GF 2                 
CALL	nullsetzer.sc
;CALL	setRF2_weak.sc
CALL	setDC_pi_half.sc
SC	PXIDAQ2/ao2	Amps	1.0709
R	20 deg Spin Rot RF 2
G	GPIB1::8	fn2=SIN	*am2=0.111	fr2=55500	mo2=BURS	bt2=16	ph2=0	of2=0
G	GPIB1::8	start2
;********Start
;MS	S18.auxiliary.pathblocker	1	pos	absolute	0.00000
F	D:\data\Cycle 187\exp_CRG-2716\rawdata\sc\TOF\20deg_RF2\*.dat
TOF	4	time=3600	period=304µ	binwidth=1µ	delay=35µ
E
END
