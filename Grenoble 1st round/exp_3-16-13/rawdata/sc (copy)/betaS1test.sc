START	standard		
CALL	ifg_1s_3p.sc
MS	S18.axis2.phase shifter	1	deg	FIT	H	COS0
CALL	betaS1pi8test.sc
MS	S18.auxiliary.LinearBlue	1	mm	FIT	O	COS0
CALL	DC2Y.sc
CALL	ifg_1s_3p.sc
MS	S18.axis2.phase shifter	1	deg	FIT	H	COS60
CALL	betaS1pi8test.sc
END
