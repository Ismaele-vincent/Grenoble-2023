START	standard
CALL	setDCpi2.sc
SC	AO_H	Amps	-0.011285
;CALL	ifg_test.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-1.00000
;CALL	beta_alpi8off_test.sc
CALL	gamma_test.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.91000
;CALL	beta_alpi8off_test.sc
CALL	gamma_test.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.82000
;CALL	beta_alpi8off_test.sc
END
