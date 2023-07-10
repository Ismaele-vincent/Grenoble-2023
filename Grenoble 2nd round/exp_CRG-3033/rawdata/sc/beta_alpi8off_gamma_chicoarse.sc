START	standard
CALL	setDCpi2.sc
I	AO_H	Amps	-0.01
I	AO_A1Z	Amps	-0.026
CALL	ifg_20s_3p.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-1.50000
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-1.00000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.87000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.74000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.61000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.48000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
I	AO_H	Amps	-0.01
I	AO_A1Z	Amps	-0.026
CALL	ifg_20s_3p.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-1.00000
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.35000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.22000
CALL	beta_alpi8off_90s.sc
CALL	gamma_90s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	-0.09000
CALL	beta_alpi8off_120s.sc
CALL	gamma_120s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.04000
CALL	beta_alpi8off_120s.sc
CALL	gamma_120s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.17000
CALL	beta_alpi8off_120s.sc
CALL	gamma_120s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.30000
CALL	beta_alpi8off_90s.sc
CALL	gamma_90s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.43000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.56000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.69000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.82000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	0.95000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	1.08000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	1.21000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	1.34000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	1.47000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
MS	S18.axis2.phaseshifter 1	1	mm	absolute	1.60000
CALL	beta_alpi8off_60s.sc
CALL	gamma_60s.sc
I	AO_H	Amps	-0.01
I	AO_A1Z	Amps	-0.026
CALL	ifg_20s_3p.sc
CALL	setA1Zoff.sc
END
