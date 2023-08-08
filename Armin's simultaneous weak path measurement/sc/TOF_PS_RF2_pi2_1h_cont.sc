START	standard
;R	RF GF 2      
CALL	RobotZ.sc
CALL	nullsetzer.sc

CALL	setRF2_pi_contin.sc
CALL	rocking.sc
CALL	setRF2_pi_half.sc
CALL	setRF2_amp_DC_off.sc
CALL	movePS_200.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=3600	period=304µ	binwidth=1µ	delay=95µ
E

CALL	setRF2_pi_contin.sc
CALL	rocking.sc
CALL	setRF2_pi_half.sc
CALL	setRF2_amp_DC_off.sc
CALL	movePS_160.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=3600	period=304µ	binwidth=1µ	delay=95µ
E
END
