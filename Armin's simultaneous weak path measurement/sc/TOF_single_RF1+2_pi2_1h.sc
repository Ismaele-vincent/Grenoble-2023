START	standard
;R	RF 1+2  Pi half    
CALL	RobotZ.sc
CALL	nullsetzer.sc
CALL	setDC_pi_half.sc

CALL	setRF1_pi_contin.sc
CALL	rocking.sc
CALL	setRF1_pi_half_contin_bothGF.sc
CALL	setRF2_pi_half_contin_bothGF.sc
CALL	movePS_0.sc
CALL	setRF1_pi_half_bothGF.sc
CALL	setRF2_pi_half_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF1_pi2\*.tof
TOF	4	time=3600	period=304µ	binwidth=1µ	delay=95µ
E


END
