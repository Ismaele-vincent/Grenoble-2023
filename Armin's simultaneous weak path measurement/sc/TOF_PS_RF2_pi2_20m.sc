START	standard
;R	RF GF 2      
CALL	RobotZ.sc
CALL	nullsetzer.sc
CALL	setRF2_pi_half.sc

CALL	setRF2_amp_DC_off.sc
CALL	movePS_240.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
CALL	setRF2_amp_DC_off.sc
CALL	movePS_280.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
CALL	setRF2_amp_DC_off.sc
CALL	movePS_320.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E


CALL	setRF2_amp_DC_off.sc
CALL	movePS_0.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
CALL	setRF2_amp_DC_off.sc
CALL	movePS_40.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
CALL	setRF2_amp_DC_off.sc
CALL	movePS_80.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
CALL	setRF2_amp_DC_off.sc
CALL	movePS_120.sc
CALL	setRF2_amp_DC_on.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_pi2\*.tof
TOF	4	time=1200	period=304µ	binwidth=1µ	delay=95µ
E
END
