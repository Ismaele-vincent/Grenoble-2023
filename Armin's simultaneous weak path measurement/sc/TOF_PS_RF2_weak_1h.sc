START	standard
;R	RF GF 2      
CALL	RobotZ.sc
CALL	nullsetzer.sc
CALL	setRF1_weak_bothGF.sc
CALL	setRF1_amp_off.sc
CALL	setDC_pi_half.sc

CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_200.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_240.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_280.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_320.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E


CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_0.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_40.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_80.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_120.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E
CALL	setRF2_pi_contin_bothGF.sc
CALL	rocking.sc
CALL	setRF2_weak_contin_bothGF.sc
CALL	movePS_160.sc
CALL	setRF2_weak_bothGF.sc
;********Start
F	D:\data\Cycle 190\exp_CRG-2887\rawdata\sc\TOF\RF2_weak\*.tof
TOF	4	time=3600	period=304�	binwidth=1�	delay=95�
E

END