function getdcmcsd(DCM)

csd = zeros(length(DCM.Hz),length(DCM.Sname),length(DCM.Sname));

y = spm_unvec( (spm_vec(DCM.Hc) + spm_vec(DCM.Rc)) , csd);
x = spm_unvec(  spm_vec(DCM.Hc) , csd );

