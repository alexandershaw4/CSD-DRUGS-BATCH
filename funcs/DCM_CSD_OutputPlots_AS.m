% Diagnostic output plots from a DCM for complex valued data
% AS

% Power spectra (Act & Pred)
%---------------------------------------
figure,PCSD(DCM.Hc{:}+DCM.Rc{:},DCM.Hc{:},'real',DCM.Hz,'Data & Predicted Spectra (real)');
figure,PCSD(DCM.Hc{:}+DCM.Rc{:},DCM.Hc{:},'imag',DCM.Hz,'Data & Predicted Spectra (imaginary)');


% Coherence & phase-delay (Act & Pred)
%---------------------------------------
[coh fsd]   = spm_csd2coh(DCM.Hc{:}+DCM.Rc{:},DCM.M.Hz);

figure,PCSD(coh,DCM.coh{:},[],DCM.Hz,'Coherence (Data & Prediction)');
figure,PCSD(fsd,DCM.fsd{:},[],DCM.Hz,'Phase-delay (Data & Prediction)');


% Directed transfer functions (Pred)
%---------------------------------------
figure, PCSD(DCM.dtf{:},[],'real',DCM.Hz,'Directed Transfer Function (real)');
figure, PCSD(DCM.dtf{:},[],'imag',DCM.Hz,'Directed Transfer Function (imag)');