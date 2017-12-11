function [data,model,hz] = DCM_CSD_Plots(P)
% Diagnostic output plots from a GROUP of DCMs for complex valued data
% AS

for i = 1:length(P)
    load(P{i});
    
    Hz     = DCM.Hz;
    Hc{i}  = DCM.Hc{:};
    Rc{i}  = DCM.Rc{:};
    coh{i} = DCM.coh{:};
    fsd{i} = DCM.fsd{:};
    dtf{i} = DCM.dtf{:};
    
end

clear DCM

DCM.Hz  = Hz;
DCM.Hc  = {squeeze(( cat(4,Hc{:})  ))};
DCM.Rc  = {squeeze(( cat(4,Rc{:})  ))};
DCM.coh = {squeeze(( cat(4,coh{:}) ))};
DCM.fsd = {squeeze(( cat(4,fsd{:}) ))};
DCM.dtf = {squeeze(( cat(4,dtf{:}) ))};

if nargout;
    data  = DCM.Hc{:}+DCM.Rc{:};
    model = DCM.Hc{:};
    hz    = DCM.Hz;
    return;
end

makeplots(DCM);

end

function makeplots(DCM)


% Power spectra (Act & Pred)
%---------------------------------------
figure,PCSD(DCM.Hc{:}+DCM.Rc{:},DCM.Hc{:},'real',DCM.Hz,'Data & Predicted Spectra (real)');
figure,PCSD(DCM.Hc{:}+DCM.Rc{:},DCM.Hc{:},'imag',DCM.Hz,'Data & Predicted Spectra (imag)');


% Coherence & phase-delay (Act & Pred)
%---------------------------------------
[coh fsd]   = spm_csd2coh(DCM.Hc{:}+DCM.Rc{:},DCM.Hz);

figure,PCSD(coh,DCM.coh{:},[],DCM.Hz,'Coherence (Data & Prediction)');
figure,PCSD(fsd,DCM.fsd{:},[],DCM.Hz,'Phase-delay (Data & Prediction)');


% Directed transfer functions (Pred)
%---------------------------------------
figure, PCSD(DCM.dtf{:},[],'real',DCM.Hz,'Directed Transfer Function (real)');
figure, PCSD(DCM.dtf{:},[],'imag',DCM.Hz,'Directed Transfer Function (imag)');

end