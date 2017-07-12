function RefitCSD_AS(P,varargin)
% Refit a bunch of CSD DCMs, for example, post reduction / optimisation
%
% AS

for i = 1:length(P);
    
    % sort input class
    if     iscell(P);   load(P{i}); DCM.name = P{i};
    elseif isstruct(P);             DCM      = P;
    elseif ischar(P);   load(P);    DCM.name = P;    DoStop = 1;
    end
    
    
    % naming
    [fp fn fe] = fileparts(DCM.name);
    fn         = strrep(fn, 'opt','refit');
    try     fn = [varargin{1} '_' fn];   end
    DCM.name   = [fn fe];
    
    % standard [neural] priors & posteriors
    pE = DCM.M.pE;
    pC = DCM.M.pC;
    Ep = DCM.Ep;
    Cp = DCM.Cp;
    pF = DCM.F;
    
    % package
    M   = DCM.M;
    xU  = DCM.xU;
    xY  = DCM.xY;
    
    M.pE = Ep;
    
    %if ~isstruct(Cp); M.pC = spm_unvec(diag(Cp),Ep); else M.pC = Cp; end
    %if ~isstruct(Cg); M.gC = spm_unvec(diag(Cg),Eg); else M.gC = Cg; end
    
    % dimensions
    %--------------------------------------------------------------------------
    Nt      = length(xY.y);                  % number of trials
    Nr      = size(DCM.C,1);                 % number of sources
    Nu      = size(DCM.C,2);                 % number of exogenous inputs
    Ns      = size(xY.y{1},1);               % number of time bins
    Nc      = size(xY.y{1},2);               % number of channels
    Nx      = size(xU.X,2);                  % number of trial-specific effects
    
    
    
    [Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    
    
    % Data ID
    %--------------------------------------------------------------------------
    try
        try
            ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
        catch
            ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
        end
    catch
        ID = spm_data_id(DCM.xY.y);
    end
    
    pE = DCM.M.pE;
    pC = DCM.M.pC;
    
    % Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
    %==========================================================================
    warning('off','SPM:negativeVariance');
    dp  = spm_vec(Qp) - spm_vec(pE);
    Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
    warning('on', 'SPM:negativeVariance');
    
    
    % predictions (csd) and error (sensor space)
    %--------------------------------------------------------------------------
    Hc  = spm_csd_mtf(Qp,DCM.M,DCM.xU);                      % prediction
    Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
    
    
    % predictions (source space - cf, a LFP from virtual electrode)
    %--------------------------------------------------------------------------
    M             = rmfield(DCM.M,'U');
    M.dipfit.type = 'LFP';
    
    
    M.U         = 1;
    M.l         = length(DCM.Sname); % Ns
    qp          = Qp;
    qp.L        = ones(1,M.l);             % set virtual electrode gain to unity
    qp.b        = qp.b - 32;              % and suppress non-specific and
    qp.c        = qp.c - 32;              % specific channel noise
    
    [Hs Hz dtf] = spm_csd_mtf(qp,M,DCM.xU);
    [ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
    [coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
    DCM.dtf     = dtf;
    DCM.ccf     = ccf;
    DCM.coh     = coh;
    DCM.fsd     = fsd;
    DCM.pst     = pst;
    DCM.Hz      = Hz;
    
    
    % store estimates in DCM
    %--------------------------------------------------------------------------
    DCM.Ep = Qp;                   % conditional expectation
    DCM.Cp = Cp;                   % conditional covariance
    DCM.Pp = Pp;                   % conditional probability
    DCM.Hc = Hc;                   % conditional responses (y), channel space
    DCM.Rc = Ec;                   % conditional residuals (y), channel space
    DCM.Hs = Hs;                   % conditional responses (y), source space
    DCM.Ce = exp(-Eh);             % ReML error covariance
    DCM.F  = F;                    % Laplace log evidence
    DCM.ID = ID;                   % data ID
    
    DCM.options.Nmodes = length(DCM.Sname);
    
    
    % and save
    %--------------------------------------------------------------------------
    [fp fn fe] = fileparts(DCM.name);
    DCM.name = [fn fe];
    save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
    
    
    % do not loop for single subject DCMs:
    try if DoStop; return; end; end
    
    
end