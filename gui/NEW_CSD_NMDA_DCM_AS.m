function NEW_CSD_NMDA_DCM_AS(Config)

Files  = Config.Files;  % input files [array]
ModelF = Config.MODELS; % models function 
Design = Config.DESIGN;
addpath(Config.mPath);
addpath(Config.dPath);

[M,C,L,B] = feval(strrep(ModelF,'.m',''));
p         = feval(strrep(Design,'.m',''));

% loop over models & subjects
%-----------------------------
for j = 1:length(M)
    for i = 1:length(Files)
        Do(M{j},C{j},L{j},B{j},j,Files{i},Config,p);
    end
end


end

function Do(M,C,L,B,nm,File,Config,p)
% Make the DCM structure.
% This will contain empirical data in DCM.xY, experiment [design] info in
% DCM.xU, function handles to the generative [DCM.M.f], forward [DCM.M.G/g],
% integrator [DCM.M.IS] and feature selection [DCM.M.FS] as well as the
% prior estimate and variances for those functions in DCM.M.pE/gE/pC/gC.


% loop over subjects
%----------------------------------------------------------------------
fprintf('\nrunning subject %d: %s\n',nm,File);


% data naming & desing matrix
%--------------------------------------------
DCM          = [];
[fp fn fe]   = fileparts(File);
DCM.name     = (['CSD_Mod_',num2str(nm),'_', fn, fe]);

if exist([DCM.name '.mat']) == 2;
    fprintf('found = skipping sub %d',s);
else
    
    DCM.xY.Dfile = File;
    Ns           = length(M.F);
    DCM.xU.X     = p.xU.X;       ... design matrix
    DCM.xU.name  = p.xU.name;    ... condition names
    tCode        = p.tCode;      ... condition index (in SPM)
        
    
    % Extrinsic connectivity - model spaces
    %---------------------------------------------
    DCM.B    = {B};                        ... trial specific
    DCM.A{3} = L;                          ... lateral [modulatory]
    DCM.C    = C;                          ... [exogenous] inputs
        
    DCM.A{1} = M.F;                  ... forward
    DCM.A{2} = M.B;                  ... backward
        
    DCM.B(2:length(DCM.xU.X))=DCM.B;
    
    % Functions
    DCM = CustomPriors(DCM,Ns,Config);    % anything non-default / built in
    DCM = PrepData(DCM,Ns,tCode,Config);  % evaluate empirical data
    DCM = MakeAndInvert(DCM);      % invert the model & save it
end
end




function DCM = PrepData(DCM,Ns,tCode,Config)
% Sets options for the empirical data, with new options including
% baselining, filtering and using a set % of trials.

DCM.M.U            = sparse(diag(ones(Ns,1)));  ... ignore [modes]
DCM.options.trials = tCode;                     ... trial code [GroupDataLocs]

try   DCM.options.Tdcm = Config.pst;
catch DCM.options.Tdcm   = [1 2000];       ... peristimulus time
end

try   DCM.options.Fdcm = Config.FreqWin;
catch DCM.options.Fdcm = [13 30];          ... frequency window
end

try   DCM.options.D  = Config.Downsample;
catch DCM.options.D  = 0;                  ... downsample
end

try   DCM.options.han = Config.Han;
catch DCM.options.han = 1;                 ... apply hanning window
end 

DCM.options.h      = 4;             ... can't remember
DCM.options.DoData = 1;             ... leave on [custom]

try Config.BaseTimes;
    DCM.options.baseTdcm = Config.BaseTimes;
end

%DCM.options.baseTdcm = [-100 0];    ... baseline times [new!]
%DCM.options.Fltdcm = [1 45];        ... bp filter [new!]
%DCM.options.Window = 8;             ... sliding window [new!]

try   DCM.options.analysis = Config.Anal;
catch DCM.options.analysis = 'CSD';  ... analyse type
end

try   DCM.xY.modality = Config.Type;
catch DCM.xY.modality = 'LFP';       ... ECD or LFP data? [LFP]
end

try   DCM.options.spatial = Config.Spatial;
catch DCM.options.spatial = 'LFP';   ... spatial model [LFP]
end

try   DCM.options.model   = Config.Neural;
catch DCM.options.model   = 'cmm_nmda';  ... neural model
end


DCM.options.Nmodes  = length(DCM.M.U);
DCM.xY.name         = {'Angular_L','Angular_R','L_Paracentral_Lob','R_Paracentral_Lob'};
DCM.Sname           = DCM.xY.name';
DCM.xY.Ic           = [1:Ns];

try   DCM.options.UseButterband = Config.Filter;
catch DCM.options.UseButterband = [1 100];
end

if isempty(DCM.options.UseButterband);
    DCM.options = rmfield(DCM.options,'UseButterband');
end

DCM.options.UseWelch = 0;
DCM.options.Smooth   = 0;

%DCM = spm_dcm_csd_data(DCM);
DCM = spm_dcm_csd_data_as(DCM);
DCM.options.DONE = 1;         ... don't force dcm_csd to recopmute this


end



function DCM = CustomPriors(DCM,Ns,Config)

% This will put a bunch of custom (non-DCM/SPM) functions or priors into a
% new structure in the DCM. These will be copied over in the call to
% spm_dcm_csd. As such this will require a custom version of spm_dcm_csd
% which includes a line to detect whether this struct has been provided and
% to copy over it's contents prior to inversion.

DCM.CUSTOM      = [];
DCM.CUSTOM.pE   = [];                           ... intrinsic priors
%DCM.CUSTOM.pE.G = zeros(Ns,13);                 ... [local coupling]
%DCM.CUSTOM.pC.G = zeros(Ns,13);                 ... variance [off]

DCM.CUSTOM.pE.H = zeros(4,4,Ns);
DCM.CUSTOM.pC.H = DCM.CUSTOM.pE.H;

Self = find(diag(speye(Ns).*( DCM.A{1}+DCM.A{2} ))); % SP gain only if in model
DCM.CUSTOM.pC.H(2,2,Self) = 1/8;

DCM.CUSTOM.pE.T = zeros(Ns,3);                  ... population time const
DCM.CUSTOM.pC.T = zeros(Ns,3)+1/8;              ... variances AMP,GAB,NMDA

DCM.CUSTOM.pC.J = full(sparse(1:4,1,1,16,1))';
DCM.CUSTOM.pE.J = full(sparse(1:4,1,[.2 .8 0 .2],16,1))';

%DCM.CUSTOM.gE.J = zeros(Ns,4,4);
%DCM.CUSTOM.gE.J(:,:,1) = -70;
%
%
%DCM.CUSTOM.gE.J = sparse([1 3 7],1,[.2 .8 .2],8,1)'; ... contributing states
%DCM.CUSTOM.gC.J = sparse([1 3 7],1,1/8       ,8,1)'; ... variances
%DCM.CUSTOM.gC.J = repmat(DCM.CUSTOM.gC.J,[3 1]);     ... contrib sources
%DCM.CUSTOM.gE.J = repmat(DCM.CUSTOM.gE.J,[3 1]);     ... variance


end


function DCM = MakeAndInvert(DCM)

DCM.CUSTOM.nograph   = 1;
DCM = spm_dcm_csd(DCM);

end