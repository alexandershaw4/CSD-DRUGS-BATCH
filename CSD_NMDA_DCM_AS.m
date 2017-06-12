function CSD_NMDA_DCM_AS(m,StNo,G)
% Top level function for CSD DCM. Calls a  number of updated SPM functions,
% stored in ----. These include hacked versions of spm_dcm_csd and spm_dcm_csd_data so add them
% to your paths. The hacks enable things like spectral baselining,
% filtering and sliding-window smoothing. Doing these on-the-fly gives the
% ability to try out different values without having to go back to SPM or
% save changes to the dataset.
%
% This script is accompannied by 2 further functions:
% 1) MODELS.m, which contains all the model architectures you want to run
% 2) GroupDataLocs.m, which contains the locations of the data you want to
% run. This can handle multi-group data.
%
% Also note that you will want to scroll down and look at:
% 1) CustomPriors - here are some custom priors [currently for CMC model]
% 2) PrepData     - IMPORTANT: options for the empirical data
% 3) checktrialcodes - special routine for checking trial codes
%
% USAGE:
% CSD_NMDA_DCM_AS(m, StNo,G) where: 
% m    = the models you want to run (int or vector) defined in MODELS.m
% StNo = subjects
% G    =  Group (corresponding to GroupDataLocs.m)
% 
% e.g. CSD_DCM_AS(1:21,1:16,1); will run models 1:21, subs 1:16 in
% group1.
%
% AS2016 [DCM]


% ensure SPM12 + custom functions


Mst       = m(end);         % Stop at model m(end)
[M,C,L,B] = MODELS;         % Model spaces
[s,p]     = GroupDataLocs;  % Data locations

s = s{G};
p = p(G);


% IGNORE... [calls subs functions]
lM         = length(M);
try if Mst ~=0;
       lM  = Mst;
    end;
end
for i = m(1):lM
    
    fprintf('\nrunning model %d\n',i);
    try kount = kount + 1; catch kount = 1; end
    if  kount > 1; StNo(1) = 1; end
    Do(M{i},C{i},L{i},B{i},i,s,p,StNo,p.d,p.f);
    
end

end

function DCM = PrepData(DCM,Ns,tCode)
% Sets options for the empirical data, with new options including
% baselining, filtering and using a set % of trials.

DCM.M.U            = sparse(diag(ones(Ns,1)));  ... ignore [modes]
DCM.options.trials = tCode;                     ... trial code [GroupDataLocs]
DCM.options.Tdcm   = [1 2000];       ... peristimulus time
DCM.options.Fdcm   = [13 30];        ... frequency window
%DCM.options.D      = 2;             ... downsample
DCM.options.han    = 1;             ... apply hanning window
DCM.options.h      = 4;             ... can't remember
DCM.options.DoData = 1;             ... leave on [custom]
%DCM.options.baseTdcm = [-100 0];      ... baseline times [new!]
%DCM.options.Fltdcm = [1 45];        ... bp filter [new!]
%DCM.options.Window = 8;             ... sliding window [new!]

DCM.options.analysis      = 'CSD';  ... analyse type
DCM.xY.modality           = 'LFP';  ... ECD or LFP data? [LFP]
DCM.options.spatial       = 'LFP';  ... spatial model [LFP]
DCM.options.model         = 'cmm_nmda';  ... neural model
DCM.options.Nmodes        = length(DCM.M.U);

DCM.xY.name     = {'Angular_L','Angular_R','L_Paracentral_Lob','R_Paracentral_Lob'};
DCM.Sname       = DCM.xY.name';
DCM.xY.Ic       = [1:Ns];


DCM.options.UseButterband = [13 30];
DCM.options.UseWelch = 0;
DCM.options.Smooth = 0;

%DCM = spm_dcm_csd_data(DCM);
DCM = spm_dcm_csd_data_as(DCM);
DCM.options.DONE = 1;         ... don't force dcm_csd to recopmute this


end

function Do(M,C,L,B,nm,s,p,StNo,d,f)
% Make the DCM structure. 
% This will contain empirical data in DCM.xY, experiment [design] info in 
% DCM.xU, function handles to the generative [DCM.M.f], forward [DCM.M.G/g], 
% integrator [DCM.M.IS] and feature selection [DCM.M.FS] as well as the
% prior estimate and variances for those functions in DCM.M.pE/gE/pC/gC.


% IGNORE ME
h     = pwd; s = s; cd(s); % [for save]
DDir  = dir(d); DDir = {DDir.name};
for n = 1:length(DDir); try fi{n} = [s ls(deblank([DDir{n} '/' f]))]; end;end; N = n; cd(h);
Data  = fi;


% Start from subject n, or 1:
if isempty(StNo) || exist('StNo') == 0
     ns = 1;
else ns = StNo;
     fprintf('starting from subject %d',ns(1));
end


    % loop over subjects
    %----------------------------------------------------------------------
    for s = ns(1):ns(end)
    fprintf('\nrunning subject %d of %d\n',s,length(Data));
    
    
        % data naming & desing matrix
        %--------------------------------------------
        DCM          = [];
        [fp fn fe]   = fileparts(Data{s});
        DCM.name     = genvarname(['CSD_Mod_',num2str(nm),'_', fn, fe]);
        
        if exist([DCM.name '.mat']) == 2; 
            fprintf('found = skipping sub %d',s);
        else
            
        DCM.xY.Dfile = Data{s};
        Ns           = length(M);
        DCM.xU.X     = p.xU.X;       ... design matrix
        DCM.xU.name  = p.xU.name;    ... condition names
        tCode        = p.tCode;      ... condition index (in SPM)
        
        
        % Extrinsic connectivity - model spaces
        %---------------------------------------------
        ALL      = M;                          ... full connectivity
        DCM.B    = {B};                        ... trial specific
        DCM.A{3} = L;                          ... lateral [modulatory]
        DCM.C    = C;                          ... [exogenous] inputs
        
        DCM.A{1} = triu(ALL);                  ... forward
        DCM.A{2} = tril(ALL);                  ... backward
        
        DCM.B(2:length(DCM.xU.X))=DCM.B;
        
        % Functions
        DCM = CustomPriors(DCM,Ns);    % anything non-default / built in
        DCM = PrepData(DCM,Ns,tCode);  % evaluate empirical data
        DCM = MakeAndInvert(DCM);      % invert the model & save it
        end
    end

end

function DCM = CustomPriors(DCM,Ns)

% This will put a bunch of custom (non-DCM/SPM) functions or priors into a
% new structure in the DCM. These will be copied over in the call to
% spm_dcm_csd. As such this will require a custom version of spm_dcm_csd
% which includes a line to detect whether this struct has been provided and
% to copy over it's contents prior to inversion.

DCM.CUSTOM      = [];
DCM.CUSTOM.f    = 'spm_fx_cmm_NMDA';            ... generative model  (.m)
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