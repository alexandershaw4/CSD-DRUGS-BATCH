function DCM = spm_dcm_csd_data_as(DCM)
% gets cross-spectral density data-features using a VAR model
% FORMAT DCM = spm_dcm_csd_data_as(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile        - name of data file
%    DCM.M.U             - channel subspace
%    DCM.options.trials  - trial to evaluate
%    DCM.options.Tdcm    - time limits
%    DCM.options.Fdcm    - frequency limits
%    DCM.options.D       - Down-sampling
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] sampled
%    DCM.xY.dt      - sampling in seconds [s] (down-sampled)
%    DCM.xY.U       - channel subspace
%    DCM.xY.y       - cross spectral density over sources
%    DCM.xY.csd     - cross spectral density over sources
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins
%    DCM.xY.code    - trial codes evaluated
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ssr_data.m 4096 2010-10-22 19:40:34Z karl $
 
% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% ensure spatial modes have been computed (see spm_dcm_ssr)
%-------------------------------------------------------------------------
try
    DCM.M.U;
catch
    errordlg('Please estimate this model first');
    error('')
end

% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch
    try
        [p,f]        = fileparts(Dfile);
        D            = spm_eeg_load(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch
        try
            [f,p]        = uigetfile('*.mat','please select data file');
            name         = fullfile(p,f);
            D            = spm_eeg_load(name);
            DCM.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
            return
        end
    end
end

 
% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    [mod, list] = modality(D, 0, 1);

    if isequal(mod, 'Multimodal')
        qstr = 'Only one modality can be modelled at a time. Please select.';
        if numel(list) < 4
            % Nice looking dialog. Will usually be OK
            options = [];
            options.Default = list{1};
            options.Interpreter = 'none';
            DCM.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
        else
            % Ugly but can accomodate more buttons
            ind = menu(qstr, list);
            DCM.xY.modality = list{ind};
        end
    else
        DCM.xY.modality = mod;
    end
end



if ~isfield(DCM.xY, 'Ic')
    DCM.xY.Ic  = setdiff(D.meegchannels(DCM.xY.modality), D.badchannels);
end

Ic        = DCM.xY.Ic;
Nc        = length(Ic);
Nm        = size(DCM.M.U,2);
DCM.xY.Ic = Ic;

% options
%--------------------------------------------------------------------------
try
    DT    = DCM.options.D;
catch
    DT    = 1;
end
try
    trial = DCM.options.trials;
catch
    trial = D.nconditions;
end

% KRISH HAS COMMENTED OUT THIS CODE AS I WANT TO SPECIFY MY OWN DOWNSAMPLE RATE - THIS ACTUALLY HAS A BIG EFFECT ON THE SPECTRAL SHAPE
%
% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
%if DT/D.fsample < 0.004
%    DT            = ceil(0.004*D.fsample);
%    DCM.options.D = DT;
%end


% get peristimulus times
%--------------------------------------------------------------------------
%try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    DCM.xY.Time = 1000*D.time; % ms
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    [i, T1]     = min(abs(DCM.xY.Time - T1));
    [i, T2]     = min(abs(DCM.xY.Time - T2));
    
    
    %KRISH HACK - SEE IF BASELINE HAS BEEN SPECIFIED
    DOBASE=0;
    try
           baseT1          = DCM.options.baseTdcm(1);
           baseT2          = DCM.options.baseTdcm(2);
           [i, baseT1]     = min(abs(DCM.xY.Time - baseT1));
           [i, baseT2]     = min(abs(DCM.xY.Time - baseT2));
           baseIt=[baseT1:DT:baseT2];
           DCM.xY.baseIt=baseIt;
           DOBASE=1;
    end
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';               % indices - bins
    DCM.xY.pst  = DCM.xY.Time(It);           % PST
    DCM.xY.It   = It;                        % Indices of time bins
    DCM.xY.dt   = DT/D.fsample;              % sampling in seconds
    Nb          = length(It);                % number of bins
    fprintf(1,'Time Indices=%d:%d in steps of %d\n',T1,T2,DT);
% catch
%     errordlg('Please specify time window');
%     error('')
% end
 
% get frequency range
%--------------------------------------------------------------------------
%try
    Hz1     = DCM.options.Fdcm(1);          % lower frequency
    Hz2     = DCM.options.Fdcm(2);          % upper frequency
    
    % Krish Hack to allow finer frequency resolution (mostly NOT for DCM
    % but for using this function to make CSDs...
    Fstep=1;
    try
        Fstep=DCM.options.FrequencyStep;
        fprintf(1,'Using fstep=%f\n',Fstep);
    catch
        Fstep=1;
    end
% catch
%     pst     = DCM.xY.pst(end) - DCM.xY.pst(1);
%     Hz1     = max(ceil(2*1000/pst),4);
%     if Hz1 < 8;
%         Hz2 = 48;
%     else
%         Hz2 = 128;
%     end
% end

 
% Frequencies
%--------------------------------------------------------------------------
DCM.xY.Hz  = fix(Hz1:Fstep:Hz2);             % Frequencies
DCM.xY.Hz  = (Hz1:Fstep:Hz2);             % Frequencies: KRISH HACK TO ALLOW HIGHER RESOLUTION
Nf         = length(DCM.xY.Hz);        % number of frequencies
Ne         = length(trial);            % number of ERPs
 
% get induced responses (use previous CSD results if possible). COMMENTED OUT BY KRISH! 
%==========================================================================
% try
%     if size(DCM.xY.csd,2) == Ne;
%         if size(DCM.xY.csd{1},1) == Nf;
%             if size(DCM.xY.csd{1},2) == Nm;
%                 DCM.xY.y  = spm_cond_units(DCM.xY.csd,'csd');
%                 return
%             end
%         end
%     end
% end
 
% Cross spectral density for each trial type
%==========================================================================



condlabels = D.condlist;

try
    DCM.xY.y=DCM.Kinitial.ForceSpectra;
    fprintf('Using user-specified spectra...\n');
    DCM.xY.U    = DCM.M.U;
    DCM.xY.code = condlabels(trial);
    return;
end

for i = 1:Ne;
   clear Pfull
    
    % trial indices
    %----------------------------------------------------------------------
    %c = D.pickconditions(condlabels{trial(i)});
    c = D.indtrial(condlabels(trial(i)), 'GOOD');
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    try c = c(1:512); end
    Nt    = length(c);

    
    % Get data
    %----------------------------------------------------------------------
    P     = zeros(Nf,Nm,Nm);
    % KRISH HACK
    
    if DOBASE==1,
       Pbase=zeros(Nf,Nm,Nm);
    end
    
    Fs = 1000/(DCM.xY.Time(2)-DCM.xY.Time(1));
  
     fprintf('\n');
%      % KRISH MOD Catch if there is mistakenly more than one channel in the  data...take the first channel...
%         if length(Ic)>1,
%             fprintf(1,'Found multichannel data, setting to first channel');
%             Ic=1;
%         end
        
    for j = 1:Nt
        
        if (j==1)
             fprintf('%s --> evaluating condition %i (trial %i) scale=%f', DCM.xY.Dfile,i,j,full(DCM.M.U));
        end
        
        Y   = full(double(D(Ic,It,c(j))'));%*DCM.M.U));
        
        
        %KRISH HACK
            if DOBASE==1,
                Ybase=full(double(D(Ic,baseIt,c(j))'*DCM.M.U));
            end
            
        % KRISH MOD Catch if there is mistakenly more than one channel in the  data...take the first channel...
%         if size(Y,2)>1,
%                 if (j==1), 
%                     fprintf(1,'  taking first channel from multichannel data');
%                 end
%                 Y=squeeze(Y(:,1));
%         end
        try
            DoButterband=DCM.options.UseButterband;
            %Y=butterband(Y,DoButterband(1),DoButterband(2),Fs);
            Y = bandpassfilter(Y',Fs,DoButterband)';
            %Ybase=butterband(Ybase,DoButterband(1),DoButterband(2),Fs);
            if (j==1), 
                fprintf(1,'   filtering from %f-%fHz Fs=%f',DoButterband(1),DoButterband(2),Fs);
            end
        end
        UseWelch=0;
        try % krish Test to use Welch method for spectral estimation...
            UseWelch=DCM.options.UseWelch;
        catch
            UseWelch=0;
        end
        if(UseWelch~=0),
        if (j==1),
            if(UseWelch==6000),
                fprintf(1,'....Using pyulear'); 
            end
            if(UseWelch==8000),
                fprintf(1,'....Using PMCOV64 '); 
            end
           if(UseWelch==5000),
                fprintf(1,'....Using NewAR 3 16'); 
            end
           if(UseWelch==7000),
                fprintf(1,'....Using NewAR 1 50'); 
            end
             if(UseWelch==20000),
                fprintf(1,'....Using Gauss smooth Peridogram'); 
            end
             if(UseWelch==4000),
                fprintf(1,'....Using Peridogram'); 
            end
           if(UseWelch==3000),
                fprintf(1,'....Using Welch'); 
            end
            if(UseWelch==1000),
                fprintf(1,'....Using PMTM'); 
            end
           if(UseWelch==1999),
                fprintf(1,'....Using RobustSpectra'); 
            end
           if(UseWelch<1000),
                fprintf(1,'....Using smoothed periodogram, kernel=+/- %f Hz',UseWelch); 
           end
        end 
           
           w=ones(length(Y),1);
           winpower=1;
        
            if (DCM.options.han ~= 0),
                w=tukeywin(length(Y),1);
                winpower=sum(w)/length(w);
                if (j==1), 
                    fprintf(' with tukeywin');
                end
            end

            %Ymod=Y-mean(Y);
            Ymod = Y - repmat(mean(Y),[size(Y,1),1]);
            
            %Ymod=Ymod.*w;
            Ymod = Ymod.*repmat(w,[1,size(Ymod,2)]);
            
            if UseWelch~=1999 && UseWelch ~= 1e6 && UseWelch ~=1e7
                [F,Pfull(j,:)]=fftsinglespectra(DCM.xY.Time/1000,Ymod,UseWelch);
            end
            
            if UseWelch==1999,
                for nch = 1:size(Ymod,2)
                    [F,Pfull(j,nch,:),comp(j,nch,:)]=fftsinglespectra(DCM.xY.Time(It)/1000,Ymod(:,nch)',2);
                end
            end
            
            if UseWelch==1e6
                Pfull(j,:,:) = fft(HighResMeanFilt(Ymod,1,4));
                
                F = 1/DCM.xY.dt*(0:(size(Pfull,2)/2))/size(Pfull,2);
            end
            
            if UseWelch == 1e7
                for xx = 1:Nc
                    for yy = 1:Nc
                        [Pfull(j,:,xx,yy),F]=cpsd(Ymod(:,xx),Ymod(:,yy),[],[],[],1/DCM.xY.dt);
                    end
                end
                
            end
            
            try dF = F;end
               Pfull(j,:)=Pfull(j,:)/winpower;
            
            if DOBASE==1 && UseWelch ~=1e6                
                if (j==1),
                    fprintf(' with baselining');
                end
                
                %Ymod=Ybase-mean(Ybase);
                
                Ymod = Ybase-repmat(mean(Ybase),[size(Ybase,1),1]);
                
                %if SpectCondition==2,
                %    Ymod=Ybase;
                %end

                if size(Ymod,1) ~= size(Y,1)
                    Ymod = mean(Ymod);
                    Ymod = repmat(Ymod,[size(Y,1),1]);
                end
                    
                
                Ymod=Ymod.*repmat(w,[1,size(Y,2)]);
                for nch = 1:size(Ymod,2)
                    [F,Pfullbase(j,nch,:)]=fftsinglespectra(DCM.xY.Time(It)/1000,Ymod(:,nch)',UseWelch);
                end
                Pfullbase(j,:)=Pfullbase(j,:)/winpower;
            end
            
        end
        if (UseWelch==0),
            %mar = spm_mar(Y,8); % ORIGINAL
            if (j==1), fprintf(1,'....Using AR spectra');end
            
            %Y = HighResMeanFilt(Y,1,2);
            %Y = Y - repmat(min(Y),[size(Y,1),1]); % make positive
            
            fprintf('Calculating CSD for trial %d of %d\n',j,Nt);
            mar = spm_mar(Y,3);%16 % CHANGE BY KRISH
            mar = spm_mar_spectra(mar,DCM.xY.Hz,1/DCM.xY.dt);
            P   = P + (mar.P);
            %Pfull = P;
            % KRISH HACK...
            if DOBASE==1,
                if (j==1), fprintf(' with baselining');end
                
                Ybase = HighResMeanFilt(Ybase,1,2);
                Ybase = Ybase - repmat(min(Ybase),[size(Ybase,1),1]); % 
                
                marbase = spm_mar(Ybase,8);%16 % CHANGE BY KRISH
                marbase = spm_mar_spectra(marbase,DCM.xY.Hz,1/DCM.xY.dt);
                Pbase   = Pbase + (marbase.P);
            end
        end
        if (j==1), 
            fprintf('\n');
        end
    end
    
    if UseWelch==0 && DCM.options.Smooth > 0
        % then apply smoothing to fourier space, too
        fprintf('Smoothing\n');
        for ii = 1:size(P,2)
            for jj = 1:size(P,3)
                P(:,ii,jj) = HighResMeanFilt(P(:,ii,jj),1,DCM.options.Smooth);
            end
        end
        
        
        if DOBASE
            for ii = 1:size(Pbase,2)
                for jj = 1:size(Pbase,3)
                    Pbase(:,ii,jj) = HighResMeanFilt(Pbase(:,ii,jj),1,DCM.options.Smooth);
                end
            end
            
        end
    end
    
    
    % store
    %----------------------------------------------------------------------
    % 
%         Pfull     = (mean(Pfull,1));
%     if DOBASE 
%         Pfullbase = (mean(Pfullbase,1)); 
%     end
%     
%     if UseWelch == 1e6
%         % plain fft
%         for ycsd = 1:size(Pfull,3)
%             for xcsd = 1:size(Pfull,3)
%                 y = squeeze(Pfull(:,:,ycsd));
%                 x = squeeze(Pfull(:,:,xcsd));
%                 CSD{i}(:,ycsd,xcsd) = y.*conj(x);
%                 
%                 %by = squeeze(Pfullbase(:,ycsd,:));
%                 %bx = squeeze(Pfullbase(:,xcsd,:));
%                 %bCSD{i}(:,ycsd,xcsd) = by.*conj(bx);
%                 
%             end
%         end
%         
%     elseif UseWelch == 1e7
%         % Pfull is trials by freqs by node by node
%         Pfull = squeeze(Pfull); % reduce over trials
%         CSD{i} = Pfull;
%         
%     elseif UseWelch ~= 0
%         % other methods
%         for ycsd = 1:size(Pfull,2)
%             for xcsd = 1:size(Pfull,2)
%                 y = squeeze(Pfull(:,ycsd,:));
%                 x = squeeze(Pfull(:,xcsd,:));
%                 CSD{i}(:,ycsd,xcsd) = y.*conj(x);
%                 
%                 %by = squeeze(Pfullbase(:,ycsd,:));
%                 %bx = squeeze(Pfullbase(:,xcsd,:));
%                 %bCSD{i}(:,ycsd,xcsd) = by.*conj(bx);
%                 
%             end
%         end
%     end
%     
%     %based{i} = CSD{i} - bCSD{i};
%     if UseWelch ~= 0
%         
%         if UseWelch == 1e6
%             newCSD{i} = CSD{i}(1:end/2+1,:,:);
%         else
%             newCSD{i} = CSD{i};
%         end
%         
%         % Note: not actually using difference spectra
%         for ycsd = 1:size(CSD{i},2)
%             for xcsd = 1:size(CSD{i},3)
%                 nCSD{i}(:,ycsd,xcsd) = spline(F,newCSD{i}(:,xcsd,ycsd),DCM.xY.Hz);
%             end
%         end
%         
%         DCM.xY.csd{i} = nCSD{i};
%         
%     else
%         % if mar model
%         P = (Pfull-Pbase);
%         DCM.xY.csd{i} = P;
%     end
%     
% end
%     DCM.xY.y      = spm_cond_units(DCM.xY.csd,'csd')   ;     
    
   % P = (P-Pbase);
   % DCM.xY.csd{i} = P/Nt;
    
    
    
    if(UseWelch==0),
        if DOBASE==1,
            P=(P-Pbase);
        end       
        DCM.xY.csd{i} = P/Nt;
        %return;
    end
    if (UseWelch~=0)
        
        if UseWelch~=1999,
            Pfull=mean(Pfull,1);
        end
        
        if UseWelch==1999,
                AlternativeFit=0;
                try
                    AlternativeFit=DCM.options.AlternativeFit;
                end
                
                if AlternativeFit==0,
                    [mnewspectra,LogFreqs,~,~,~,FitPar]=RobustSpectraFit(F,Pfull,2);
                end
                if AlternativeFit~=0,
                    [mnewspectra,LogFreqs,~,~,~,FitPar]=RobustSpectraFitTest(F,Pfull,2);
                    fprintf(1,'Using alternative fit');
                end
                F=exp(LogFreqs);
                mnewspectra=exp(mnewspectra);
                mnewspectra=mnewspectra-min(mnewspectra);
                 Pfull=mnewspectra;
                 DCM.RobustSpectra.Spectra=Pfull;
                 DCM.RobustSpectra.Freqs=F;
                 DCM.RobustSpectra.FitPars=FitPar;
        end
               
 
        P=spline(F,Pfull,DCM.xY.Hz);
        comp = spline(dF,comp,DCM.xY.Hz);

        if DOBASE==1,
            Pfullbase=mean(Pfullbase,1);
            P=P-spline(F,Pfullbase,DCM.xY.Hz);
        end
        
        P(find(isinf(P)))=0;
        P(find(isnan(P)))=0;
        P(find(P<=0))=0;
        P=P-min(P);
        DCM.xY.csd{i}=P';
        
        try DCM.xY.complex = comp; end
        
    end
    
end

 
% place cross-spectral density in xY.y
%==========================================================================
DCM.xY.y    = spm_cond_units(DCM.xY.csd,'csd'); 
% 
% %KRISH HACK! If No normalisation is selected, do not normalise by power
   DoNormalise=1;
    try
        DoNormalise=DCM.options.DoNormalise;
        fprintf(1,'\nDoNormalise=%f\n',DoNormalise);
    catch
        DoNormalise=1;
    end

 if DoNormalise==0,
     DCM.xY.y    = krish_cond_unitsNONORMALISE(DCM.xY.csd,'csd'); 
     fprintf(1,'Not normalising to power....\n');
 end
 
 if DoNormalise==2,
     DCM.xY.y    = DCM.xY.csd; 
     fprintf(1,'Simple Spectral Copy....\n');
 end
 
%  % make pos
%  Q = @squeeze;
%  for i = 1:length(DCM.xY.y)
%      y = DCM.xY.y{i};
%      for ii = 1:size(y,2)
%          for jj = 1:size(y,3)
%              y(:,ii,jj) = y(:,ii,jj) + min(y(:,ii,jj));
%          end
%      end
%      DCM.xY.y{i} = y;
%  end
 
 
 
DCM.xY.U    = DCM.M.U;
DCM.xY.code = condlabels(trial);


