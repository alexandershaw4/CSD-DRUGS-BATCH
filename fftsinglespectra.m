function [f,mx,varargout]=fftsinglespectra(t,x,MoveAverageKernel,Fs)
    
    if(nargin<4),
        Fs=1/(t(2)-t(1));
    end
    Nyquist=Fs/2;
    FRange=0:0.1:Nyquist-0.1;
    
    if MoveAverageKernel==1000,
        %x=detrend(x);
        [mx,f]=pmtm(x,5,[],Fs);
         f=f';
        return;
    end
        
   if MoveAverageKernel==3000,
        %x=detrend(x);
        [mx,f]=pwelch(x,[],[],[],Fs);
         f=f';
        return;
   end
   
   if MoveAverageKernel==4000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
         f=f';
        return;
   end
   if MoveAverageKernel==14000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=fastsmooth(mx,4,3,1);
         f=f';
        return;
   end
   if MoveAverageKernel==15000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=fastsmooth(mx,8,3,1);
         f=f';
        return;
   end
   if MoveAverageKernel==16000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=savitzkyGolayFilt(mx,4,0,9);
         f=f';
        return;
   end
   if MoveAverageKernel==17000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=savitzkyGolayFilt(mx,2,0,7);
         f=f';
        return;
   end
    if MoveAverageKernel==18000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=savitzkyGolayFilt(mx,1,0,5);
         f=f';
        return;
    end
    if MoveAverageKernel==19000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=savitzkyGolayFilt(mx,4,0,21);
         f=f';
        return;
    end
    if MoveAverageKernel==20000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        deltaf=f(2)-f(1); %fstep size
        
        sd=(2.0/deltaf); %2HZ is best so far!
        
        %width=7; %4Hz wide window
        width=12*sd; % Make the window width to be +/-6 standard deviations
        if width<2,
            width=2;
        end
        
        mx=ksmoothgauss(mx,width,sd);
        f=f';
        return;
   end
   if MoveAverageKernel==9000,
        %x=detrend(x);
        [mx,f]=periodogram(x,[],[],Fs);
        mx=medfilt1(mx,7);
         f=f';
        return;
   end
    if MoveAverageKernel==6000,
        %x=detrend(x);
        Order=64; %Estimate from AR model analysis one  dataset
        [mx,f]=pyulear(x,Order,FRange,Fs);
         f=f';
        return;
   end
     if MoveAverageKernel==8000,
        %x=detrend(x);
        Order=64; %Estimate from AR model analysis one  dataset
        [mx,f]=pmcov(x,Order,FRange,Fs);
         f=f';
        return;
   end
   
   if MoveAverageKernel==5000,
        WindowSize=max(t)-min(t);
        MinF=2/WindowSize;
        % Make this an integer frequency
        MinF=ceil(MinF);
        DT=3;
        Nyquist=Fs/(2*DT);
        f=MinF:0.1:Nyquist;
        vedata1=x;
        vedata1=detrend(x);
        vedata1=vedata1(1:DT:length(vedata1));
        mar=NEWbarf_spm_ar(vedata1,16);
        mx=abs(NEWbarf_spm_ar_freq(mar,f,Fs/DT))';
         %mx=sqrt(mx); % Want amplitude rather than power
        return;
   end
   
   if MoveAverageKernel==7000,
        WindowSize=max(t)-min(t);
        MinF=2/WindowSize;
        % Make this an integer frequency
        MinF=ceil(MinF);
        DT=1;Order=50;
        Nyquist=Fs/(2*DT);
        f=MinF:0.1:Nyquist;
        vedata1=x;
        vedata1=detrend(x);
        vedata1=vedata1(1:DT:length(vedata1));
        mar=NEWbarf_spm_ar(vedata1,Order);
        mx=abs(NEWbarf_spm_ar_freq(mar,f,Fs/DT))';
         %mx=sqrt(mx); % Want amplitude rather than power
        return;
   end
        
%     %x=detrend(x);
%     
%     %[mx,f]=pmtm(detrend(x),[],[],Fs);
%     % Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
%     nfft= 2^(nextpow2(length(x)));
% 
%     % Take fft, padding with zeros so that length(fftx) is equal to nfft 
%     fftx = fft(x,nfft); 
% 
%     % Calculate the numberof unique points
%     NumUniquePts = ceil((nfft+1)/2); 
% 
%     % FFT is symmetric, throw away second half 
%     fftx = fftx(1:NumUniquePts); 
% 
%     % Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
%     mx = abs(fftx)/length(x); %need to change from length(x) - need to divide by 2 then this effectively scales by Fs and turns it into true spectral density to match the periodogram... 
% 
%     % Take the square of the magnitude of fft of x. 
%     % KRISH MOD - I WANT THE AMPLITUDE SPECTRUM, RATHER THAN THE POWER, so do not square!
%     %mx = mx.^2; 
%     
%     % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
%     % The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.
%     if rem(nfft, 2) % odd nfft excludes Nyquist point
%         mx(2:end) = mx(2:end)*2;
%     else
%         mx(2:end -1) = mx(2:end -1)*2;
%     end
%     % Turn into power spectral density
%     
%     mx=mx.*conj(mx);
%     % Smooth the estimate by +/- MoveAverageKernel if KernelSize is 2000
%     % use an optimised spline function
%     if MoveAverageKernel==2000 
%         mx=smooth1q(mx);
%     end
%     if MoveAverageKernel~=2000,
%         mx=moving_average(mx,MoveAverageKernel);
%     end
%     
%     % This is an evenly spaced frequency vector with NumUniquePts points. 
%     f = (0:NumUniquePts-1)*Fs/nfft; 
    
    [mx,f]=periodogram(x,[],[],Fs);
    
    sudo_imag = imag(fft(x,length(mx)));
    varargout{1} = sudo_imag;
    
    f=f';    
    if MoveAverageKernel==2000 
        mx=smooth1q(mx);
    end
    if MoveAverageKernel~=2000,
        mx=moving_average(mx,MoveAverageKernel);
    end
    
  return

  function [ar] = NEWbarf_spm_ar (Z,p,verbose)
% Bayesian autoregressive modelling
% FORMAT [ar] = spm_ar (Z,p,verbose)
%
% y_pred (t) = -\sum_{i=1}^p a_i y (t-i) + e (t)
% Note the sign and ordering 
%
% The noise, e(t), is Gaussian
%
% Z             [N x 1] univariate time series 
% p             (scalar) order of model
% verbose       1=print out fitting progress (default=0)
%
% ar            data structure
% ----------------------------------
% ar.a_mean     AR coefficients
% ar.a_cov
% ar.mean_beta  error precision 
% ar.b_beta
% ar.c_beta
% ar.mean_alpha weight precision 
% ar.b_alpha
% ar.c_alpha
% ar.y          targets
% ar.y_pred     predictions
% ar.r2         proportion of variance explained
% ar.p          model order
% ar.fm         negative free energy
%
% For algorithmic details see:
%
% W.D. Penny and S.J. Roberts. Bayesian Methods for Autoregressive Models.
% In IEEE Workshop on Neural Networks for Signal Processing, Sydney Australia, 2000
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_ar.m 1276 2008-03-28 18:29:19Z guillaume $

if nargin < 2, 
   disp('spm_ar.m needs at least two arguments'); 
   return
end

if nargin < 3 | isempty(verbose)
    verbose=0;
end

Z=Z(:);
y=Z(p+1:end);
for i=1:p,
    x(:,i)=Z(p-i+1:end-i);
end

N=length(Z);
if abs(mean(Z)) > 3*(std(Z)/sqrt(N))
  % If mean is greater than 3SE away from 0
  % disp('Warning from vbar: mean subtracted from data');
  Z=Z-mean(Z);
end

% Initialise coefficients to maximum likelihood solution
% if p > 1
%     % In case columns of x are collinear
%     [ux,dx,vx]=svd(x);
%     ddx=diag(dx);
%     svd_tol=max(ddx)*eps*p;
%     rank_X=sum(ddx > svd_tol);
%     ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
%     ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
%     Xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';
%     X2=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)';
%     
%     a_mean= Xp*y;
%     y_pred= x*a_mean;
%     v=mean((y-y_pred).^2);
%     a_cov = v*X2;
%     xtx=X2;
% else
a_mean = pinv(x)*y;
y_pred = x*a_mean;
v=mean((y-y_pred).^2);
xtx=x'*x;
a_cov = v*inv(xtx);
    

% Setting to these values gives updates for mean_alpha, mean_beta
% approx to evidence framework 
b_alpha_prior=1000;
c_alpha_prior=0.001;
mean_alpha_prior=b_alpha_prior*c_alpha_prior;
b_beta_prior=1000;
c_beta_prior=0.001;

xty=x'*y;
xt=x';
yty=y'*y;

max_iters=32;
lik=[];
tol=0.0001;
for it=1:max_iters,

  E_w=a_mean'*a_mean;
  % Update weight precision
  b_alpha=0.5*E_w+0.5*trace(a_cov)+(1/b_alpha_prior);
  b_alpha=1/b_alpha;
  c_alpha=0.5*p+c_alpha_prior;
  mean_alpha=b_alpha*c_alpha;
  
  E_d_av=0.5*yty-a_mean'*xty;
  E_d_av=E_d_av+0.5*a_mean'*xtx*a_mean;
  E_d_av=E_d_av+0.5*trace(a_cov*xtx);

  % Update noise precision
  b_beta=1/(E_d_av+(1/b_beta_prior));
  c_beta=0.5*N+c_beta_prior;
  mean_beta=b_beta*c_beta;
  
  % Update weights
  a_cov=inv(mean_beta*xtx+mean_alpha*eye(p));
  a_mean=mean_beta*a_cov*xty;
  
  % Calculate f_m (negative free energy)
  l_av=0.5*N*(psi(c_beta)+log(b_beta))-0.5*N;
  kl_weights=spm_kl_normal(a_mean,a_cov,zeros(1,p),(1/mean_alpha)*eye(p));
  kl_alpha=spm_kl_gamma(b_alpha,c_alpha,b_alpha_prior,c_alpha_prior);
  kl_beta=spm_kl_gamma(b_beta,c_beta,b_beta_prior,c_beta_prior);
  f_m=l_av-kl_weights-kl_alpha-kl_beta;
  if verbose
      disp(sprintf('Iteration %d L_av=%1.3f KL_w=%1.3f KL_alpha=%1.3f KL_beta=%1.3f Fm=%1.3f',it,l_av,kl_weights,kl_alpha,kl_beta,f_m));
  end
  % Convergence criterion
  oldlik=lik;
  lik=f_m;
  
  if (it<=2)
    likbase=lik;
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase))
    break;
  end;

end

% Now reverse sign of coefficients to keep format
% with ar_spec etc. (covariances will be unchanged)
ar.a_mean=-a_mean;

% Load up data structure
ar.y_pred = x*a_mean;
ar.y=y;
vy=std(y)^2;
ey=std(y-ar.y_pred)^2;
ar.r2=(vy-ey)/vy; 

ar.p=p;
ar.a_cov=a_cov;
ar.mean_beta=mean_beta;
ar.mean_alpha=mean_alpha;
ar.b_beta=b_beta;
ar.c_beta=c_beta;
ar.b_alpha=b_alpha;
ar.c_alpha=c_alpha;
ar.fm=f_m;
ar.l_av=l_av;
ar.iterations=it;
return
  
  function [p] = NEWbarf_spm_ar_freq (ar, freq, fs)
% Compute spectra from AR coefficients
% FORMAT [p] = spm_ar_freq (ar, freq, fs)
%
% ar    AR model data structure (see spm_ar.m)
% freq  [Nf x 1] vector containing list of frequencies
% fs    sample rate
%
% p     [Nf x 1] vector containing power estimates
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_ar_freq.m 1143 2008-02-07 19:33:33Z spm $

Nf=length(freq);
Np=length(ar.a_mean);
T=1/fs;
freq=freq(:);
    
if Np < Nf
    % Usual case - fewer AR coeffs than freqs to evaluate at
    exponent=-i*2*pi*freq*T;
    
    denom=ones(Nf,1);
    for j=1:Np,
        denom=denom+ar.a_mean(j)*exp(j*exponent);
    end
    denom=abs(denom).^2;
    p=1./denom;
else
    p_order=[1:Np]';
    for f=1:Nf,
        denom=abs(1+sum(ar.a_mean.*exp(-i*2*pi*freq(f)*T*p_order)));
        p(f)=1/abs(denom)^2;
    end
end

noise_dev=sqrt(1/ar.mean_beta);
p=p*noise_dev*T;
return
