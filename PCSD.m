function PCSD(X,Y,varargin)
% Plot A/CSD from 3D double input: x(nf,ns,ns)
%            
% Examples:
%           PCSD(X);          % plot 3D double  x [real part]
%           PCSD(X,Y);        % plot 3D doubles x & y 
%           PCSD(X,Y,'imag'); % imaginary part only
%           PCSD(X,Y,'real'); % real part only 
%           PCSD(X,[],[],Hz); % provide freq/time steps
% AS

try Y                                  ; catch Y  = X*0;         end
try s  = eval(['@' lower(varargin{1})]); catch s  = @real;       end
try Hz = varargin{2}                   ; catch Hz = 1:size(X,1); end 

if  isempty(Y); DoY = 0;  else   DoY    = 1; end
try ttitle = varargin{3}; catch; ttitle = []; end

f = @(x,s,m,n)squeeze(s(x(:,m,n)));

p = size(X,2);
b = 0;
l = 3;

plottype = 'plot';%'scatter'; % plot / scatter

for m = 1:p
    for n = 1:p
        if n >= m
            
        if m == n; cl = 'g'; 
        else       cl = 'b';
        end
        b = b + 1;
        subplot(p,p,b);
        switch plottype
            case 'plot'
                plot(Hz,f(X,s,m,n),cl,'LineWidth',l); hold on;
                if DoY;
                    plot(Hz,f(Y,s,m,n),'r','LineWidth',l);
                end
                hold off;
                xlim([Hz(1) Hz(end)]);
            
            case 'scatter'
                if DoY;
                    scatter(f(X,s,m,n),f(Y,s,m,n));
                end
        end
        
        if ~isempty(ttitle) && b == 1;
            title(ttitle,'fontsize',16);
        end
        
        else
            b = b + 1;
        end
    end
end
        

% Esimates:
%--------------------------------------------------------------------------
% DCM.dtf                   - directed transfer functions (source space)
% DCM.ccf                   - cross covariance functions (source space)
% DCM.coh                   - cross coherence functions (source space)
% DCM.fsd                   - specific delay functions (source space)
% DCM.pst                   - peristimulus time
% DCM.Hz                    - frequency
%
% DCM.Ep                    - conditional expectation
% DCM.Cp                    - conditional covariance
% DCM.Pp                    - conditional probability
% DCM.Hc                    - conditional responses (y), channel space
% DCM.Rc                    - conditional residuals (y), channel space
% DCM.Hs                    - conditional responses (y), source space
% DCM.Ce                    - eML error covariance
% DCM.F                     - Laplace log evidence
% DCM.ID                    -  data ID