function [y,Ep,W] = dcm_contributions_csd(DCM,Param,k,dofit)



% step
%------------------------------------------------------------------------
f = eval(['@' DCM.M.IS]);
x = {DCM.Ep,DCM.M,DCM.xU};


% opts & step
%------------------------------------------------------------------------
if nargin < 4
    dofit = 0;    % use this approach to actualy fit the model
end
if nargin < 3
    k     = 2;
end

% fields / indices 
%------------------------------------------------------------------------
if nargin < 2;
    ind   = find(spm_vec(DCM.M.pC)); % all whose prior variance is > 0
else
    if iscell(Param); % find specified parameters
        ind     = [];
        for p   = 1:length(Param);
            ind = [ind ; spm_fieldindices(DCM.Ep,Param{p})];
        end
    elseif ischar(Param)
           ind = spm_fieldindices(DCM.Ep,Param);
    end
end
try ind;
catch
    ind = find(spm_vec(DCM.M.pC));
end

% initial (current) state & data
%------------------------------------------------------------------------
Y     = f(x{:});

if dofit;
    H = DCM.Hc{:} + DCM.Rc{:};
    E = err_f(H,Y);
end

% do increases
%------------------------------------------------------------------------
for i = 1:length(ind)
    P = spm_vec(DCM.Ep);
    P(ind(i)) = P(ind(i)) + k;
    x{1} = spm_unvec(P , DCM.Ep);
    y{i} = f(x{:});
    
    if dofit;
        % score
        e(i)  = err_f(y{i},H);
        pE(i) = x{1}; 
        plot(1:i,log(e),':r');
        title(sprintf('Iteration %d',i));drawnow;
        
        % criteria
        if ( i == 1 && e(i) < E ) || ( (i > 1) && e(i) < e(i-1) )
            DCM.Ep = x{1};
            fprintf('Model fit improved on iteration %d,updating...\n',i);
        end
    
    end
end

if dofit;
    
    [v,W] = min(e);
    Ep    = pE(W);
    
    PCSD(H,y{W}{:},[],DCM.Hz);

else
    Ep = [];
    W  = [];
end


end

function e = err_f(y,H)

% Error function
v = @spm_vec;
e = (sum (v(H) - v(y) ).^2 );

end
