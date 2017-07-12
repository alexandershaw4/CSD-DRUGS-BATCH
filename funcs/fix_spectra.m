

nc = length(DCM.xY.y);
Q  = @squeeze;
Hz = DCM.xY.Hz;

for ic = 1:nc
    Y  = DCM.xY.y{ic};

    % make some autospectra
    for ns = 1:size(Y,2)
        Y(:,ns,ns) = TSNorm(Q(Y(:,ns,ns)));
        Y(:,ns,ns) = Q(Y(:,ns,ns)) + (1./Hz.^2)';
        Y(:,ns,ns) = HighResMeanFilt(Q(Y(:,ns,ns)),1,4);
    end
    
    % then calc cross funcs
    for nsx = 1:size(Y,2)
        for nsy = 1:size(Y,3)
            if nsx ~= nsy
                Y(:,nsx,nsy) = Y(:,nsx,nsx) .* conj(Y(:,nsy,nsy));
            end
        end
    end
    DCM.xY.y{ic} = Y;

end

