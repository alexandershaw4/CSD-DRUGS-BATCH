
xY = DCM.xY.y{1};

xY = VecRetainDim(xY,1);

y  = PEig(xY,1:6);

% dY = PEig(xY,1:size(xY,2));

% dY = spm_unvec(dY',xY);

% dY = spm_unvec(dY,DCM.xY.y{1});