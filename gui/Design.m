function p = Design
% For use with ERP_DCM_ALEX.m
% This func for designing DCM experimental setup
%
% AS2016


p(1).xU.X     = 1;            % (ns x ns) design matrix - e.g. eye(8)
p(1).xU.name  = {'2s'};       % condition names
p(1).tCode    = 1;            % condition codes in SPM

% p(1).xU.X     = eye(8);            % roving
% p(1).xU.name  = {'Deviant' 'r1','2','3','4','5','6','7'}; % condition names
% p(1).tCode    = 1:8;             % condition codes in SPM


