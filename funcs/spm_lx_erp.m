function [L] = spm_lx_erp(P,dipfit)
% observer matrix for a neural mass model: y = G*x
% FORMAT [G] = spm_lx_erp(P,dipfit)
% FORMAT [G] = spm_lx_erp(P,M)
%
% M.dipfit - spatial model specification
%
% G        - where y = L*x; G = dy/dx
% x        - state vector
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_erp.m 5939 2014-04-06 17:13:50Z karl $

% extract dipfit from model if necessary
%--------------------------------------------------------------------------
if  isfield(dipfit,'dipfit'), dipfit = dipfit.dipfit; end
if ~isfield(dipfit,'type'),   dipfit = 'LFP';         end

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L       = spm_erp_L(P,dipfit);               % lead field per source

% L       = kron(P.J,L);                       % lead-field per state

try
    
l = kron(P.J(1,:),L);

L1 = kron(P.J(1,:),L); % IFG's
L2 = kron(P.J(2,:),L); % STG's
L3 = kron(P.J(3,:),L); % A1's

L = l*0;

L(1,1) = L1(1,1); % SS,lIFG
L(2,2) = L2(2,2); % SS lSTG
L(3,3) = L3(3,3);  % SS lA1
L(4,4) = L1(4,4);  % SS rIFG
L(5,5) = L2(5,5);  % SS rSTG
L(6,6) = L3(6,6);  % SS rA1
L(1,7) = L1(1,7);  % SP
L(2,8) = L2(2,8);  % SP 
L(3,9) = L3(3,9);  % SP
L(4,10) = L1(4,10); % SP
L(5,11) = L2(5,11); % SP
L(6,12) = L3(6,12); % SP
L(1,19) = L1(1,19); % DP  
L(2,20) = L2(2,20); % DP
L(3,21) = L3(3,21); % DP
L(4,22) = L1(4,22); % DP
L(5,23) = L2(5,23); % DP
L(6,24) = L3(6,24); % DP

catch
    L       = kron(P.J,L);
end
    




% % Alex code: we want sep J [contributions] for each node but symmatrically
% % confined, so 3 nodes (IFG, STG & A1). The order of nodes is:
% % 1 = LIFG, 2 = LSTG, 3=LA1, 4 = RIFG, 5 = RSTG, 6 = RA1
% 
% L1 = kron(P.J(1,:),L); % IFG's
% L2 = kron(P.J(2,:),L); % STG's
% L3 = kron(P.J(3,:),L); % A1's
% 
% L       = l*0; % get sparse size
% L(1,1)  = L1(1,1);
% L(2,2)  = L2(2,2);
% L(3,3)  = L3(3,3);
% L(4,4)  = L1(4,4);
% L(5,5)  = L2(5,5);
% L(6,6)  = L3(6,6);
% L(1,13) = L1(1,13);
% L(2,14) = L2(2,14);
% L(3,15) = L3(3,15);
% L(4,16) = L1(4,16);
% L(5,17) = L2(5,17);
% L(6,18) = L3(6,18);
% L(1,37) = L1(1,37);
% L(2,38) = L2(2,38);
% L(3,39) = L3(3,39);
% L(4,40) = L1(4,40);
% L(5,41) = L2(5,41);
% L(6,42) = L3(6,42);