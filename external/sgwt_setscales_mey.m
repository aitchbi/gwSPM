% sgwt_setscales_mey : Compute a set of wavelet scales adapted to spectrum bounds
%
% function nsgwt_setscales_mey(lmax,Nscales,varargin)
%
% returns a set of dyadically spaced scales btw 0 and max
%
% Inputs:
% lmax - maximum eigenvalue of Laplacian.
% Nscales - # of wavelet scales
% t1 - upper boundary of wavelets
%
% Outputs:
% s - wavelet scales
%
% 2011, Nora Leonardi, Dimitri Van De Ville


%% OLD ...
% function s=sgwt_setscales_mey(lmax,Nscales,varargin)
%     control_params={'t1',1};
%     argselectAssign(control_params);
%     argselectCheck(control_params,varargin);
%     argselectAssign(varargin);
%
%     t11=4*t1/3;
%
%     smin=t11/lmax;
%
%     s=zeros(1,Nscales);
%     for k=1:Nscales
%         s(k)=smin*2^(Nscales-k);
%     end
% end
%%

function s=sgwt_setscales_mey(lmax,Nscales,hamShift,varargin)
control_params={'t1',1};
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

%%
t11=2*t1/3;

smin=hamShift*t11/lmax;
%fprintf('Changed XXX sgwt_setscales_mey\n');
% smax=t11/lmax*17; % the larger the more towards center

s=zeros(1,Nscales);
for k=1:Nscales
    %s(k)=smax*2^(-k); % start from bottom
    s(k)=smin*2^(Nscales+1-k); % start from top
end
end
