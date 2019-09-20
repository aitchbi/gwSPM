% This function is part of the toolbox:
%       gwSPM: Graph-based, Wavelet-based Statistical Parametric Mapping
%       (v1.00)
%
% 	Author: Hamid Behjat
% 
%   Biomedical Signal Processing Group, 
%   Dept. of Biomedical Engineering,
%   Lund University, Sweden
% 
%   June 2016
%
function [tMx,tMn,option] = gwspm_estimate_computeTime(trans,n1,n2,option)

if isempty(n1)
    n1 = 60;
end

if isempty(n2)
    n2 =100; 
end

if isempty(option)
    if gwspm_check_par([])
        option = 'parallel';
    else
        option = 'sequential';
    end
    switchOff = 1;
else
    switchOff = 0;
end

[~, t1] = gwspm_construct_atoms(trans,1,n1,'cbr',option);
[~, t2] = gwspm_construct_atoms(trans,1,n2,'cbl',option);

t = (1+trans.wav_scales)*...
    (t1.single*numel(trans.cbr.indices)+...
    t2.single*numel(trans.cbl.indices));
tMx = t/3600;
tMn = 0.9*t/3600; 

if strcmp(option,'parallel') && switchOff
    delete(gcp('nocreate'))
end
