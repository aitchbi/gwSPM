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
function [mask_adj, indAdj]=gwspm_adjust_mask(mask,numClean,connClean,name) 
% numClean : max size of isolated componnets to be removed
% connClean: connectivity to use for clean operation

indG=find(mask);
maskBW=zeros(size(mask));
maskBW(indG)=1;

% Remove isolated componets of size 'numClean' voxels or less
GMBW2=bwareaopen(maskBW,numClean,connClean);

% Change back to grayscale
indAdj=find(GMBW2);
mask_adj=zeros(size(mask));
mask_adj(indAdj)=mask(indAdj);
sprintf('%d voxels were removed for adjusting the %s mask.',numel(indG)-numel(indAdj),name)

