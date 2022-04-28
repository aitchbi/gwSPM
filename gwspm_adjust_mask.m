function [mask_out,ind_out]=gwspm_adjust_mask(mask,p,conn,name) 
% Remove isolated componets of size p voxels or less from input mask based
% on neighborhood connectivity as specified by conn.

ind_in = find(mask);

v0 = zeros(size(mask));

mask_bw = v0;
mask_bw(ind_in) = 1;

mask_bw2 = bwareaopen(mask_bw,p,conn);
ind_out = find(mask_bw2);

mask_out = v0;
mask_out(ind_out) = mask(ind_out);

sprintf('%d voxels removed for cleaning %s mask.',...
    length(ind_in)-length(ind_out),name)

