function [tMx,tMn,opts] = gwspm_estimate_computeTime(trans,n1,n2,opts)

if isempty(n1)
    n1 = 60;
end

if isempty(n2)
    n2 =100; 
end

if isempty(opts)
    if gwspm_check_par([])
        opts = 'parallel';
    else
        opts = 'sequential';
    end
    SwitchOff = 1;
else
    SwitchOff = 0;
end

[~,t1] = gwspm_construct_atoms(trans,1,n1,'cbr',opts);
[~,t2] = gwspm_construct_atoms(trans,1,n2,'cbl',opts);

t = (1+trans.wav_scales)*...
    (t1.single*numel(trans.cbr.indices)+...
    t2.single*numel(trans.cbl.indices));
tMx = t/3600;
tMn = 0.9*t/3600; 

if strcmp(opts,'parallel') && SwitchOff
    delete(gcp('nocreate'))
end
