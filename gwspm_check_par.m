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
function check = gwspm_check_par(option)
% option: 'dlt' to shut down opened pool, [] to keep open

fprintf('Checking whether paralle pool is avaialble.. \n');
try
    obj = parpool();
    check = obj.Connected;
    fprintf('Yes! \n');
    if strcmp(option,'dlt')
        delete(gcp('nocreate'))
    end
catch
    try
        delete(gcp('nocreate'))
        obj = parpool();
        check = obj.Connected;
        fprintf('Yes! \n');
        if strcmp(option,'dlt')
            delete(gcp('nocreate'))
        end
    catch
        fprintf('No :( \n');
    end
end

