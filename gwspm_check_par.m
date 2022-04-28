function chk = gwspm_check_par(opts)
% opts=='dlt': shut down opened pool
% opts=[]: keep opened pool open

if isempty(opts)
    fprintf('Starting paralle pool if avaialble.. \n');
else
    fprintf('Checking if paralle pool is avaialble.. \n');
end
try
    obj = parpool();
    chk = obj.Connected;
    if isequal(opts,'dlt')
        delete(gcp('nocreate'))
    end
catch
    try
        delete(gcp('nocreate'))
        obj = parpool();
        chk = obj.Connected;
        if isequal(opts,'dlt')
            delete(gcp('nocreate'))
        end
    catch
        fprintf(' parpool not available.\n');
    end
end

