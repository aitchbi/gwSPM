function chk = gwspm_check_xCon(SPM,name,wCon,descrip)

chk=0;
for iter=1:length(SPM.xCon)
    tmp=strfind(SPM.xCon(iter).name,' T('); 
    if numel(tmp)==0
        continue
    elseif tmp-1~=length(name)
        continue
    elseif ~strcmp(SPM.xCon(iter).name(1:tmp(end)-1),name)
        continue
    elseif numel(strfind(SPM.xCon(iter).name,descrip))==0
        continue
    elseif ~isfield(SPM.gWavelet.wCon(iter),'mask')
        SPM.gWavelet.wCon(iter).mask=[];
    end
    if SPM.gWavelet.wCon(iter).typeI==wCon.typeI && ...
            SPM.gWavelet.wCon(iter).siglevel==wCon.siglevel && ...
            (strcmp(SPM.gWavelet.wCon(iter).mask,wCon.mask)==1 || ...
            length([wCon.mask SPM.gWavelet.wCon(iter).mask])==0) %#ok<ISMT>
        chk=iter;
        break;
    end
end