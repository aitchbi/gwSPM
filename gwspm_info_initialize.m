function info = gwspm_info_initialize(I)
info.contrastDirName = 'gwspm_import_contrasts';
info.templateDirName = 'gwspm_templates';
info.templateDirRoot =  fileparts(deblank(I(1,:)));
info.templateName = 'w1Template_6.nii,1';
info.maskName = 'GM';
