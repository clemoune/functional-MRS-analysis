%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script forms the 13min time courses and apply appropriate BOLD correction.
% It needs the detrended time courses found in the 
% "group/metabolites/matlab_detrended_timecourse" folder and it generates 
% the 13min time courses the same folder (A4 files corrected & uncorrected).
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the detrended data

currentdir=pwd;
addpath(strcat(currentdir,filesep,'support_functions'))
%Functional
group='functional';
directory_detrended_func=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
load(strcat(directory_detrended_func,'A3_detrended_timecouse.mat'),'func_res_corr2')

%Control
group='control';
directory_detrended_cont=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
load(strcat(directory_detrended_cont,'A3_detrended_timecouse.mat'),'control_res_corr2')


%% Forming 13 min time courses summing the 5 blocks

size_block=135;

for mouse_number=1:10
for i=1:size_block
reshaped_timecourse_c=[control_res_corr2(mouse_number).tttc(:,i) control_res_corr2(mouse_number).tttc(:,i+136) control_res_corr2(mouse_number).tttc(:,i+136*2) control_res_corr2(mouse_number).tttc(:,i+136*3) control_res_corr2(mouse_number).tttc(:,i+136*4)]; 
[spec_water fid_phased_c fid_raw STD_PHI]=phaser(reshaped_timecourse_c); 
control_block(mouse_number).timecourse(i,:)=fid_phased_c;

% % uncomment to save the fid_asc control files to be used with LCModel (no need of BOLD corr)
% group='control'; 
% directory_LCModel_control=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'LCModel',filesep,'mouse',num2str(mouse_number,filesep,num2str(i),filesep);
% mkdir(directory_LCModel_control)
% save_fid(fid_phased_c',0,directory_LCModel_control);

reshaped_timecourse_f=[func_res_corr2(mouse_number).tttc(:,i) func_res_corr2(mouse_number).tttc(:,i+136) func_res_corr2(mouse_number).tttc(:,i+136*2) func_res_corr2(mouse_number).tttc(:,i+136*3) func_res_corr2(mouse_number).tttc(:,i+136*4)]; 
[spec_water fid_phased_f fid_raw STD_PHI]=phaser(reshaped_timecourse_f); 
functional_block(mouse_number).timecourse(i,:)=fid_phased_f;

end
  
end

save(strcat(directory_detrended_cont,'A4_TimeCourse13min.mat'),'control_block');
save(strcat(directory_detrended_func,'A4_TimeCourse13min_not_BOLD_corrected.mat'),'functional_block')


close all

%% Applying BOLD correction to functional time courses (line broadening on BOLD events)

bw=4000;
time=[1:2048];
lb_spec_act=[0 2 0.4 1.3 0.6 0.7 2.5 0.5 0.7 0.7];


% hemodynamic response repetition identified on water time courses
init_bold=[33 81 130 178 226];
end_bold=[50 101 149 199 246];
step=4;

vect_bold=[round(init_bold(1)/step):round(end_bold(1)/step) round(init_bold(2)/step):round(end_bold(2)/step) round(init_bold(3)/step):round(end_bold(3)/step) round(init_bold(4)/step):round(end_bold(4)/step) round(init_bold(5)/step):round(end_bold(5)/step)];
N_bol=size(vect_bold);

for mouse_number=1:10
        functional_block_bold_corr(mouse_number).timecourse=functional_block(mouse_number).timecourse;
        for k=1:N_bol(2)
        functional_block_bold_corr(mouse_number).timecourse(vect_bold(k),:)=(functional_block(mouse_number).timecourse(vect_bold(k),:)).*exp(-lb_spec_act(mouse_number)*time/bw);
        end
        
% % uncomment to save the fid_asc control files to be used with LCModel (no need of BOLD corr)
% group='functional'; 
% for j=1:size(functional_block(mouse_number).tttc(:,:),1)
% directory_LCModel_func=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'LCModel',filesep,'mouse',num2str(mouse_number,filesep,num2str(i),filesep);
% mkdir(directory_LCModel_func)
% save_fid(functional_block_bold_corr(mouse_number).timecourse(j,:)',0,directory_LCModel_func);
% end       

end 
save(strcat(directory_detrended_func,'A4_TimeCourse13min_BOLD_corrected.mat'),'functional_block_bold_corr')


