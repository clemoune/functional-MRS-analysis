%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script normalises the timecourses by block.
% It needs the total time courses found in the 
% "group/metabolites/mouse#/matlab_raw_timecourse" and it generates 
% the A2 files in "group/metabolites/matlab_norm_timecourses"
% folder.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load full time courses
% FUNCTIONAL
currentdir=pwd;
group='functional';

for mouse_number=1:10
directory_fid=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'mouse',num2str(mouse_number),filesep);
directory_serie=strcat(directory_fid,filesep,'matlab_raw_timecourse',filesep);
load(strcat(directory_serie,'serie_moving_total_timecourse.mat'))
functional_exp(mouse_number).tttc=fft(Serie_moving_av');
end 

% CONTROL
group='control';

for mouse_number=1:10
directory_fid=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'mouse',num2str(mouse_number),filesep);
directory_serie=strcat(directory_fid,filesep,'matlab_raw_timecourse',filesep);
load(strcat(directory_serie,'serie_moving_total_timecourse.mat'))
control_exp(mouse_number).tttc=fft(Serie_moving_av');
end 

%% Normalising each block

block_rep=136 

for mouse_number=1:10
Norm_Acq_func=[];
Norm_Acq_control=[];
for block=1:4
Norm_Acq_func=[Norm_Acq_func functional_exp(mouse_number).tttc(:,[1:block_rep]+(block-1)*block_rep)/mean(mean(abs(functional_exp(mouse_number).tttc(:,[1:block_rep]+(block-1)*block_rep))))];
Norm_Acq_control=[Norm_Acq_control control_exp(mouse_number).tttc(:,[1:block_rep]+(block-1)*block_rep)/mean(mean(abs(control_exp(mouse_number).tttc(:,[1:block_rep]+(block-1)*block_rep))))];

end
Norm_Acq_func=[Norm_Acq_func functional_exp(mouse_number).tttc(:,[1:block_rep-1]+4*block_rep)/mean(mean(abs(functional_exp(mouse_number).tttc(:,[1:block_rep-1]+4*block_rep))))];
Norm_Acq_control=[Norm_Acq_control control_exp(mouse_number).tttc(:,[1:block_rep-1]+4*block_rep)/mean(mean(abs(control_exp(mouse_number).tttc(:,[1:block_rep-1]+4*block_rep))))];

functional_norm(mouse_number).tttc=Norm_Acq_func;
control_norm(mouse_number).tttc=Norm_Acq_control;

end


%% Outlier correction (baseline)
%Correcting for very drifted acquisitions (putting another block of the same mouse to replace it)
control_norm(6).tttc(:,[1:block_rep]+2*block_rep)=control_norm(6).tttc(:,[1:block_rep]+1*block_rep);
control_norm(8).tttc(:,[1:block_rep])=control_norm(8).tttc(:,[1:block_rep]+1*block_rep);

%% Saving .mat
group='functional';
directory_norm=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_norm_timecourse',filesep)
mkdir(directory_norm);
save(strcat(directory_norm,'A2_normalised_timecouse.mat'), 'functional_norm', '-mat')

group='control';
directory_norm=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_norm_timecourse',filesep)
mkdir(directory_norm);
save(strcat(directory_norm,'A2_normalised_timecouse.mat'), 'control_norm', '-mat')

