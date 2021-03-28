%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the mean difference spectra for each group (Figure 5)
% It needs the 13min time courses (A4 files) found in the
% "group/metabolites/matlab_detrended_timecourse". 
% You can jump to l.61 to skip the summing.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the detrended + BOLD corrected data

currentdir=pwd;
addpath(strcat(currentdir,filesep,'support_functions'))
%Functional
group='functional';
directory_detrended_func=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
load(strcat(directory_detrended_func,'A4_TimeCourse13min_not_BOLD_corrected.mat'))
load(strcat(directory_detrended_func,'A4_TimeCourse13min_BOLD_corrected.mat'))

%Control
group='control';
directory_detrended_cont=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
load(strcat(directory_detrended_cont,'A4_TimeCourse13min.mat'))

%% Generating the average spectra

mouse_number_functional=[1 3 5 6 8 9 10];
mouse_number_control=[1 2 3 4 6 7 9];
size_block=135;

for j=1:size_block
FID_tot_uncorr=[];
FID_tot_corr=[];
FID_tot_control=[];

%Functional
    for mouse_number=mouse_number_functional   
    FID_tot_uncorr=[FID_tot_uncorr functional_block(mouse_number).timecourse(j,:)']; 
    FID_tot_corr=[FID_tot_corr functional_block_bold_corr(mouse_number).timecourse(j,:)'];
    end
    
    [spec_water fid_phased_uncorr fid_raw STD_PHI]=phaser(FID_tot_uncorr);
    Average_functional_uncorr(j,:)=fid_phased_uncorr;
    
   
    [spec_water fid_phased_corr fid_raw STD_PHI]=phaser(FID_tot_corr);
    Average_functional_BOLD_corr(j,:)=fid_phased_corr;  
    
%Control
    for mouse_number=mouse_number_control    
    FID_tot_control=[FID_tot_control control_block(mouse_number).timecourse(j,:)'];  
    end
    [spec_water fid_phased_control fid_raw STD_PHI]=phaser(FID_tot_control);
    Average_control(j,:)=fid_phased_control;
   
    
end

save(strcat(directory_detrended_cont,'A6_Average_control.mat'),'Average_control','-mat')
save(strcat(directory_detrended_func,'A6_Average_functional_BOLD_corr.mat'),'Average_functional_BOLD_corr','-mat')
save(strcat(directory_detrended_func,'A6_Average_functional_uncorr.mat'),'Average_functional_uncorr','-mat')

%% FUNCTIONAL: generating the difference spectra
  
load(strcat(directory_detrended_func,'A6_Average_functional_uncorr.mat'))
load(strcat(directory_detrended_func,'A6_Average_functional_BOLD_corr.mat'))

bw=4000;
time=[1:2048];
freq=1;
phi_spec=-0.5;
lb_diff=4;
lb_act=0.00; % to play with residual BOLD / shim drift
vector=-10*time/2048+9.69;%%ppm_scaling

% Generating average spectra during ACTIVE and RECOVERY periods
[spec_water TOTAL_ACT fid_raw STD_PHI]=phaser(Average_functional_BOLD_corr(1:67,:)');
[spec_water TOTAL_ACT_NC fid_raw STD_PHI]=phaser(Average_functional_uncorr(1:67,:)');
[spec_water TOTAL_REST fid_raw STD_PHI]=phaser(Average_functional_BOLD_corr(69:135,:)');

diff_spec_bold_corr=real(fftshift(fft((((TOTAL_ACT).*exp(1i*freq*time/bw).*exp(-lb_act*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));
diff_spec_no_corr=real(fftshift(fft((((TOTAL_ACT_NC).*exp(1i*freq*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));

vector=-10*time/2048+9.69;%%ppm_scaling

figure
hold on 
plot(flip(vector),zeros(2048,1)+10,'--','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(flip(vector),diff_spec_bold_corr/sqrt(7)+10,'k','LineWidth',0.8)
plot(flip(vector),zeros(2048,1)+60+10,'--','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(flip(vector),diff_spec_no_corr/sqrt(7)+60+10,'k','LineWidth',0.8)
plot(flip(vector),real(fftshift(fft(TOTAL_REST'*exp(1i*phi_spec))))/100+70+60+10,'Color',[0.75 0.75 0.75],'LineWidth',1)
plot(flip(vector),real(fftshift(fft(TOTAL_ACT'*exp(1i*phi_spec))))/100+100+60+10,'Color',[0.5 0.5 0.5],'LineWidth',1)

ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [-40 300+10];
ax.XLim = [0.5 4.];
ax.XDir = 'reverse';
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 0.9;
title('Mean Difference Functional')
v = [0.5 -30;  0.5 220+10; 1.5 220+10; 1.5 -30];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.6,'EdgeColor','none')

%% CONTROL: generating the difference spectra
  
load(strcat(directory_detrended_cont,'A6_Average_control.mat'))
bw=4000;
time=[1:2048];
freq=1.55;
phi_spec=-0.5;
lb_diff=4;
lb_act=0.00; % to play with residual BOLD / shim drift

% Generating average spectra during NO STIM 1 and NO STIM 2 periods
[spec_water TOTAL_ACT fid_raw STD_PHI]=phaser(Average_control(1:67,:)');
[spec_water TOTAL_REST fid_raw STD_PHI]=phaser(Average_control(69:135,:)');

diff_spec_control=real(fftshift(fft((((TOTAL_ACT).*exp(1i*freq*time/bw).*exp(-lb_act*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));

vector=-10*time/2048+9.69;%%ppm_scaling

figure
hold on 
plot(flip(vector),zeros(2048,1)+60+10,'--','LineWidth',1,'Color',[0.5 0.5 0.5])
plot(flip(vector),diff_spec_control/sqrt(7)+60+10,'k','LineWidth',0.8)
plot(flip(vector),real(fftshift(fft(TOTAL_REST'*exp(1i*phi_spec))))/100+70+60+10,'Color',[0.75 0.75 0.75],'LineWidth',1)
plot(flip(vector),real(fftshift(fft(TOTAL_ACT'*exp(1i*phi_spec))))/100+100+60+10,'Color',[0.5 0.5 0.5],'LineWidth',1)

ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [-40 300+10];
% ax.YLim = [-150 170];
ax.XLim = [0.5 4.];
ax.XDir = 'reverse';
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 0.9;
title('Mean Difference Control')
v = [0.5 -30;  0.5 220+10; 1.5 220+10; 1.5 -30];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.6,'EdgeColor','none')