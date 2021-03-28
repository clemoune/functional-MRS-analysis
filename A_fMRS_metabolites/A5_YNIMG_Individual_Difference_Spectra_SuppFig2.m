%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the individual difference spectra, used to
% determine outliers (Supplementary Figure 2). 
% Functional exp: outliers = mice 2, 4, 7
% Control exp: outliers = mice 5, 8, 10
% It needs the 13min time courses (A4 files) found in the
% "group/metabolites/matlab_detrended_timecourse".
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
%% FUNCTIONAL: individual difference spectra with and without BOLD correction

bw=4000;
time=[1:2048];

% A line broadening is applied to the difference spectra for visualisation
lb_diff=4;

% A frequency correction is applied to minimise the residuals (not to be mistaken by a BOLD effect)
freq_func=-[3.2 3.7 0.4 1.9 -7.2 1.5 2.6 -5 1.1 3.45];

vector=-10*time/2048+9.69;%%ppm_scaling


for mouse_number=1:10
        
[spec_water TOTAL_ACT fid_raw STD_PHI]=phaser(functional_block_bold_corr(mouse_number).timecourse(1:67,:)');
[spec_water TOTAL_ACT_NC fid_raw STD_PHI]=phaser(functional_block(mouse_number).timecourse(1:67,:)');

[spec_water TOTAL_REST fid_raw STD_PHI]=phaser(functional_block_bold_corr(mouse_number).timecourse(69:135,:)');

diff_spec_bold_corr=real(fftshift(fft((((TOTAL_ACT).*exp(1i*freq_func(mouse_number)*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));
diff_spec_no_corr=real(fftshift(fft((((TOTAL_ACT_NC).*exp(1i*freq_func(mouse_number)*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));

figure
hold on 
plot(flip(vector),zeros(2048,1),'--','LineWidth',2,'Color',[0 0 1])
plot(flip(vector),diff_spec_bold_corr,'k','LineWidth',0.8)
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [-40 150];
ax.XDir = 'reverse';
ax.XLim = [0.5 4.];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title(strcat('Difference spectrum with BOLD corr mouse',num2str(mouse_number)))
v = [0.5 -35;  0.5 65+65; 1.5 65+65; 1.5 -35];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.8,'EdgeColor','none')

figure
hold on 
plot(flip(vector),zeros(2048,1),'--','LineWidth',2,'Color',[0 0 1])
plot(flip(vector),diff_spec_no_corr','k','LineWidth',0.8)
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [-40 150];
ax.XLim = [0.5 4.];
ax.XDir = 'reverse';
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title(strcat('Difference spectrum without BOLD corr mouse',num2str(mouse_number)))
v = [0.5 -35;  0.5 65+65; 1.5 65+65; 1.5 -35];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.8,'EdgeColor','none')
end
%% CONTROL: individual difference spectra with and without BOLD correction

bw=4000;
time=[1:2048];

% A line broadening is applied to the difference spectra for visualisation
lb_diff=4;

% A frequency correction is applied to minimise the residuals (not to be mistaken by a BOLD effect)
freq_control=-[0.0 -2.2 1.0 1.1 13.9 -0.7 -4.2 4 -7.8 2.5];
freq_func=-[3.2 3.7 0.4 1.9 -7.2 1.5 2.6 -5 1.1 3.45];

vector=-10*time/2048+9.69;%%ppm_scaling

for mouse_number=1:10
             
[spec_water TOTAL_ACT_NC fid_raw STD_PHI]=phaser(control_block(mouse_number).timecourse(1:67,:)');

[spec_water TOTAL_REST fid_raw STD_PHI]=phaser(control_block(mouse_number).timecourse(69:135,:)');

diff_spec=real(fftshift(fft((((TOTAL_ACT_NC).*exp(1i*freq_control(mouse_number)*time/bw)-TOTAL_REST).*exp(-lb_diff*time/bw))')));

figure
hold on 
plot(flip(vector),zeros(2048,1),'--','LineWidth',2,'Color',[0 0 1])
plot(flip(vector),diff_spec,'k','LineWidth',0.8)
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [-40 150];
ax.XDir = 'reverse';
ax.XLim = [0.5 4.];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title(strcat('Difference spectrum for control mouse',num2str(mouse_number)))
v = [0.5 -35;  0.5 65+65; 1.5 65+65; 1.5 -35];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.8,'EdgeColor','none')

end

