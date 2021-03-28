%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates illustrative singlets time courses in Figure 2.
% It computes linewidth and amplitude for the sum of the 3 singlet from the
% average time courses.
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


%% Computing linewidth and amplitude with and without BOLD correction 

bw=4000;
time=[1:2048];
mouse_number_functional=[1 3 5 6 8 9 10];
size_timecourse=135;
addpath(strcat(currentdir,filesep,'support_functions'))

% To generate the corrected BOLD courses, uncomment the following line to
% apply a line broadening (as in Figure 3)
% lb_spec_act=0.5;
lb_spec_act=0.0;


% hemodynamic response repetition identified on uncorrected water time courses
step=4;
init_bold=round([33 81 130 178 226]/step);
end_bold=round([50 101 149 199 246]/step);


vect_bold=[init_bold(1):end_bold(1) init_bold(2):end_bold(2) init_bold(3):end_bold(3) init_bold(4):end_bold(4) init_bold(5):end_bold(5)];
N_bol=size(vect_bold);

for mouse_number=mouse_number_functional
        func_corrected(mouse_number).timecourse=functional_block(mouse_number).timecourse;
        for k=1:N_bol(2)
        func_corrected(mouse_number).timecourse(vect_bold(k),:)=(functional_block(mouse_number).timecourse(vect_bold(k),:)).*exp(-lb_spec_act*time/bw);
        end
        
end        

% Computation of linewidth and amplitude



for mouse_number = mouse_number_functional

[linewidth_f,height_f] = course_lw(func_corrected(mouse_number).timecourse',size_timecourse);

linewidth_metab_f(mouse_number,:)=linewidth_f;
height_metab_f(mouse_number,:)=height_f;
end


N_linewidth_f=linewidth_metab_f(:,:)./mean(linewidth_metab_f(:,round(size_timecourse/2):end),2);
N_height_f=height_metab_f(:,:)./mean(height_metab_f(:,round(size_timecourse/2):end),2);
std_N_linewidth_f=std(N_linewidth_f,1);
std_N_height_f=std(N_height_f,1);

%% Generating the figures
% change lb_spec_act for BOLD correction

% Grey areas to identify the BOLD events
pos1=[round(33/step) 0.97 round(50/step+1)-round(33/step) 0.06];
pos2=[round(81/step) 0.97 round(101/step+1)-round(81/step) 0.06];
pos3=[round(130/step) 0.97 round(149/step+1)-round(130/step) 0.06];
pos4=[round(178/step) 0.97 round(199/step)-round(178/step) 0.06];
pos5=[round(226/step) 0.97 round(246/step)-round(226/step) 0.06];

%Generating plot for linewidth

subplot(2,1,2)
hold on 
rectangle('Position',pos1,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos2,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos3,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos4,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos5,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])

plot(mean(N_linewidth_f(mouse_number_functional,:),1),'Color',[0.8 0 0.2],'Linewidth',1.8)
shadedErrorBar([],mean(N_linewidth_f(mouse_number_functional,:),1),std_N_linewidth_f/sqrt(7),'lineProps',{'Color',[0.8 0 0.2]})

ax = gca; % current axes
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.96 1.04];
ax.XLim = [-20/4 550/4];
ax.XColor ='none';
ax.YColor ='k';
ax.FontSmoothing = 'on';
title(strcat('Mean linewidth of metabolites time course - lb applied=',num2str(lb_spec_act)))

%Generating plot for amplitude
subplot(2,1,1)
hold on
rectangle('Position',pos1,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos2,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos3,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos4,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos5,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
plot(mean(N_height_f(mouse_number_functional,:),1),'Color','k','Linewidth',1.8)
shadedErrorBar([],mean(N_height_f(mouse_number_functional,:),1),std_N_height_f/sqrt(7),'lineProps','k')


ax = gca; % current axes
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.98 1.04];
ax.XLim = [-20/4 550/4];
ax.XColor ='none';
ax.YColor ='k';
ax.FontSmoothing = 'on';
title(strcat('Mean amplitude of metabolites time course - lb applied=',num2str(lb_spec_act)))
