%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates water time courses upon stimulation and generates
% the material of Figure 2 (for identification of BOLD events).
% It needs the raw data for water found each mouse# folder.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading water data

currentdir=pwd;
addpath(currentdir,filesep,'support_functions')
ACQ=[19 19 20 20 20 21 20 19 19 21];
Nb_Repet=272;
WATER_TOT=[];
for mouse_number=1:10
[WATER]=load_array_FID2(strcat(currentdir,filesep,'mouse',num2str(mouse_number),filesep,num2str(ACQ(mouse_number)),filesep,'fid'),Nb_Repet);
Water_Functional(mouse_number).timecourse=WATER;
WATER_TOT=[WATER_TOT WATER];
end

%% Computing linewidth and amplitude with and without BOLD correction 
% it also generates the figures. If you want to check individual
% timecourses see last section.
bw=4000;
time=[1:2048];

% Mice 3 and 10 have signal intensity variation, and mouse 7's first point is too intense. 
mouse_water_functional=[1 2 4 5 6 8 9];
% To generate the corrected BOLD courses, uncomment the following line to
% apply a line broadening (as in Figure 2)
% lb_spec_act=0.4;
lb_spec_act=0.0;


% hemodynamic response repetition identified on uncorrected water time courses
init_bold=[33 81 130 178 226];
end_bold=[50 101 149 199 246];
step=4;

vect_bold=[init_bold(1):end_bold(1) init_bold(2):end_bold(2) init_bold(3):end_bold(3) init_bold(4):end_bold(4) init_bold(5):end_bold(5)];
N_bol=size(vect_bold);
WATER_TOT_c=[];

for mouse_number=mouse_water_functional
        Water_corrected(mouse_number).timecourse=Water_Functional(mouse_number).timecourse;
        for k=1:N_bol(2)
        Water_corrected(mouse_number).timecourse(:,vect_bold(k))=(Water_Functional(mouse_number).timecourse(:,vect_bold(k))').*exp(-lb_spec_act*time/bw);
        end
        WATER_TOT_c=[WATER_TOT_c Water_corrected(mouse_number).timecourse];
        
end        

% Computation of linewidth and amplitude
linewidth_c=zeros(size(mouse_water_functional,2),Nb_Repet);
height_c=zeros(size(mouse_water_functional,2),Nb_Repet);
integration_c=zeros(size(mouse_water_functional,2),Nb_Repet);
linewidth_c_norm=zeros(size(mouse_water_functional,2),Nb_Repet);
height_c_norm=zeros(size(mouse_water_functional,2),Nb_Repet);

for i =1:size(mouse_water_functional,2)
    mat_file_c=WATER_TOT_c(:,(i-1)*Nb_Repet+1:(i-1)*Nb_Repet+Nb_Repet);
   for j=1:Nb_Repet
    zero_fill_c=abs(fftshift(fft([mat_file_c(:,j)' zeros(4096,1)']))); 
    
    [M I]=max(zero_fill_c);
    
    half_width=M/2;
    k=I;
    while zero_fill_c(1,k)>half_width
        k=k+1;
    end 
    
    k_2=I;
    while zero_fill_c(1,k_2)>half_width
        k_2=k_2-1;
    end
linewidth_c(i,j)=k-k_2;    
height_c(i,j)=M;
integration_c(i,j)=sum(zero_fill_c);
   end
linewidth_c_norm(i,:)=linewidth_c(i,:)/mean(linewidth_c(i,:),2);    
height_c_norm(i,:)=height_c(i,:)/mean(height_c(i,:),2);
integration_c_norm(i,:)=integration_c(i,:)/mean(integration_c(i,:),2);    
   
end


M_linewidth_c=mean(linewidth_c_norm(1:size(mouse_water_functional,2),:),1);
M_height_c=mean(height_c_norm(1:size(mouse_water_functional,2),:),1);
std_M_linewidth_c=std(linewidth_c_norm(1:size(mouse_water_functional,2),:),1);
std_M_height_c=std(height_c_norm(1:size(mouse_water_functional,2),:),1);
M_integration_c=mean(integration_c(1:size(mouse_water_functional,2),:),1);

% Average over a phase-cycle (= 4 repetitions)
linewidth_phasecycle_c=zeros;
height_phasecycle_c=zeros;
integration_phasecycle_c=zeros;
step=4;

for i=1:1:272/step
   linewidth_phasecycle_c(1,i)=mean(M_linewidth_c(1,(i-1)*step+1:(i-1)*step+step));
   height_phasecycle_c(1,i)=mean(M_height_c(1,(i-1)*step+1:(i-1)*step+step)); 
   std_linewidth_phasecycle_c(1,i)=mean(std_M_linewidth_c(1,(i-1)*step+1:(i-1)*step+step));
   std_height_phasecycle_c(1,i)=mean(std_M_height_c(1,(i-1)*step+1:(i-1)*step+step));
   integration_phasecycle_c(1,i)=mean(M_integration_c(1,(i-1)*step+1:(i-1)*step+step)); 
end

% Grey areas to identify the BOLD events
pos1=[round(33/step) 0.982 round(50/step+1)-round(33/step) 1.018-0.982];
pos2=[round(81/step) 0.982 round(101/step+1)-round(81/step) 1.018-0.982];
pos3=[round(130/step) 0.982 round(149/step+1)-round(130/step) 1.018-0.982];
pos4=[round(178/step) 0.982 round(199/step)-round(178/step) 1.018-0.982];
pos5=[round(226/step) 0.982 round(246/step)-round(226/step) 1.018-0.982];

%Generating plot for linewidth
figure
hold on 
rectangle('Position',pos1,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos2,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos3,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos4,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos5,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])

plot(linewidth_phasecycle_c,'Color',[0.8 0 0.2],'Linewidth',1.8)
shadedErrorBar([],linewidth_phasecycle_c,std_linewidth_phasecycle_c/sqrt(7),'lineProps',{'Color',[0.8 0 0.2]})

ax = gca; % current axes
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.98 1.02];
ax.XLim = [-20/4 290/4];
ax.XColor ='none';
ax.YColor ='k';
ax.FontSmoothing = 'on';
title(strcat('Mean linewidth of water time course - lb applied=',num2str(lb_spec_act)))

%Generating plot for amplitude
figure
hold on 
rectangle('Position',pos1,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos2,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos3,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos4,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
rectangle('Position',pos5,'Curvature', 0.2,'EdgeColor','none','FaceColor',[0.9 0.9 0.9])
plot(height_phasecycle_c,'Color','k','Linewidth',1.8)
shadedErrorBar([],height_phasecycle_c,std_height_phasecycle_c/sqrt(7),'lineProps','k')


ax = gca; % current axes
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.98 1.02];
ax.XLim = [-20/4 290/4];
ax.XColor ='none';
ax.YColor ='k';
ax.FontSmoothing = 'on';
title(strcat('Mean amplitude of water time course - lb applied=',num2str(lb_spec_act)))


%% Check independent time courses
linewidth=zeros(10,Nb_Repet);
height=zeros(10,Nb_Repet);
integration=zeros(10,Nb_Repet);
for i=1:10
    mat_file=WATER_TOT(:,(i-1)*Nb_Repet+1:(i-1)*Nb_Repet+Nb_Repet);
   for j=1:Nb_Repet
    zero_fill=abs(fftshift(fft([mat_file(:,j)' zeros(4096,1)']))); 
    
    [M I]=max(zero_fill);
    
    half_width=M/2;
    k=I;
    while zero_fill(1,k)>half_width
        k=k+1;
    end 
    
    k_2=I;
    while zero_fill(1,k_2)>half_width
        k_2=k_2-1;
    end
linewidth(i,j)=k-k_2;    
height(i,j)=M;
integration(i,j)=sum(zero_fill);
    end
end
% Mice 3 and 10 have signal intensity variation, and mouse 7's first point is too intense. 
mouse_water_functional=[1 2 4 5 6 8 9];
figure
plot(detrend(height(mouse_water_functional,:)',0)+1,'LineWidth',1.5)

