%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes mean time courses and histograms found in Figure 6.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Loading time courses

currentdir=pwd;
load(strcat(currentdir,filesep,'functional',filesep,'C1_Concentrations.mat'))
load(strcat(currentdir,filesep,'control',filesep,'C1_Concentrations.mat'))
load(strcat(currentdir,filesep,'C1_Metabolites_List.mat'))

addpath(strcat(currentdir,filesep,'support_functions'))

mouse_number_functional=[1 3 5 6 8 9 10];
mouse_number_control=[1 2 3 4 6 7 9];

size_block=135;
%% Plotting the time courses

colours=zeros(22,3);
colours(1,:)=[0 0 1];
colours(2,:)=[0.8 0 0.8];
colours(3,:)=[0 0.9 0.0];
colours(4,:)=[1 0.2 0.1];
colours(5,:)=[0.95 0.75 0.02];
colours(6,:)=[0.35 0 0];
colours(8,:)=[0.4 0.4 0];
colours(12,:)=[1 0 1];
colours(13,:)=[0.5 0 0.5];
colours(14,:)=[0 0 1];
colours(15,:)=[0 0.5 0.7]; 
nb_indiv=7;
sizemat=67;

Mean_TimeCourse_F=zeros;
Mean_TimeCourse_C=zeros;
Mean_CRLB_F=zeros;
Mean_CRLB_C=zeros;

for Individu=mouse_number_functional
Mean_TimeCourse_F=Mean_TimeCourse_F+functional(Individu).table;
Mean_CRLB_F=Mean_CRLB_F+functional(Individu).crlb/nb_indiv;
end

for Individu=mouse_number_control
Mean_TimeCourse_C=Mean_TimeCourse_C+control(Individu).table;
Mean_CRLB_C=Mean_CRLB_C+control(Individu).crlb/nb_indiv;
end

Mean_TimeCourse_F_Norm=Mean_TimeCourse_F./mean(Mean_TimeCourse_F(round(size_block/2):size_block,:),1);
Mean_TimeCourse_C_Norm=Mean_TimeCourse_C./mean(Mean_TimeCourse_C(round(size_block/2):size_block,:),1);

CRLB_F=mean(Mean_CRLB_F,1);
CRLB_C=mean(Mean_CRLB_C,1);

for number_metab=[1:22]

        for it_std=1:size_block
            vec_func=[];
            vec_control=[];            
            for Individu=mouse_number_functional
            vec_func=[vec_func functional(Individu).table(it_std,number_metab)./functional(Individu).table(1,number_metab)];
            end
            for Individu=mouse_number_control
            vec_control=[vec_control control(Individu).table(it_std,number_metab)./control(Individu).table(1,number_metab)];
            end
    % Standard error         
    ser_control(it_std,number_metab)=std(vec_control)/sqrt(nb_indiv);
    ser_functional(it_std,number_metab)=std(vec_func)/sqrt(nb_indiv);
        end
        
A=ones(size_block,1);
Activation_onset=([32 80 128 176 224])/4;  
StimOn=NaN(size_block,1);
StimOn(Activation_onset)=1*max(Mean_TimeCourse_F_Norm(1:size_block,number_metab)+ ser_functional(1:size_block,number_metab));
StimOn(Activation_onset+1)=1*max(Mean_TimeCourse_F_Norm(1:size_block,number_metab)+ ser_functional(1:size_block,number_metab));
StimOn(Activation_onset+2)=1*max(Mean_TimeCourse_F_Norm(1:size_block,number_metab)+ ser_functional(1:size_block,number_metab));
StimOn(Activation_onset+3)=1*max(Mean_TimeCourse_F_Norm(1:size_block,number_metab)+ ser_functional(1:size_block,number_metab));
StimOn(Activation_onset+4)=1*max(Mean_TimeCourse_F_Norm(1:size_block,number_metab)+ ser_functional(1:size_block,number_metab));

time_axis=[0.1:0.1:13.5];

figure
shadedErrorBar(time_axis,Mean_TimeCourse_F_Norm(:,number_metab)',ser_functional(:,number_metab)','lineProps',{'Color',colours(number_metab,:)})
hold on 
plot(time_axis,Mean_TimeCourse_F_Norm(:,number_metab),'Color',colours(number_metab,:),'LineWidth',1.2)
plot(time_axis,A,'Color',[0 0 0],'LineStyle','--')
plot(time_axis,StimOn,'Color',[0 1 1],'LineWidth',5)
ax = gca; 
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.88 1.12];
ax.XLim = [-0.5 14];
ax.XColor ='k';
ax.YColor ='k';
ax.FontSmoothing = 'on';


title(strcat('Functional',{' '},metab(number_metab).name(1,:)))

figure
shadedErrorBar(time_axis,Mean_TimeCourse_C_Norm(:,number_metab)',ser_control(:,number_metab)','lineProps',{'Color',colours(number_metab,:)})
hold on 
plot(time_axis,Mean_TimeCourse_C_Norm(:,number_metab),'Color',colours(number_metab,:),'LineWidth',1.2)
plot(time_axis,A,'Color',[0 0 0],'LineStyle','--')
plot(time_axis,StimOn,'Color',[.5 .5 .5],'LineWidth',5)
ax = gca; 
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.88 1.12];
ax.XLim = [-0.5 14];
ax.XColor ='k';
ax.YColor ='k';
ax.FontSmoothing = 'on';


title(strcat('Control',{' '},metab(number_metab).name(1,:)))

end

save(strcat(currentdir,filesep,'control',filesep,'C2_Mean_Norm_Timecourse.mat'), 'Mean_TimeCourse_C_Norm', '-mat')
save(strcat(currentdir,filesep,'functional',filesep,'C2_Mean_Norm_Timecourse.mat'), 'Mean_TimeCourse_F_Norm', '-mat')
save(strcat(currentdir,filesep,'C2_Colour_Codes.mat'), 'colours', '-mat')
%% Plotting histograms


for number_metab=[1:6]
    functional_pool_act=[];
    functional_pool_rest=[];
    control_pool_act=[];
    control_pool_rest=[];    

    for t=1:sizemat

     for Individu=mouse_number_functional  
        functional_pool_act= [functional_pool_act functional(Individu).table(t,number_metab)./mean(functional(Individu).table(sizemat+1:end,number_metab))];
        functional_pool_rest=[functional_pool_rest functional(Individu).table(sizemat+1+t,number_metab)./mean(functional(Individu).table(sizemat+1:end,number_metab))];
     end

          for Individu=mouse_number_control  
        control_pool_act= [control_pool_act control(Individu).table(t,number_metab)./mean(control(Individu).table(sizemat+1:end,number_metab))];
        control_pool_rest=[control_pool_rest control(Individu).table(sizemat+1+t,number_metab)./mean(control(Individu).table(sizemat+1:end,number_metab))];
     end
end


max_b=max([functional_pool_act functional_pool_rest]);
min_b=min([functional_pool_act functional_pool_rest]);
bin_size = 0.01;


figure
bin_vec = min_b:bin_size:max_b;
h1=histogram(functional_pool_rest,bin_vec,'normalization','probability', 'FaceAlpha', 0.8, 'FaceColor', [ 0.5  0.5 0.5]); hold on
h2=histogram(functional_pool_act,bin_vec,'normalization','probability', 'FaceColor', colours(number_metab,:)); hold off
ylabel('Counts (/total)'); xlabel('Normalized NAA concentration'); axis([0.85 1.15 0 0.21])
legend([h1,h2],'Recovery','BP')
title(strcat('Functional',{' '},metab(number_metab).name(1,:)))

max_b=max([control_pool_act control_pool_rest]);
min_b=min([control_pool_act control_pool_rest]);
bin_size = 0.01;

figure
bin_vec = min_b:bin_size:max_b;
h1=histogram(control_pool_rest,bin_vec,'normalization','probability', 'FaceAlpha', 0.8, 'FaceColor', [ 0.5  0.5 0.5]); hold on
h2=histogram(control_pool_act,bin_vec,'normalization','probability', 'FaceColor', colours(number_metab,:)); hold off
ylabel('Counts (/total)'); xlabel('Normalized concentration'); axis([0.85 1.15 0 0.21])
legend([h1,h2],'No Stim 2','No Stim 1')
title(strcat('Control',{' '},metab(number_metab).name(1,:)))

end


