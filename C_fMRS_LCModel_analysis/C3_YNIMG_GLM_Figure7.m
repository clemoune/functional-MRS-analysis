%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs 2 types of GLMs on the mean time courses.
% 1) it generates the results for the Table 1 that are saved in the "group" 
% folders (C3 files).
% 2) it generates the metabolic response function with a finite impulse
% response approach (Figure 7)
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the mean time courses
currentdir=pwd;
load(strcat(currentdir,filesep,'functional',filesep,'C2_Mean_Norm_Timecourse.mat'))
load(strcat(currentdir,filesep,'control',filesep,'C2_Mean_Norm_Timecourse.mat'))

load(strcat(currentdir,filesep,'C1_Metabolites_List.mat'))
%% 1) Regressors 
time_axis=[0.1:0.1:13.5];

% lock to visual stimulation [Block Paradigm]
Reg_BP=[zeros(1,8) ones(1,4) zeros(1,8) ones(1,4) zeros(1,8) ones(1,4) zeros(1,8) ones(1,4) zeros(1,8) ones(1,4) zeros(1,8) zeros(1,67)];

% lock to "active" [Active/Recovery]
Reg_ActRec=[ones(1,68) zeros(1,67)];

% linear drift
Reg_Linear=flip([1:135]/135);

subplot(3,1,1)
plot(time_axis,Reg_BP, 'k')
axis([-0.5 14 -0.3 1.3])
title('Block Paradigm - stim events')

subplot(3,1,2)
plot(time_axis,Reg_ActRec, 'k')
axis([-0.5 14 -0.3 1.3])
title('Active - Recovery')

subplot(3,1,3)
plot(time_axis,Reg_Linear, 'k')
axis([-0.5 14 -0.3 1.3])
title('Linear Drift')



%% 1) Applying the glm to mean time courses - generating values in table 1 (without Bonferonni correction)

% Functional
for number_metab=[1:6]
[b, dev, stats]=glmfit([Reg_BP;Reg_ActRec;Reg_Linear]', Mean_TimeCourse_F_Norm(:,number_metab),'normal')
p_value_F(number_metab).name=metab(number_metab).name(1,:);
p_value_F(number_metab).Reg_BP=stats.p(2);
p_value_F(number_metab).Reg_ActRec=stats.p(3);
p_value_F(number_metab).Reg_Linear=stats.p(4);
end

% Control
for number_metab=[1:6]

[b, dev, stats]=glmfit([Reg_BP;Reg_ActRec;Reg_Linear]', Mean_TimeCourse_C_Norm(:,number_metab),'normal')
p_value_C(number_metab).name=metab(number_metab).name(1,:);
p_value_C(number_metab).Reg_BP=stats.p(2);
p_value_C(number_metab).Reg_ActRec=stats.p(3);
p_value_C(number_metab).Reg_Linear=stats.p(4);
end

save(strcat(currentdir,filesep,'functional',filesep,'C3_GLM_Table.mat'),'p_value_F','-mat')
save(strcat(currentdir,filesep,'control',filesep,'C3_GLM_Table.mat'),'p_value_C','-mat')

%% 2) Metabolic response function 

% Regressor
time_axis=[0.1:0.1:13.5];
Reg_FIR=zeros(13,135);
Reg_FIR(1,:)=[ones(1,10) zeros(1,125)];

plot(time_axis,Reg_FIR(1,:), 'k')
axis([-0.5 14 -1 18])
title('Finite impulse response')
hold on
for i=2:13
Reg_FIR(i,:)=[zeros(1,(i-1)*10) ones(1,10) zeros(1,135-i*10)];
plot(time_axis,Reg_FIR(i,:)+i*1.3-1, 'k')
end
hold off

% Metabolic Response Function
figure
for number_metab=[14 3 4 5]
[b, dev, stats]=glmfit(Reg_FIR', Mean_TimeCourse_F_Norm(:,number_metab),'normal')
f(number_metab)=plot(b(2:14),'LineWidth', 1.5, 'Color', colours(number_metab,:));
hold on
end
legend([f(14) f(3) f(4) f(5)],strcat(metab(14).name(1,:)), strcat(metab(3).name(1,:)), strcat(metab(4).name(1,:)), strcat(metab(5).name(1,:)))
hold off