%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script detrends the time courses.
% It needs the normalised time courses found in the A2 files in 
% "group/metabolites/matlab_norm_timecourses"
% and it generates the total detrended time courses in the 
% "group/metabolites/matlab_detrended_timecourse" folder (A3 files).
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading data 
currentdir=pwd;

%Functional
group='functional';
directory_norm=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_norm_timecourse',filesep)
load((strcat(directory_norm,'A2_normalised_timecouse.mat')))

%Control
group='control';
directory_norm=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_norm_timecourse',filesep)
load((strcat(directory_norm,'A2_normalised_timecouse.mat')))

%% Transforming back to FID

for mouse_number=1:10
func_res(mouse_number).tttc=ifft(functional_norm(mouse_number).tttc);
control_res(mouse_number).tttc=ifft(control_norm(mouse_number).tttc);
end


%% Calculating time courses mean linewidth & height based on the 3 singlets
% note that the course_lw function is calibrated for these acquisitions (ranges for singlets are hardly coded)
size_timecourse=679;
addpath(strcat(currentdir,filesep,'support_functions'))


for mouse_number = 1:10

[linewidth_f,height_f] = course_lw(func_res(mouse_number).tttc,size_timecourse);
[linewidth_c,height_c] = course_lw(control_res(mouse_number).tttc,size_timecourse);

linewidth_metab_c(mouse_number,:)=linewidth_c;
height_metab_c(mouse_number,:)=height_c;
linewidth_metab_f(mouse_number,:)=linewidth_f;
height_metab_f(mouse_number,:)=height_f;
end

%% CONTROL: checking mean drift

time_vector = [0.1000:0.1000:67.9];

figure
hold on 
plot(time_vector,mean(linewidth_metab_c([1:10],:),1)/mean(mean(linewidth_metab_c([1:10],1:size_timecourse),1))-0.4,'Color',[0/256 128/256 0/256],'Linewidth',1.5)
plot(time_vector,mean(height_metab_c([1:10],:),1)/mean(mean(height_metab_c([1:10],1:size_timecourse),1)),'k','LineWidth',2)
A=ones(size_timecourse,1);
A_2=ones(size_timecourse,1)*0.6;
plot(time_vector,A,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)
plot(time_vector,A_2,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)

ax = gca; % current axes
ax.FontSize = 14;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.2 1.2];
ax.XLim = [-0.15 68];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title('Mean height (black) and linewidth (green) - control')

% check the mean ratio between the first and last block
control_lw_ratio=mean(mean(linewidth_metab_c([1:10],1:136),1)/mean(mean(linewidth_metab_c([1:10],4*136:679),1)))
control_ampl_ratio=mean(mean(height_metab_c([1:10],1:136),1)/mean(mean(height_metab_c([1:10],4*136:679),1)))

%% FUNCTIONAL: checking mean drift

time_vector = [0.1000:0.1000:67.9];

figure
hold on 
plot(time_vector,mean(linewidth_metab_f([1:10],:),1)/mean(mean(linewidth_metab_f([1:10],1:size_timecourse),1))-0.4,'Color',[0/256 128/256 0/256],'Linewidth',1.5)
plot(time_vector,mean(height_metab_f([1:10],:),1)/mean(mean(height_metab_f([1:10],1:size_timecourse),1)),'k','LineWidth',2)
A=ones(679,1);
A_2=ones(679,1)*0.6;
plot(time_vector,A,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)
plot(time_vector,A_2,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)

ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.2 1.2];
ax.XLim = [-0.15 68];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title('Mean height (black) and linewidth (green) - functional')

% check the mean ratio between the first and last block
func_lw_ratio=mean(mean(linewidth_metab_f([1:10],1:136),1)/mean(mean(linewidth_metab_f([1:10],4*136:679),1)))
func_ampl_ratio=mean(mean(height_metab_f([1:10],1:136),1)/mean(mean(height_metab_f([1:10],4*136:679),1)))


%% CONTROL: Fit the slope for total time course linewidth correction
% figures for individual visual check (blue > linewidth)

control_slope_height=[];
control_intercept_height=[];
control_slope_lw=[];
control_intercept_lw=[];

for mouse_number=[1:10]

figure
f1 = fittype('a*x+b');
[fit1,gof1,fitinfo1] = fit([1:size_timecourse]',(mean(height_metab_c(mouse_number,:),1))',f1,'StartPoint',[1 1]);
plot(fit1.a*[1:size_timecourse]+fit1.b)
hold on 
plot(mean(height_metab_c(mouse_number,:),1))

control_slope_height=[control_slope_height fit1.a];
control_intercept_height=[control_intercept_height fit1.b];

figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_metab_c(mouse_number,:),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_metab_c(mouse_number,:),1),'b')
control_slope_lw=[control_slope_lw fit2.a];
control_intercept_lw=[control_intercept_lw fit2.b];


end

%% FUNCTIONAL: Fit the slope for total time course linewidth correction
% figures for individual visual check (blue > linewidth)

func_slope_height=[];
func_intercept_height=[];
func_slope_lw=[];
func_intercept_lw=[];


for mouse_number=[1:10]

figure
f1 = fittype('a*x+b');
[fit1,gof1,fitinfo1] = fit([1:size_timecourse]',(mean(height_metab_f(mouse_number,:),1))',f1,'StartPoint',[1 1]);
plot(fit1.a*[1:size_timecourse]+fit1.b)
hold on 
plot(mean(height_metab_f(mouse_number,:),1))

func_slope_height=[func_slope_height fit1.a];
func_intercept_height=[func_intercept_height fit1.b];

figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_metab_f(mouse_number,:),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_metab_f(mouse_number,:),1),'b')
func_slope_lw=[func_slope_lw fit2.a];
func_intercept_lw=[func_intercept_lw fit2.b];


end
%% Block-outlier-linewidth detection & visualisation

for mouse_number=1:10
    control_blocks_lw(1:5,mouse_number)=[mean(linewidth_metab_c(mouse_number,1:136),2) mean(linewidth_metab_c(mouse_number,137:136+136),2) mean(linewidth_metab_c(mouse_number,1+2*136:136+2*136),2) mean(linewidth_metab_c(mouse_number,1+3*136:136+3*136),2) mean(linewidth_metab_c(mouse_number,1+4*136:135+4*136),2)];
    control_blocks_height(1:5,mouse_number)=[mean(height_metab_c(mouse_number,1:136),2) mean(height_metab_c(mouse_number,137:136+136),2) mean(height_metab_c(mouse_number,1+2*136:136+2*136),2) mean(height_metab_c(mouse_number,1+3*136:136+3*136),2) mean(height_metab_c(mouse_number,1+4*136:135+4*136),2)];
   
    functional_blocks_lw(1:5,mouse_number)=[mean(linewidth_metab_f(mouse_number,1:136),2) mean(linewidth_metab_f(mouse_number,137:136+136),2) mean(linewidth_metab_f(mouse_number,1+2*136:136+2*136),2) mean(linewidth_metab_f(mouse_number,1+3*136:136+3*136),2) mean(linewidth_metab_f(mouse_number,1+4*136:135+4*136),2)];
    functional_blocks_height(1:5,mouse_number)=[mean(height_metab_f(mouse_number,1:136),2) mean(height_metab_f(mouse_number,137:136+136),2) mean(height_metab_f(mouse_number,1+2*136:136+2*136),2) mean(height_metab_f(mouse_number,1+3*136:136+3*136),2) mean(height_metab_f(mouse_number,1+4*136:135+4*136),2)];
  
    figure 
    plot(control_blocks_lw(:,mouse_number))
    axis([1 5 20 30])
    title(strcat('Control - blocks mean linewidth - mouse',num2str(mouse_number)))
    figure 
    plot(functional_blocks_lw(:,mouse_number))
    axis([1 5 20 30])
    title(strcat('Functional - blocks mean linewidth - mouse',num2str(mouse_number)))  
end

%% CONTROL: correction for drifts and outliers 
% with manual adjustments of slopes on blocks to get flat linewidth timecourse

bw=4000;
time=[1:2048];
%%%%%%%%%%%%1
        control_res_corr1(1).tttc=control_res(1).tttc;
%outlier        
        control_res_corr1(1).tttc(:,136)=control_res(1).tttc(:,135);
%drift        
        if control_slope_lw(1)>0
            for k=1:size_timecourse
             control_res_corr2(1).tttc(:,k)=(control_res_corr1(1).tttc(:,k)).*exp(-(size_timecourse-k)*control_slope_lw(1)/1*time'/bw);
            end  
        else
            for k=1:size_timecourse
             control_res_corr2(1).tttc(:,k)=(control_res_corr1(1).tttc(:,k)).*exp(-(-k)*control_slope_lw(1)/1*time'/bw);
            end
        end

        
%%%%%%%%%%%%2   

        control_res_corr1(2).tttc=control_res(2).tttc;
%block  4      
for i=1+3*136:4*136        
control_res_corr1(2).tttc(:,i)=control_res(2).tttc(:,i).*exp(-abs(control_blocks_lw(4,2)-control_blocks_lw(5,2))/1.5*time'/bw); %manually increase
end

[linewidth_c_2,height_c_2] = course_lw(control_res_corr1(2).tttc,size_timecourse);

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_2(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_2(1,:),1))
slope_lw_2=fit2.a;

        if slope_lw_2>0
            for k=1:size_timecourse
             control_res_corr2(2).tttc(:,k)=(control_res_corr1(2).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_2/2*time'/bw);
            end  
        else
            for k=1:size_timecourse
             control_res_corr2(2).tttc(:,k)=(control_res_corr1(2).tttc(:,k)).*exp(-(-k)*slope_lw_2/2*time'/bw);
            end
        end    

        
%%%%%%%%%%%%3
        control_res_corr1(3).tttc=control_res(3).tttc;
%outlier        
        control_res_corr1(3).tttc(:,136)=control_res(3).tttc(:,135);
        
%block  1      
for i=1:135        
control_res_corr1(3).tttc(:,i)=control_res(3).tttc(:,i).*exp(-abs(control_blocks_lw(1,3)-control_blocks_lw(2,3))/1.5*time'/bw); %manually increase
end
        
        [linewidth_c_3,height_c_3] = course_lw(control_res_corr1(3).tttc,544);

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:544]',(mean(linewidth_c_3(1,1:544),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:544]+fit2.b)
hold on 
plot(mean(linewidth_c_3(1,:),1))
slope_lw_3=fit2.a;

        if slope_lw_3>0
            for k=1:544
             control_res_corr2(3).tttc(:,k)=(control_res_corr1(3).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_3/1*time'/bw);
            end  
        else
            for k=1:544
             control_res_corr2(3).tttc(:,k)=(control_res_corr1(3).tttc(:,k)).*exp(-(-k)*slope_lw_3/1*time'/bw);
            end
        end    
 control_res_corr2(3).tttc(:,545:679)=control_res_corr1(3).tttc(:,545:679);

 %%%%%%%%%%%%4   

        control_res_corr1(4).tttc=control_res(4).tttc;
%block  4      
for i=1+3*136:4*136        
control_res_corr1(4).tttc(:,i)=control_res(4).tttc(:,i).*exp(-abs(control_blocks_lw(4,4)-control_blocks_lw(5,4))/1.5*time'/bw);
end

[linewidth_c_4,height_c_4] = course_lw(control_res_corr1(4).tttc,size_timecourse)

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_4(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:679]+fit2.b)
hold on 
plot(mean(linewidth_c_4(1,:),1))
slope_lw_4=fit2.a;

        if slope_lw_4>0
            for k=1:size_timecourse
             control_res_corr2(4).tttc(:,k)=(control_res_corr1(4).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_4/1*time'/bw);
            end  
        else
            for k=1:679
             control_res_corr2(4).tttc(:,k)=(control_res_corr1(4).tttc(:,k)).*exp(-(-k)*slope_lw_4/1*time'/bw);
            end
        end      

        
%%%%%%%%%%%%5   

        control_res_corr1(5).tttc=control_res(5).tttc;

%drift        
        if control_slope_lw(5)>0
            for k=1:size_timecourse
             control_res_corr2(5).tttc(:,k)=(control_res_corr1(5).tttc(:,k)).*exp(-(size_timecourse-k)*control_slope_lw(5)/1*time'/bw);
            end  
        else
            for k=1:size_timecourse
             control_res_corr2(5).tttc(:,k)=(control_res_corr1(5).tttc(:,k)).*exp(-(-k)*control_slope_lw(5)/1*time'/bw);
            end
        end        
        
        
%%%%%%%%%%%%6   

        control_res_corr1(6).tttc=control_res(6).tttc;
%general outlier
%block  2-3-4      
for i=1:4*136        
control_res_corr1(6).tttc(:,i)=control_res(6).tttc(:,i).*exp(-abs(control_blocks_lw(1,6)-control_blocks_lw(4,6))/2.5*time'/bw);
end

[linewidth_c_6,height_c_6] = course_lw(control_res_corr1(6).tttc,size_timecourse)

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_6(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_6(1,:),1))
slope_lw_6=fit2.a;

        if slope_lw_6>0
            for k=1:size_timecourse
             control_res_corr2(6).tttc(:,k)=(control_res_corr1(6).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_6*time'/bw);
            end  
        else
            for k=1:679
             control_res_corr2(6).tttc(:,k)=(control_res_corr1(6).tttc(:,k)).*exp(-(-k)*slope_lw_6/2*time'/bw);
            end
        end      

        
%%%%%%%%%%%%7
        control_res_corr1(7).tttc=control_res(7).tttc;
%outlier        
        control_res_corr1(7).tttc(:,136)=control_res(7).tttc(:,135);

%block 1
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:136]',(mean(linewidth_metab_c(7,1:136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:136]+fit2.b)
hold on 
plot(mean(linewidth_metab_c(7,:),1))
slope_lw_7a=fit2.a;

        if slope_lw_7a>0
            for k=1:136
             control_res_corr1(7).tttc(:,k)=(control_res(7).tttc(:,k)).*exp(-(136-k)*slope_lw_7a/1*time'/bw);
            end  
        else
            for k=1:136
             control_res_corr1(7).tttc(:,k)=(control_res(7).tttc(:,k)).*exp(-(-k)*slope_lw_7a/1*time'/bw);
            end
        end 
        
[linewidth_c_7,height_c_7] = course_lw(control_res_corr1(7).tttc,size_timecourse);
control_res_corr2(7).tttc=control_res_corr1(7).tttc;
%drift 2 3 4 5       
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([137:679]',(mean(linewidth_c_7(1,137:679),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[137:679]+fit2.b)
hold on 
plot(mean(linewidth_c_7(1,:),1),'k')
slope_lw_7=fit2.a;

        if slope_lw_7>0
            for k=137:size_timecourse
             control_res_corr2(7).tttc(:,k)=(control_res_corr1(7).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_7/1*time'/bw);
            end  
        else
            for k=137:size_timecourse
             control_res_corr2(7).tttc(:,k)=(control_res_corr1(7).tttc(:,k)).*exp(-(-k)*slope_lw_7/1*time'/bw);
            end
        end        

       
%%%%%%%%%%%%8   

        control_res_corr1(8).tttc=control_res(8).tttc;

        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1+3*136:4*136]',(mean(linewidth_metab_c(8,1+3*136:4*136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1+3*136:4*136]+fit2.b)
hold on 
plot(mean(linewidth_metab_c(8,:),1))
slope_lw_8a=fit2.a;

        if slope_lw_8a>0
            for k=1+3*136:4*136
             control_res_corr1(8).tttc(:,k)=(control_res(8).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_8a/2*time'/bw);
            end  
        else
            for k=1+3*136:4*136
             control_res_corr1(8).tttc(:,k)=(control_res(8).tttc(:,k)).*exp(-(-k)*slope_lw_8a/2*time'/bw);
            end
        end   
%outlier
control_res_corr1(8).tttc(:,544)=control_res_corr1(8).tttc(:,543);

[linewidth_c_8,height_c_8] = course_lw(control_res_corr1(8).tttc,size_timecourse)

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_8(1,1:679),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_8(1,:),1))
slope_lw_8=fit2.a;

        if slope_lw_8>0
            for k=1:size_timecourse
             control_res_corr2(8).tttc(:,k)=(control_res_corr1(8).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_8/2*time'/bw);
            end  
        else
            for k=1:size_timecourse
             control_res_corr2(8).tttc(:,k)=(control_res_corr1(8).tttc(:,k)).*exp(-(-k)*slope_lw_8/2*time'/bw);
            end
        end              

%%%%%%%%%%%%9   

        control_res_corr1(9).tttc=control_res(9).tttc;
%drift        
        if control_slope_lw(9)>0
            for k=1:size_timecourse
             control_res_corr2(9).tttc(:,k)=(control_res_corr1(9).tttc(:,k)).*exp(-(size_timecourse-k)*control_slope_lw(9)/1.5*time'/bw);
            end  
        else
            for k=1:size_timecourse
             control_res_corr2(9).tttc(:,k)=(control_res_corr1(9).tttc(:,k)).*exp(-(-k)*control_slope_lw(9)/1.5*time'/bw);
            end
        end             

        
%%%%%%%%%%%%10   

        control_res_corr1(10).tttc=control_res(10).tttc;
        

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:136]',(mean(linewidth_metab_c(10,1:136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:136]+fit2.b)
hold on 
plot(mean(linewidth_metab_c(10,:),1))
slope_lw_10a=fit2.a;

f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1+136:136*2]',(mean(linewidth_metab_c(10,1+136:136*2),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1+136:136*2]+fit2.b)
hold on 
plot(mean(linewidth_metab_c(10,:),1))
slope_lw_10b=fit2.a;


        if slope_lw_10a>0
            for k=1:136
             control_res_corr2(10).tttc(:,k)=(control_res_corr1(10).tttc(:,k)).*exp(-(136-k)*slope_lw_10a/2*time'/bw);
            end  
        else
            for k=1:136
             control_res_corr2(10).tttc(:,k)=(control_res_corr1(10).tttc(:,k)).*exp(-(-k)*slope_lw_10a/2*time'/bw);
            end
        end      
        
        if slope_lw_10b>0
            for k=1+136:136*2
             control_res_corr2(10).tttc(:,k)=(control_res_corr1(10).tttc(:,k)).*exp(-(136-k)*slope_lw_10b/2*time'/bw);
            end  
        else
            for k=1+136:136*2
             control_res_corr2(10).tttc(:,k)=(control_res_corr1(10).tttc(:,k)).*exp(-(-k)*slope_lw_10b/2*time'/bw);
            end
        end      
control_res_corr2(10).tttc(:,136*2+1:size_timecourse)=(control_res_corr1(10).tttc(:,136*2+1:size_timecourse));              
        
%block 2-5      
for i=1+136:size_timecourse        
control_res_corr2(10).tttc(:,i)=control_res_corr2(10).tttc(:,i).*exp(-abs(control_blocks_lw(1,10)-control_blocks_lw(5,10))*2*time'/bw);
end

    
%% correction functional 

%%%%%%%%%%%%1
        func_res_corr1(1).tttc=func_res(1).tttc;
%drift        
        if func_slope_lw(1)>0
            for k=1:size_timecourse
             func_res_corr2(1).tttc(:,k)=(func_res_corr1(1).tttc(:,k)).*exp(-(size_timecourse-k)*func_slope_lw(1)/3*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(1).tttc(:,k)=(func_res_corr1(1).tttc(:,k)).*exp(-(-k)*func_slope_lw(1)/3*time'/bw);
            end
        end

        
%%%%%%%%%%%%2   
        func_res_corr1(2).tttc=func_res(2).tttc;

%drift 1
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:136]',(mean(linewidth_metab_f([2],1:136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([2],:),1))
slope_lw_2b=fit2.a;

        if slope_lw_2b>0
            for k=1:136
             func_res_corr1(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(130-k)*slope_lw_2b/1.5*time'/bw);
            end  
        else
            for k=1:136
             func_res_corr1(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(-k)*slope_lw_2b/1.5*time'/bw);
            end
        end   


%drift 2
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([137:272]',(mean(linewidth_metab_f([2],137:272),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[137:272]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([2],:),1))
slope_lw_2a=fit2.a;

        if slope_lw_2a>0
            for k=137:272
             func_res_corr1(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(272-k)*slope_lw_2a/1.0*time'/bw);
            end  
        else
            for k=137:272
             func_res_corr1(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(-k)*slope_lw_2a/1.0*time'/bw);
            end
        end   

        [linewidth_c_2,height_c_2] = course_lw(func_res_corr1(2).tttc,size_timecourse);
        
% total drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_2(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_2(1,:),1))
slope_lw_2=fit2.a;

        if slope_lw_2>0
            for k=1:size_timecourse
             func_res_corr2(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_2/3*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(2).tttc(:,k)=(func_res_corr1(2).tttc(:,k)).*exp(-(-k)*slope_lw_2/3*time'/bw);
            end
        end
        
%%%%%%%%%%%%3         
        func_res_corr1(3).tttc=func_res(3).tttc;
%outlier        
        func_res_corr1(3).tttc(:,136)=func_res(3).tttc(:,135);
%drift only 2 first blocks

figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:300]',(mean(linewidth_metab_f([3],1:300),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:300]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([3],:),1))
slope_lw_3=fit2.a;

        if slope_lw_3>0
            for k=1:300
             func_res_corr2(3).tttc(:,k)=(func_res_corr1(3).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_3/3*time'/bw);
            end  
        else
            for k=1:300
             func_res_corr2(3).tttc(:,k)=(func_res_corr1(3).tttc(:,k)).*exp(-(-k)*slope_lw_3/3*time'/bw);
            end
        end    
        
 func_res_corr2(3).tttc(:,301:size_timecourse)=func_res_corr1(3).tttc(:,301:size_timecourse);
        
%%%%%%%%%%%%4         
        func_res_corr1(4).tttc=func_res(4).tttc;
%block replacement 3 for 1        
%         func_res_corr1(4).tttc(:,1:136)=func_res(4).tttc(:,1+2*136:3*136);


%drift 1
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:136]',(mean(linewidth_metab_f([4],1:136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([4],:),1))
slope_lw_4b=fit2.a;

        if slope_lw_4b>0
            for k=1:136
             func_res_corr1(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(136-k)*slope_lw_4b/1.5*time'/bw);
            end  
        else
            for k=1:136
             func_res_corr1(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(-k)*slope_lw_4b/1.5*time'/bw);
            end
        end   


%drift 2-3
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([136:408]',(mean(linewidth_metab_f([4],136:408),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[136:408]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([4],:),1))
slope_lw_4a=fit2.a;

        if slope_lw_4a>0
            for k=136:408
             func_res_corr1(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(408-k)*slope_lw_4a/1.5*time'/bw);
            end  
        else
            for k=136:408
             func_res_corr1(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(-k)*slope_lw_4a/1.5*time'/bw);
            end
        end   

        [linewidth_c_4,height_c_4] = course_lw(func_res_corr1(4).tttc,size_timecourse);
        
% total drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_4(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_4(1,:),1))
slope_lw_4=fit2.a;

        if slope_lw_4>0
            for k=1:size_timecourse
             func_res_corr2(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_4/3*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(4).tttc(:,k)=(func_res_corr1(4).tttc(:,k)).*exp(-(-k)*slope_lw_4/3*time'/bw);
            end
        end
           
%%%%%%%%%%%%5         
        func_res_corr1(5).tttc=func_res(5).tttc;
%drift 1 / 234 /5
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:136]',(mean(linewidth_metab_f([5],1:136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([5],:),1))
slope_lw_5a=fit2.a;

[fit2,gof2,fitinfo2] = fit([1+136:4*136]',(mean(linewidth_metab_f([5],1+136:4*136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1+136:4*136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([3],:),1))
slope_lw_5b=fit2.a;

[fit2,gof2,fitinfo2] = fit([1+4*136:135+4*136]',(mean(linewidth_metab_f([5],1+4*136:135+4*136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1+4*136:135+4*136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([3],:),1))
slope_lw_5c=fit2.a;

        if slope_lw_5a>0
            for k=1:136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(135-k)*slope_lw_5a/3*time'/bw);
            end  
        else
            for k=1:136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(-k)*slope_lw_5a/3*time'/bw);
            end
        end       
        
        if slope_lw_5a>0
            for k=1+136:4*136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(4*136-k)*slope_lw_5a/3*time'/bw);
            end  
        else
            for k=1+136:4*136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(-k)*slope_lw_5a/3*time'/bw);
            end
        end          
        
        
        if slope_lw_5a>0
            for k=1+4*136:135+4*136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(135+4*136-k)*slope_lw_5a/3*time'/bw);
            end  
        else
            for k=1+4*136:135+4*136
             func_res_corr2(5).tttc(:,k)=(func_res_corr1(5).tttc(:,k)).*exp(-(-k)*slope_lw_5a/3*time'/bw);
            end
        end  
        
%%%%%%%%%%%%6   
        func_res_corr1(6).tttc=func_res(6).tttc;
%drift        
        if func_slope_lw(6)>0
            for k=1:size_timecourse
             func_res_corr2(6).tttc(:,k)=(func_res_corr1(6).tttc(:,k)).*exp(-(size_timecourse-k)*func_slope_lw(6)/3*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(6).tttc(:,k)=(func_res_corr1(6).tttc(:,k)).*exp(-(-k)*func_slope_lw(6)/3*time'/bw);
            end
        end        
        
%%%%%%%%%%%%7         
        func_res_corr1(7).tttc=func_res(7).tttc;
%outlier        
        func_res_corr1(7).tttc(:,136)=func_res(7).tttc(:,135);        
        func_res_corr1(7).tttc(:,272)=func_res(7).tttc(:,271);   
        func_res_corr1(7).tttc(:,214)=func_res(7).tttc(:,213); 
        func_res_corr1(7).tttc(:,502)=func_res(7).tttc(:,501); 
        func_res_corr1(7).tttc(:,507)=func_res(7).tttc(:,506); 
%drift 1-2
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:272]',(mean(linewidth_metab_f([7],1:272),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:272]+fit2.b)
hold on 
plot(mean(linewidth_metab_f([7],:),1))
slope_lw_7a=fit2.a;


        if slope_lw_7a>0
            for k=1:272
             func_res_corr1(7).tttc(:,k)=(func_res_corr1(7).tttc(:,k)).*exp(-(272-k)*slope_lw_7a/1.0*time'/bw);
            end  
        else
            for k=1:272
             func_res_corr1(7).tttc(:,k)=(func_res_corr1(7).tttc(:,k)).*exp(-(-k)*slope_lw_7a/1.0*time'/bw);
            end
        end       
        
[linewidth_c_7,height_c_7] = course_lw(func_res_corr1(7).tttc,size_timecourse);

% total drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_7(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_7(1,:),1))
slope_lw_7=fit2.a;

        if slope_lw_7>0
            for k=1:size_timecourse
             func_res_corr2(7).tttc(:,k)=(func_res_corr1(7).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_7/2*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(7).tttc(:,k)=(func_res_corr1(7).tttc(:,k)).*exp(-(-k)*slope_lw_7/2*time'/bw);
            end
        end                  
            
%%%%%%%%%%%%8         
        func_res_corr1(8).tttc=func_res(8).tttc;
%outlier        
        func_res_corr1(8).tttc(:,443)=func_res(8).tttc(:,442);
%drift        
        if func_slope_lw(8)>0
            for k=1:size_timecourse
             func_res_corr2(8).tttc(:,k)=(func_res_corr1(8).tttc(:,k)).*exp(-(size_timecourse-k)*func_slope_lw(8)/3*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(8).tttc(:,k)=(func_res_corr1(8).tttc(:,k)).*exp(-(-k)*func_slope_lw(8)/3*time'/bw);
            end
        end    
        
%%%%%%%%%%%%9         
        func_res_corr1(9).tttc=func_res(9).tttc;
%block 1     
for i=1:136       
func_res_corr1(9).tttc(:,i)=func_res(9).tttc(:,i).*exp(-abs(functional_blocks_lw(1,9)-functional_blocks_lw(3,9))/1.5*time'/bw);
end

[linewidth_c_9,height_c_9] = course_lw(func_res_corr1(9).tttc,size_timecourse);

%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_9(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_9(1,:),1))
slope_lw_9=fit2.a;

        if slope_lw_9>0
            for k=1:size_timecourse
             func_res_corr2(9).tttc(:,k)=(func_res_corr1(9).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_9/2*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(9).tttc(:,k)=(func_res_corr1(9).tttc(:,k)).*exp(-(-k)*slope_lw_9/2*time'/bw);
            end
        end          

        
%%%%%%%%%%%%10        
        func_res_corr1(10).tttc=func_res(10).tttc;
%outlier        
        func_res_corr1(10).tttc(:,136)=func_res(10).tttc(:,135);

%drift_block3
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1+2*136:3*136]',(mean(linewidth_metab_f(10,1+2*136:3*136),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1+2*136:3*136]+fit2.b)
hold on 
plot(mean(linewidth_metab_f(10,:),1))
slope_lw_10a=fit2.a;

        if slope_lw_10a>0
            for k=1:1:1+2*136:3*136
             func_res_corr1(10).tttc(:,k)=(func_res_corr1(10).tttc(:,k)).*exp(-(3*136-k)*slope_lw_10a*time'/bw);
            end  
        else
            for k=1:1:1+2*136:3*136
             func_res_corr1(10).tttc(:,k)=(func_res_corr1(10).tttc(:,k)).*exp(-(-k)*slope_lw_10a*time'/bw);
            end
        end              
[linewidth_c_10,height_c_10] = course_lw(func_res_corr1(10).tttc,size_timecourse);        
%drift        
figure
f2 = fittype('a*x+b');
[fit2,gof2,fitinfo2] = fit([1:size_timecourse]',(mean(linewidth_c_10(1,1:size_timecourse),1))',f2,'StartPoint',[1 1]);
plot(fit2.a*[1:size_timecourse]+fit2.b)
hold on 
plot(mean(linewidth_c_10(1,:),1))
slope_lw_10=fit2.a;

        if slope_lw_10>0
            for k=1:size_timecourse
             func_res_corr2(10).tttc(:,k)=(func_res_corr1(10).tttc(:,k)).*exp(-(size_timecourse-k)*slope_lw_10/2*time'/bw);
            end  
        else
            for k=1:size_timecourse
             func_res_corr2(10).tttc(:,k)=(func_res_corr1(10).tttc(:,k)).*exp(-(-k)*slope_lw_10/2*time'/bw);
            end
        end          
      
        
        
%% Checking the linear linewidth correction

size_timecourse=679;
addpath(strcat(currentdir,filesep,'support_functions'))


for mouse_number = 1:10

[linewidth_f,height_f] = course_lw(func_res_corr2(mouse_number).tttc,size_timecourse);
[linewidth_c,height_c] = course_lw(control_res_corr2(mouse_number).tttc,size_timecourse);

linewidth_metab_c_corr(mouse_number,:)=linewidth_c;
height_metab_c_corr(mouse_number,:)=height_c;
linewidth_metab_f_corr(mouse_number,:)=linewidth_f;
height_metab_f_corr(mouse_number,:)=height_f;
end

%% CONTROL: checking the linear linewidth correction

time_vector = [0.1000:0.1000:67.9];

figure
hold on 
plot(time_vector,mean(linewidth_metab_c_corr([1:10],:),1)/mean(mean(linewidth_metab_c_corr([1:10],1:size_timecourse),1))-0.4,'Color',[0/256 128/256 0/256],'Linewidth',1.5)
plot(time_vector,mean(height_metab_c_corr([1:10],:),1)/mean(mean(height_metab_c_corr([1:10],1:size_timecourse),1)),'k','LineWidth',2)
A=ones(size_timecourse,1);
A_2=ones(size_timecourse,1)*0.6;
plot(time_vector,A,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)
plot(time_vector,A_2,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)

ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.2 1.2];
ax.XLim = [-0.15 68];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title('Mean height (black) and linewidth (green) - control - corrected')

% check the mean ratio between the first and last block
control_lw_ratio_corr=mean(mean(linewidth_metab_c_corr([1:10],1:136),1)/mean(mean(linewidth_metab_c_corr([1:10],4*136:679),1)))
control_ampl_ratio_corr=mean(mean(height_metab_c_corr([1:10],1:136),1)/mean(mean(height_metab_c_corr([1:10],4*136:679),1)))

%% FUNCTIONAL: checking mean drift

time_vector = [0.1000:0.1000:67.9];

figure
hold on 
plot(time_vector,mean(linewidth_metab_f_corr([1:10],:),1)/mean(mean(linewidth_metab_f_corr([1:10],1:size_timecourse),1))-0.4,'Color',[0/256 128/256 0/256],'Linewidth',1.5)
plot(time_vector,mean(height_metab_f_corr([1:10],:),1)/mean(mean(height_metab_f_corr([1:10],1:size_timecourse),1)),'k','LineWidth',2)
A=ones(679,1);
A_2=ones(679,1)*0.6;
plot(time_vector,A,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)
plot(time_vector,A_2,'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',2)

ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.YLim = [0.2 1.2];
ax.XLim = [-0.15 68];
ax.XColor ='k';
ax.YColor ='none';
ax.FontSmoothing = 'on';
ax.LineWidth = 1.2;
title('Mean height (black) and linewidth (green) - functional - corrected')

% check the mean ratio between the first and last block
func_lw_ratio_corr=mean(mean(linewidth_metab_f_corr([1:10],1:136),1)/mean(mean(linewidth_metab_f_corr([1:10],4*136:679),1)))
func_ampl_ratio_corr=mean(mean(height_metab_f_corr([1:10],1:136),1)/mean(mean(height_metab_f_corr([1:10],4*136:679),1)))

%% save timecourses

%Functional
group='functional';
directory_detrended=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
mkdir(directory_detrended)
save(strcat(directory_detrended,'A3_detrended_timecouse.mat'),'func_res_corr2','linewidth_metab_f','height_metab_f','linewidth_metab_f_corr','height_metab_f_corr')

%Control
group='control';
directory_detrended=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'matlab_detrended_timecourse',filesep)
mkdir(directory_detrended)
save(strcat(directory_detrended,'A3_detrended_timecouse.mat'),'control_res_corr2','linewidth_metab_c','height_metab_c','linewidth_metab_c_corr','height_metab_c_corr')


