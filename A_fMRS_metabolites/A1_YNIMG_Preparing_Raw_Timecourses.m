%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script prepares the time courses, phase&frequency corrected, and
% eddy current corrected, but not yet corrected for linewidth.
% It generates the type of .mat file found in the "matlab_raw_timecourse"
% file. If an error occurs due to svdfid, check comment l.74 or l.100.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% FUNCTIONAL DATA %%
currentdir=pwd;
group='functional';

% SERIE = acquisition file numbers
ACQ = [8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;10 12 14 16 18;8 13 13 15 17;7 11 13 15 17;8 11 13 15 17];
% mouse 8: scan 11 is an outlier and is replaced by repeating scan 13(see in matlab_raw_timecourse, load serie_moving_total_timecourse_old and plot plot(mean(abs((Serie_moving_av)),2)) for visualisation)
Number_of_Scans=size(ACQ);
Nb_of_Repet=544;
Number_of_Points=2048;
Resolution=4; 
Block=8;

for mouse_number=1:10

%load water for eddy current correction
directory_water=strcat(currentdir,filesep,group,filesep,'water_ref',filesep,'mouse',num2str(mouse_number),filesep);
load(strcat(directory_water,'water_eddycurrent.mat'),'Water_tot_Rest')

directory_fid=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'mouse',num2str(mouse_number),filesep);
addpath(strcat(currentdir,filesep,'support_functions'))

%parameters for spectra


fichier = '';
fichier_serie='';
TOTAL_FID=[];
numero='';

%preparing concatenated raw FID

    for k=1:1:Number_of_Scans(2)
    numero = num2str(ACQ(mouse_number,k));

    fichier=strcat(directory_fid,filesep,numero,filesep,'fid');    

    if exist(eval('fichier')) == 2
    fichier=fichier;
    else
    fichier=strcat(directory_fid,filesep,numero,filesep,'ser');
    end

% Loading FID 
[FIDk]=load_array_FID2(fichier,Nb_of_Repet);

TOTAL_FID=[TOTAL_FID FIDk]; 

    end

%building the raw average time course (over 5 blocks)

    for  i=1:1:(Nb_of_Repet-Block)/Resolution +1 
     sliding_window=[];
     vector=[];
        for j=1:1:Number_of_Scans(2)
       sliding_window=(j-1)*Nb_of_Repet+[1:Block]+(i-1)*Resolution;
       vector = [vector sliding_window];       
        end  
      
    rebuilt_FID=TOTAL_FID(:,vector);
    [spec_water fid_phased fid_raw STD_PHI]=phaser(rebuilt_FID);   
    fid_cor=deconv_ECC(fid_phased',Water_tot_Rest');
    
    % when the water residual is too low, svd_fid brings an error: comment
    % the 2 following lines, uncomment "soustraction=fid_cor';" and proceed
    [td_synth, td_diff, params]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    soustraction=td_diff';
%     soustraction=fid_cor'; %when the water residual is too low
    Serie_moving_av(i,:)=soustraction;
    numero_moving_av=num2str(i);

    end 

directory_serie=strcat(directory_fid,filesep,'matlab_raw_timecourse2',filesep);
mkdir(directory_serie)
name_serie=strcat(directory_serie,'serie_moving_average_resolved','.mat');
save(name_serie,'Serie_moving_av','-mat')
clear rebuilt_FID
clear Serie_moving_av


%building the raw total time course (concatenated 5 blocks)
    for  i=1:1:(Nb_of_Repet*Number_of_Scans(2)-Block)/Resolution +1 
    sliding_window=[];
    sliding_window=[1:Block]+(i-1)*Resolution;    
      
    rebuilt_FID=TOTAL_FID(:,sliding_window);
    [spec_water fid_phased fid_raw STD_PHI]=phaser(rebuilt_FID);   
    fid_cor=deconv_ECC(fid_phased',Water_tot_Rest');
    % when the water residual is too low, svd_fid brings an error: comment
    % the 2 following lines, uncomment "soustraction=fid_cor';" and proceed
    [td_synth, td_diff, params]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    soustraction=td_diff';
%     soustraction=fid_cor'; %when the water residual is too low
    Serie_moving_av(i,:)=soustraction;
    numero_moving_av=num2str(i);

    end 

name_serie=strcat(directory_serie,'serie_moving_total_timecourse','.mat');
save(name_serie,'Serie_moving_av','-mat')
clearvars -except currentdir group mouse_number ACQ Number_of_Scans Nb_of_Repet Number_of_Points Resolution Block


end


%% CONTROL DATA %%
currentdir=pwd;
group='control';

% SERIE = acquisition file numbers
ACQ = [17 19 21 23 26;14 16 18 19 21;7 11 13 15 17;9 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 11 13 15 17;8 12 13 15 17;8 11 13 15 17];
Number_of_Scans=size(ACQ);
Nb_of_Repet=544;
Number_of_Points=2048;
Resolution=4; 
Block=8;

for mouse_number=1:10

%load water for eddy current correction
directory_water=strcat(currentdir,filesep,group,filesep,'water_ref',filesep,'mouse',num2str(mouse_number),filesep);
load(strcat(directory_water,'water_eddycurrent.mat'),'Water_tot_Rest')

directory_fid=strcat(currentdir,filesep,group,filesep,'metabolites',filesep,'mouse',num2str(mouse_number),filesep);
addpath(strcat(currentdir,filesep,'support_functions'))

%parameters for spectra


fichier = '';
fichier_serie='';
TOTAL_FID=[];
numero='';

%preparing concatenated raw FID

    for k=1:1:Number_of_Scans(2)
    numero = num2str(ACQ(mouse_number,k));

    fichier=strcat(directory_fid,filesep,numero,filesep,'fid');    

    if exist(eval('fichier')) == 2
    fichier=fichier;
    else
    fichier=strcat(directory_fid,filesep,numero,filesep,'ser');
    end

% Loading FID 
[FIDk]=load_array_FID2(fichier,Nb_of_Repet);

TOTAL_FID=[TOTAL_FID FIDk]; 

    end

%building the raw average time course (over 5 blocks)

    for  i=1:1:(Nb_of_Repet-Block)/Resolution +1 
     sliding_window=[];
     vector=[];
        for j=1:1:Number_of_Scans(2)
       sliding_window=(j-1)*Nb_of_Repet+[1:Block]+(i-1)*Resolution;
       vector = [vector sliding_window];       
        end  
      
    rebuilt_FID=TOTAL_FID(:,vector);
    [spec_water fid_phased fid_raw STD_PHI]=phaser(rebuilt_FID);   
    fid_cor=deconv_ECC(fid_phased',Water_tot_Rest');
    
    % when the water residual is too low, svd_fid brings an error: comment
    % the 2 following lines, uncomment "soustraction=fid_cor';" and proceed
    [td_synth, td_diff, params]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    soustraction=td_diff';
%     soustraction=fid_cor'; %when the water residual is too low
    Serie_moving_av(i,:)=soustraction;
    numero_moving_av=num2str(i);

    end 

directory_serie=strcat(directory_fid,filesep,'matlab_raw_timecourse2',filesep);
mkdir(directory_serie)
name_serie=strcat(directory_serie,'serie_moving_average_resolved','.mat');
save(name_serie,'Serie_moving_av','-mat')
clear rebuilt_FID
clear Serie_moving_av

%building the raw total time course (concatenated 5 blocks)
    for  i=1:1:(Nb_of_Repet*Number_of_Scans(2)-Block)/Resolution +1 
    sliding_window=[];
    sliding_window=[1:Block]+(i-1)*Resolution;    
      
    rebuilt_FID=TOTAL_FID(:,sliding_window);
    [spec_water fid_phased fid_raw STD_PHI]=phaser(rebuilt_FID);   
    fid_cor=deconv_ECC(fid_phased',Water_tot_Rest');
    % when the water residual is too low, svd_fid brings an error: comment
    % the 2 following lines, uncomment "soustraction=fid_cor';" and proceed
    [td_synth, td_diff, params]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    soustraction=td_diff';
%     soustraction=fid_cor'; %when the water residual is too low
    Serie_moving_av(i,:)=soustraction;
    numero_moving_av=num2str(i);

    end 

name_serie=strcat(directory_serie,'serie_moving_total_timecourse','.mat');
save(name_serie,'Serie_moving_av','-mat')
clearvars -except currentdir group mouse_number ACQ Number_of_Scans Nb_of_Repet Number_of_Points Resolution Block


end