%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads the absolute concentrations given in the .COORD files 
% of the LCModel analysis
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentdir=pwd;
addpath(strcat(currentdir,filesep,'support_functions'))

mouse_number_functional=[1 3 5 6 8 9 10];
mouse_number_control=[1 2 3 4 6 7 9];
size_block=135;

for nb_s=mouse_number_functional
group='functional';
resultat=[];
crlb=[];

total_directory = strcat(currentdir,filesep,group,filesep,'mouse_',num2str(nb_s))

for k=1:size_block

file_num=strcat(total_directory,filesep,num2str(k));

if exist(file_num, 'dir') == 7
    
file=strcat(total_directory,filesep,num2str(k),'/fid_asc.COORD');
[error_flag lcmodelresults]=readcoord(file);

for j=1:1:27
if size(lcmodelresults.metabconc(j).name) == size('NAA+NAAG')
    if lcmodelresults.metabconc(j).name == 'NAA+NAAG'
       resultat(k,1) = lcmodelresults.metabconc(j).absconc;
       crlb(k,1) = lcmodelresults.metabconc(j).SD;
       metab(1).name='NAA+NAAG';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Cr+PCr')
    if lcmodelresults.metabconc(j).name == 'Cr+PCr'
       resultat(k,2) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,2) = lcmodelresults.metabconc(j).SD;
       metab(2).name='Cr+PCr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GPC+Cho+PCho')
    if lcmodelresults.metabconc(j).name == 'GPC+PCho+Cho'
       resultat(k,3) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,3) = lcmodelresults.metabconc(j).SD;
       metab(3).name='GPC+PCho+Cho';
    end
end


if size(lcmodelresults.metabconc(j).name) == size('Glu')
    if lcmodelresults.metabconc(j).name == 'Glu'
       resultat(k,4) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,4) = lcmodelresults.metabconc(j).SD;
       metab(4).name='Glu';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Ins')
    if lcmodelresults.metabconc(j).name == 'Ins'
       resultat(k,5) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,5) = lcmodelresults.metabconc(j).SD;
       metab(5).name='Ins';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Tau')
    if lcmodelresults.metabconc(j).name == 'Tau'
       resultat(k,6) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,6) = lcmodelresults.metabconc(j).SD;
       metab(6).name='Tau';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GABA')
    if lcmodelresults.metabconc(j).name == 'GABA'
       resultat(k,7) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,7) = lcmodelresults.metabconc(j).SD;
       metab(7).name='GABA';
    end
end
    
if size(lcmodelresults.metabconc(j).name) == size('MM_fMR')
    if lcmodelresults.metabconc(j).name == 'MM_fMR'
       resultat(k,8) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,8) = lcmodelresults.metabconc(j).SD;
       metab(8).name='MM_fMR';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Gln')
    if lcmodelresults.metabconc(j).name == 'Gln'
       resultat(k,9) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,9) = lcmodelresults.metabconc(j).SD;
       metab(9).name='Gln';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Lac')
    if lcmodelresults.metabconc(j).name == 'Lac'
       resultat(k,10) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,10) = lcmodelresults.metabconc(j).SD;
       metab(10).name='Lac';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Ace')
    if lcmodelresults.metabconc(j).name == 'Ace'
       resultat(k,11) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,11) = lcmodelresults.metabconc(j).SD;
       metab(11).name='Ace';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('PCr')
    if lcmodelresults.metabconc(j).name == 'PCr'
       resultat(k,12) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,12) = lcmodelresults.metabconc(j).SD;
       metab(12).name='PCr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Cr')
    if lcmodelresults.metabconc(j).name == 'Cr'
       resultat(k,13) = lcmodelresults.metabconc(j).absconc;
       crlb(k,13) = lcmodelresults.metabconc(j).SD;
       metab(13).name='Cr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('NAA')
    if lcmodelresults.metabconc(j).name == 'NAA'
       resultat(k,14) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,14) = lcmodelresults.metabconc(j).SD;
       metab(14).name='NAA';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('NAAG')
    if lcmodelresults.metabconc(j).name == 'NAAG'
       resultat(k,15) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,15) = lcmodelresults.metabconc(j).SD;
       metab(15).name='NAAG';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GPC')
    if lcmodelresults.metabconc(j).name == 'GPC'
       resultat(k,16) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,16) = lcmodelresults.metabconc(j).SD;
       metab(16).name='GPC';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GSH')
    if lcmodelresults.metabconc(j).name == 'GSH'
       resultat(k,17) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,17) = lcmodelresults.metabconc(j).SD;
       metab(17).name='GSH';
    end
end


if size(lcmodelresults.metabconc(j).name) == size('Gly')
   if lcmodelresults.metabconc(j).name == 'Gly'
      resultat(k,18) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,18) = lcmodelresults.metabconc(j).SD;
      metab(18).name='Gly';
   end
end 

if size(lcmodelresults.metabconc(j).name) == size('Asp')
   if lcmodelresults.metabconc(j).name == 'Asp'
      resultat(k,19) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,19) = lcmodelresults.metabconc(j).SD;
      metab(19).name='Asp';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('Glc')
   if lcmodelresults.metabconc(j).name == 'Glc'
      resultat(k,20) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,20) = lcmodelresults.metabconc(j).SD;
      metab(20).name='Glc';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('Glu+Gln')
   if lcmodelresults.metabconc(j).name == 'Glu+Gln'
      resultat(k,21) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,21) = lcmodelresults.metabconc(j).SD;
      metab(21).name='Glu+Gln';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('PCho')
   if lcmodelresults.metabconc(j).name == 'PCho'
      resultat(k,22) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,22) = lcmodelresults.metabconc(j).SD;
      metab(22).name='PCho';
   end
end


end
end

end
functional(nb_s).table=resultat;
functional(nb_s).crlb=crlb;
end



for nb_s=mouse_number_control
group='control';
resultat=[];
crlb=[];

total_directory=strcat(currentdir,filesep,group,filesep,'mouse_',num2str(nb_s))

for k=1:size_block

file_num=strcat(total_directory,filesep,num2str(k));

if exist(file_num, 'dir') == 7
    
file=strcat(total_directory,filesep,num2str(k),filesep,'fid_asc.COORD');
[error_flag lcmodelresults]=readcoord(file);

for j=1:1:27
if size(lcmodelresults.metabconc(j).name) == size('NAA+NAAG')
    if lcmodelresults.metabconc(j).name == 'NAA+NAAG'
       resultat(k,1) = lcmodelresults.metabconc(j).absconc;
       crlb(k,1) = lcmodelresults.metabconc(j).SD;
       metab(1).name='NAA+NAAG';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Cr+PCr')
    if lcmodelresults.metabconc(j).name == 'Cr+PCr'
       resultat(k,2) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,2) = lcmodelresults.metabconc(j).SD;
       metab(2).name='Cr+PCr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GPC+Cho+PCho')
    if lcmodelresults.metabconc(j).name == 'GPC+PCho+Cho'
       resultat(k,3) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,3) = lcmodelresults.metabconc(j).SD;
       metab(3).name='GPC+PCho+Cho';
    end
end


if size(lcmodelresults.metabconc(j).name) == size('Glu')
    if lcmodelresults.metabconc(j).name == 'Glu'
       resultat(k,4) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,4) = lcmodelresults.metabconc(j).SD;
       metab(4).name='Glu';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Ins')
    if lcmodelresults.metabconc(j).name == 'Ins'
       resultat(k,5) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,5) = lcmodelresults.metabconc(j).SD;
       metab(5).name='Ins';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Tau')
    if lcmodelresults.metabconc(j).name == 'Tau'
       resultat(k,6) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,6) = lcmodelresults.metabconc(j).SD;
       metab(6).name='Tau';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GABA')
    if lcmodelresults.metabconc(j).name == 'GABA'
       resultat(k,7) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,7) = lcmodelresults.metabconc(j).SD;
       metab(7).name='GABA';
    end
end
    
if size(lcmodelresults.metabconc(j).name) == size('MM_fMR')
    if lcmodelresults.metabconc(j).name == 'MM_fMR'
       resultat(k,8) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,8) = lcmodelresults.metabconc(j).SD;
       metab(8).name='MM_fMR';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Gln')
    if lcmodelresults.metabconc(j).name == 'Gln'
       resultat(k,9) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,9) = lcmodelresults.metabconc(j).SD;
       metab(9).name='Gln';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Lac')
    if lcmodelresults.metabconc(j).name == 'Lac'
       resultat(k,10) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,10) = lcmodelresults.metabconc(j).SD;
       metab(10).name='Lac';
    end
end 

if size(lcmodelresults.metabconc(j).name) == size('Ace')
    if lcmodelresults.metabconc(j).name == 'Ace'
       resultat(k,11) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,11) = lcmodelresults.metabconc(j).SD;
       metab(11).name='Ace';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('PCr')
    if lcmodelresults.metabconc(j).name == 'PCr'
       resultat(k,12) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,12) = lcmodelresults.metabconc(j).SD;
       metab(12).name='PCr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('Cr')
    if lcmodelresults.metabconc(j).name == 'Cr'
       resultat(k,13) = lcmodelresults.metabconc(j).absconc;
       crlb(k,13) = lcmodelresults.metabconc(j).SD;
       metab(13).name='Cr';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('NAA')
    if lcmodelresults.metabconc(j).name == 'NAA'
       resultat(k,14) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,14) = lcmodelresults.metabconc(j).SD;
       metab(14).name='NAA';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('NAAG')
    if lcmodelresults.metabconc(j).name == 'NAAG'
       resultat(k,15) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,15) = lcmodelresults.metabconc(j).SD;
       metab(15).name='NAAG';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GPC')
    if lcmodelresults.metabconc(j).name == 'GPC'
       resultat(k,16) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,16) = lcmodelresults.metabconc(j).SD;
       metab(16).name='GPC';
    end
end

if size(lcmodelresults.metabconc(j).name) == size('GSH')
    if lcmodelresults.metabconc(j).name == 'GSH'
       resultat(k,17) = lcmodelresults.metabconc(j).absconc; 
       crlb(k,17) = lcmodelresults.metabconc(j).SD;
       metab(17).name='GSH';
    end
end


if size(lcmodelresults.metabconc(j).name) == size('Gly')
   if lcmodelresults.metabconc(j).name == 'Gly'
      resultat(k,18) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,18) = lcmodelresults.metabconc(j).SD;
      metab(18).name='Gly';
   end
end 

if size(lcmodelresults.metabconc(j).name) == size('Asp')
   if lcmodelresults.metabconc(j).name == 'Asp'
      resultat(k,19) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,19) = lcmodelresults.metabconc(j).SD;
      metab(19).name='Asp';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('Glc')
   if lcmodelresults.metabconc(j).name == 'Glc'
      resultat(k,20) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,20) = lcmodelresults.metabconc(j).SD;
      metab(20).name='Glc';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('Glu+Gln')
   if lcmodelresults.metabconc(j).name == 'Glu+Gln'
      resultat(k,21) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,21) = lcmodelresults.metabconc(j).SD;
      metab(21).name='Glu+Gln';
   end
end

if size(lcmodelresults.metabconc(j).name) == size('PCho')
   if lcmodelresults.metabconc(j).name == 'PCho'
      resultat(k,22) = lcmodelresults.metabconc(j).absconc; 
      crlb(k,22) = lcmodelresults.metabconc(j).SD;
      metab(22).name='PCho';
   end
end


end
end

end
control(nb_s).table=resultat;
control(nb_s).crlb=crlb;
end

save(strcat(currentdir,filesep,'functional',filesep,'C1_Concentrations.mat'),'functional')
save(strcat(currentdir,filesep,'control',filesep,'C1_Concentrations.mat'),'control')
save(strcat(currentdir,filesep,'C1_Metabolites_List.mat'),'metab')