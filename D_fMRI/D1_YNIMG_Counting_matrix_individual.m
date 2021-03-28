%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes the voxel activity in the voxel for each mouse. It
% generates the D1 files.
% by C. Ligneul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentdir=pwd;
ACQ_EPI1=[4 4 4 4 4 4 3];
ACQ_EPI2=[18*ones(1,7)];
ACQ_VOX=[12*ones(1,7)];
mouse_number_functional=[1 3 5 6 8 9 10];
addpath(strcat(currentdir,filesep,'support_functions'))

for ki=1:7
directoryEPI=strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_EPI1(ki)), filesep)
directoryspectro=strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_VOX(ki)), filesep)

if ki==7
    [spec_vox_EPI] = Spectro_Vox_EPI_dirinv(directoryEPI,directoryspectro);
else
    [spec_vox_EPI] = Spectro_Vox_EPI(directoryEPI,directoryspectro);
end
[image_EPI im] = read_2dseqbis(strcat(directoryEPI,'pdata\1\'));

spmfile=niftiread(strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_EPI1(ki)), filesep,'Processed\Simple individual analysis\spmT_0001.nii'));
Percentage_Active_Voxel=zeros;
for i = 1:8
    if sum(sum(spec_vox_EPI(:,:,i)))>0
%         Fifi=figure
%         imagesc(rot90(double(image_EPI(:,:,1,1,i)).*spec_vox_EPI(:,:,i)+double(image_EPI(:,:,1,1,i)),3))
%         savefig(Fifi,strcat(directoryspectro,filesep,'Superimposition_Vox_EPI_end',num2str(i),'.fig'))
%         
        dd=spmfile(:,:,i);
        ssv=spec_vox_EPI(:,:,i);
        dd_thresh_indices=find(dd>2.34);


        [row_act,col_act] = ind2sub([110 96],dd_thresh_indices);
 
        matrix_act=zeros(110,96);
        for k=1:size(row_act,1)
            matrix_act(row_act(k),col_act(k))=1;
        end
  Percentage_Active_Voxel(i)=sum(sum(matrix_act.*ssv))/sum(sum(ssv));   
    end
end
Mean_Fraction=mean(nonzeros(Percentage_Active_Voxel));
EPI1_Activity_Voxel(mouse_number_functional(ki)).Active_Fraction=Mean_Fraction;
end
save(strcat(currentdir,filesep,'D1_EPI1_Activity_Vox.mat'),'EPI1_Activity_Voxel')

for ki=1:7
directoryEPI=strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_EPI2(ki)), filesep)
directoryspectro=strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_VOX(ki)), filesep)

if ki==7
    [spec_vox_EPI] = Spectro_Vox_EPI_dirinv(directoryEPI,directoryspectro);
else
    [spec_vox_EPI] = Spectro_Vox_EPI(directoryEPI,directoryspectro);
end
[image_EPI im] = read_2dseqbis(strcat(directoryEPI,'pdata\1\'));

spmfile=niftiread(strcat(currentdir,filesep,'mouse',num2str(mouse_number_functional(ki)),filesep,num2str(ACQ_EPI2(ki)), filesep,'Processed\Simple individual analysis\spmT_0001.nii'));
Percentage_Active_Voxel=zeros;
for i = 1:8
    if sum(sum(spec_vox_EPI(:,:,i)))>0
%         Fifi=figure
%         imagesc(rot90(double(image_EPI(:,:,1,1,i)).*spec_vox_EPI(:,:,i)+double(image_EPI(:,:,1,1,i)),3))
%         savefig(Fifi,strcat(directoryspectro,filesep,'Superimposition_Vox_EPI_end',num2str(i),'.fig'))
        
        dd=spmfile(:,:,i);
        ssv=spec_vox_EPI(:,:,i);
        dd_thresh_indices=find(dd>2.34);


        [row_act,col_act] = ind2sub([110 96],dd_thresh_indices);
 
        matrix_act=zeros(110,96);
        for k=1:size(row_act,1)
            matrix_act(row_act(k),col_act(k))=1;
        end
  Percentage_Active_Voxel(i)=sum(sum(matrix_act.*ssv))/sum(sum(ssv));
    end
end
Mean_Fraction=mean(nonzeros(Percentage_Active_Voxel));
EPI2_Activity_Voxel(mouse_number_functional(ki)).Active_Fraction=Mean_Fraction;
end
save(strcat(currentdir,filesep,'D1_EPI2_Activity_Vox.mat'),'EPI1_Activity_Voxel')