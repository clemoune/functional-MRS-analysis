
function [Matrix_Vox_Spectro_EPI] = Spectro_Vox_EPI(directory_EPI,directory_spectro)

% This function creates a 256*256*10 matrix (T2_TurboRARE size)
% and fills with 1 voxels belonging to the spectroscopic voxel. The
% center of the matrix correponds to the center of the anatomical image. 

method_T2=strcat(directory_EPI,filesep,'method');
method_spectro=strcat(directory_spectro,filesep,'method');

T2_center=zeros; %coordinates in mm extracted from method T2 file
T2_center_1=''; 
T2_center_2='';
T2_center_3='';

spectro_vox_center=zeros; %coordinates in mm extracted from method spectro file
spectro_vox_center_1='';
spectro_vox_center_2='';
spectro_vox_center_3='';

spectro_vox_dimension=zeros; %coordinates in mm extracted from method spectro file
spectro_vox_dimension_1='';
spectro_vox_dimension_2='';
spectro_vox_dimension_3='';

% matrix=zeros(256,256,10);
% matrix_center=[128 5 128];
% matrix=zeros(200,160,20);
% matrix_center=[100 10 80];
matrix=zeros(110,96,8);
matrix_center=[55.5 48.5 4.5];
spectro_center_matrix=zeros; %center of spectroscopic voxel in matrix
spectro_corner_matrix=zeros; %initial corner of spectroscopic voxel in matrix

% Anatomical image center
fileid=fopen(method_T2,'r');

p=[];
p=fscanf(fileid,'%s',1);
while strcmp(p,'##$PVM_SPackArrReadOffset=(')==0
            p=fscanf(fileid,'%s',1);
end
   
for j=1:2
  p=fscanf(fileid,'%s',1);
end
T2_center_1=fscanf(fileid,'%s',1);

while strcmp(p,'##$PVM_SPackArrPhase1Offset=(')==0
            p=fscanf(fileid,'%s',1);
end
   
for j=1:2
  p=fscanf(fileid,'%s',1);
end
T2_center_2=fscanf(fileid,'%s',1);

while strcmp(p,'##$PVM_SPackArrSliceOffset=(')==0
            p=fscanf(fileid,'%s',1);
end
   
for j=1:2
  p=fscanf(fileid,'%s',1);
end
T2_center_3=fscanf(fileid,'%s',1);


fclose(fileid);

% Spectroscopic voxel size 
fileid2=fopen(method_spectro,'r');

p=[];
p=fscanf(fileid2,'%s',1);
while strcmp(p,'##$PVM_VoxArrSize=(')==0
            p=fscanf(fileid2,'%s',1);
end
   
for j=1:3
  p=fscanf(fileid2,'%s',1);
end

spectro_vox_dimension_1=fscanf(fileid2,'%s',1);
spectro_vox_dimension_2=fscanf(fileid2,'%s',1);
spectro_vox_dimension_3=fscanf(fileid2,'%s',1);

fclose(fileid2);


% Spectroscopic voxel center
fileid3=fopen(method_spectro,'r');

p=[];
p=fscanf(fileid3,'%s',1);
while strcmp(p,'##$PVM_VoxArrPosition=(')==0
            p=fscanf(fileid3,'%s',1);
end
for j=1:3
  p=fscanf(fileid3,'%s',1);
end

spectro_vox_center_1=fscanf(fileid3,'%s',1);
spectro_vox_center_2=fscanf(fileid3,'%s',1);
spectro_vox_center_3=fscanf(fileid3,'%s',1);

fclose(fileid3);

% Matrix creation

T2_center(1)=str2num(T2_center_1)/0.145;
T2_center(2)=str2num(T2_center_2)/0.145;
T2_center(3)=str2num(T2_center_3)/0.5;

spectro_vox_center(1)=str2num(spectro_vox_center_1)/0.145;
spectro_vox_center(2)=str2num(spectro_vox_center_2)/0.145;
spectro_vox_center(3)=str2num(spectro_vox_center_3)/0.5;

spectro_vox_dimension(1)=str2num(spectro_vox_dimension_1)/0.145;
spectro_vox_dimension(2)=str2num(spectro_vox_dimension_2)/0.145;
spectro_vox_dimension(3)=str2num(spectro_vox_dimension_3)/0.5;

spectro_center_matrix(1)=matrix_center(1)+(spectro_vox_center(1)-T2_center(1)); 
spectro_center_matrix(2)=matrix_center(2)-(spectro_vox_center(2)-T2_center(2)); 
spectro_center_matrix(3)=matrix_center(3)-(spectro_vox_center(3)+T2_center(3)); 

spectro_corner_matrix(1)=spectro_center_matrix(1)-spectro_vox_dimension(1)/2;
spectro_corner_matrix(2)=spectro_center_matrix(2)-spectro_vox_dimension(2)/2;
spectro_corner_matrix(3)=spectro_center_matrix(3)-spectro_vox_dimension(3)/2;

    %integer
spectro_vox_dimension=round(spectro_vox_dimension);
spectro_corner_matrix=round(spectro_corner_matrix);
spectro_center_matrix=round(spectro_center_matrix);

matrix(spectro_corner_matrix(1)+1:spectro_corner_matrix(1)+spectro_vox_dimension(1),spectro_corner_matrix(2)+1:spectro_corner_matrix(2)+spectro_vox_dimension(2),spectro_corner_matrix(3)+1:spectro_corner_matrix(3)+spectro_vox_dimension(3))=ones(spectro_vox_dimension(1),spectro_vox_dimension(2),spectro_vox_dimension(3));
Matrix_Vox_Spectro_EPI=matrix;


[image_EPI im] = read_2dseqbis(strcat(directory_EPI,'pdata\1\'));
image_EPI_db=double(image_EPI(:,:,1,1,spectro_center_matrix(3)+1));
imagesc(rot90((image_EPI_db + image_EPI_db.*Matrix_Vox_Spectro_EPI(:,:,spectro_center_matrix(3)+1)),3))
end
% 
