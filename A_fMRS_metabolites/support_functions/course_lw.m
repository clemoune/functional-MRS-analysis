function [linewidth_c,height_c] = course_lw(time_spec,nb_time_points)
%UNTITLED3 Summary of this function goes here
% 

linewidth=zeros(1,nb_time_points);
height=zeros(1,nb_time_points);
integration=zeros(1,nb_time_points);

   for j=1:nb_time_points
    spec_abs_c=abs(fftshift(fft([time_spec(:,j)' zeros(1,4096)]))); 
    spec_abs_naa_c=spec_abs_c(1,4600:4800);   %NAA
    spec_abs_tcr_c=spec_abs_c(1,4020:4170);   %tCr
    spec_abs_tcho_c=spec_abs_c(1,3850:4040);   %tCho
    [M1 I1]=max(spec_abs_naa_c);
    [M2 I2]=max(spec_abs_tcr_c);
    [M3 I3]=max(spec_abs_tcho_c);
    
    half_width1=M1/2;
    k1=I1;
    half_width2=M2/2;
    k2=I2;
    half_width3=M3/2;
    k3=I3;

    
    while spec_abs_naa_c(1,k1)>half_width1
        k1=k1+1;
    end     
    k1_2=I1;
    while spec_abs_naa_c(1,k1_2)>half_width1
        k1_2=k1_2-1;
    end

    
    while spec_abs_tcr_c(1,k2)>half_width2
        k2=k2+1;
    end     
    k2_2=I2;
    while spec_abs_tcr_c(1,k2_2)>half_width2
        k2_2=k2_2-1;
    end

    
    while spec_abs_tcho_c(1,k3)>half_width3
        k3=k3+1;
    end    
    k3_2=I3;
    while spec_abs_tcho_c(1,k3_2)>half_width3
        k3_2=k3_2-1;
    end
    
    linewidth_c(1,j)=(k1+k2+k3-(k1_2+k2_2+k3_2))/3;
    height_c(1,j)=(M1+M2+M3)/3;
    integration_c(1,j)=(sum(spec_abs_naa_c)+sum(spec_abs_tcr_c)+sum(spec_abs_tcho_c))/3;
    end

end

