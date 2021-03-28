%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasage des spectres individuels
% Julien Valette september 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec_water fid_phased fid_raw STD_PHI]=phaser(fid_matrix)

global spec_ref spec_p_lb fid_p_lb time factor_for_phi NAA_inf NAA_sup

nb_pts_cplx=size(fid_matrix,1);


dw=1/4000;

lb=3; % lb for phase correction only, not on the final fid


% with 2048 points, BW=4000Hz 
NAA_inf=550;
NAA_sup=800;

NS=size(fid_matrix,2);
%NS=45;
%%% Generating the raw sum FID %%%%%%%%%%%%

fid_raw=sum(fid_matrix,2);

%%% Generating the reference spectrum for correction %%%%%%%%%%%%

time=((0:nb_pts_cplx-1)*dw)';

fid_ref=fid_raw.*exp(-lb*time)/NS;
% fid_ref=fid_matrix(:,2).*exp(-lb*time);
spec_ref=fftshift(fft(fid_ref));
spec_ref=spec_ref(NAA_inf:NAA_sup);
spec_water=zeros(NS,length(spec_ref));
%NS_real=0;
%fidtotale=zeros(nb_pts_cplx,1);
%fidpost=zeros(nb_pts_cplx,1);

fid_phased=zeros(nb_pts_cplx,1);
list_phase=zeros(NS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual spectrum correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options=optimset;
optjv=optimset(options,'MaxFunEvals',1e6,'display', 'off');

for p=1:NS
        
        fid_p=squeeze(fid_matrix(:,p));
        
    	fid_p_lb=fid_p.*exp(-lb*time);
	
    	%%% Frequency correction %%% 
        warning off all    
        y=lsqnonlin('cout_nu',[0 1],[-80 0],[80 3],optjv);
        warning on all   
        delta_nu=y(1);
        factor_for_phi=y(2);
	
        correction=exp(2*pi*1i*delta_nu*time);
        
        fid_p=fid_p.*correction;

    	%%% Phase correction %%% 
        
        fid_p_lb=fid_p.*exp(-lb*time);	
        spec_p_lb=fftshift(fft(fid_p_lb));
        spec_p_lb=spec_p_lb(NAA_inf:NAA_sup);
        
        warning off all 
        phi=lsqnonlin('cout_phase_spec',0,-pi,pi,optjv);
    	warning on all 
        
        list_phase(p)=phi*180/pi;
	
        fid_p=exp(1i*phi)*fid_p;
        
        %%% Frequency correction %%% 
        warning off all     
        y=lsqnonlin('cout_nu',[0 1],[-80 0],[80 3],optjv);
        warning on all 
        delta_nu=y(1);
        factor_for_phi=y(2);
	
        correction=exp(2*pi*1i*delta_nu*time);
        
        fid_p=fid_p.*correction;
        
%         figure
%         plot(real(fftshift(fft(fid_p))))
%         hold on
        
        fid_phased=fid_phased + fid_p;

        spec_water(p,:)=spec_ref;

end

% g=figure
% figure(g)
% plot(list_phase);

MEAN_PHI=mean(list_phase);
STD_PHI=std(list_phase);

%  f=figure;
%  figure(f);
%  hold on;
%  plot(real(fftshift(fft(fid_raw'))));
%  plot(real(fftshift(fft(fid_phased'))),'r');
%  hold off

fid_raw=fid_raw';
fid_phased=fid_phased';
