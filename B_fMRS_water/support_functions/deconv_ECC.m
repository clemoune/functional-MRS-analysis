% Julien Valette, Eddy Current correction
% reference Klose, 1990, MRM

% Inputs
% fid: the FID of the spectrum to correct
% fid_ref: the FID of water acquired in the exact same conditions (except WS off)

% Outputs
% fid_cor: eddy current corrected FID


function fid_cor=deconv_ECC(fid,fid_ref)

cor_ECC=-angle(fid_ref);

fid_cor=fid.*exp(1i*cor_ECC);

% f=figure;
% figure(f);
% hold on;
% plot(real(fftshift(fft(fid))));
% plot(real(fftshift(fft(fid_cor))),'r');