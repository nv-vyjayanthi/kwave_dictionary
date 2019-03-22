% filter and reconstruct from the PA response 

% run this to get data
kwave_BuildFromImpulse_PAresp_v1

%% reconstruction cont.

% NFFT = round((2+1/2)*kgrid.Nt); % number of pts to do fft with
% freq = linspace(-1,1,NFFT)*fs/2;
% psf = impulse_resp_int.';
% impulse_resp_ft = ifftshift(fft(fftshift(impulse_resp_int)));
% psf_fft_abs = abs(impulse_resp_ft);
% % psf to feed to reconstruction algo fcn (along time)
% % psf_hat = psf_fft_abs(round(NFFT/2)+1:round(NFFT/2)+numel(t_array)).';
% psf = psf./max(psf);

% figure; plot(psf_fft_abs); title('Entire FT')

% % case1: using filteres&integrated impulse response
% NFFT = 2*numel(t_array) + 1;
% psf = impulse_resp_int.';
% psf = psf./max(psf);
% psf_dct = dct(psf);
% % psf_ft = abs(ifftshift(fft(psf,NFFT)));
% % psf_ft = real(ifftshift(fft(fftshift(psf))));
% psf_ft = ifftshift(fft(psf,NFFT));
% psf_ft = abs(psf_ft)/2;
% psf_hat = psf_ft/max(psf_ft);
% psf_hat = psf_hat./max(psf_hat);
% psf_hat = psf_hat(round(numel(psf_ft)/2)+1:end);

% case2: using only integrated impulse response
NFFT = 2*numel(t_array) + 1;
% psf = abs(cumtrapz(t_array,impulse_resp)).';
psf = impulse_resp_int.';
% psf(psf<0) = 0;
psf = psf./max(psf);
% psf_dct = dct(psf);
psf_ft = abs(ifftshift(fft(fftshift([zeros(round(NFFT/4)+1,1); psf; zeros(round(NFFT/4),1)]))));
psf_hat = psf_ft;
psf_hat = psf_hat./max(psf_hat);
psf_hat = psf_hat(round(numel(psf_ft)/2)+1:end);

freq = linspace(-1,1,NFFT).*fs/2; % fs in MHz
figure;subplot(1,2,1); plot(t_array/1e-9,psf); title('PSF'); xlabel('t(ns)')
subplot(1,2,2); plot(psf_hat); title('PSF FT')

%% check if psf is actually kernel

% obj = squeeze(p00_new(det_pos(1),det_pos(2),:));
% obj = [zeros(1,round(numel(t_array)/2 - numel(obj)/2)).'; obj; zeros(1,round(numel(t_array)/2 - numel(obj)/2)).'];
% img = idct(dct(obj).*psf_hat);
% figure;plot(img);
% hold on

obj_pos_t = round((kgrid.dz*obj_pos_z/medium.sound_speed)/kgrid.dt) + zero_pad - round(det_pos(3)*kgrid.dz./medium.sound_speed/kgrid.dt);
obj = zeros([length(t_array),1]);
obj(obj_pos_t) = 1;
img_t_dct = idct(dct(obj).*psf_hat);

obj_t_ft = ifftshift(fft(fftshift(obj)));
psf_ft_conv = ifftshift(fft(fftshift(psf)));
img_t_ft = ifftshift(ifft(fftshift(obj_t_ft.*psf_ft_conv)));

figure; 
hold on
plot(img_t_dct/max(img_t_dct));
plot((obj_resp_int/max(obj_resp_int)))
plot(img_t_ft./max(img_t_ft));
legend('DCT convlution','True','Conv')
title('Comparing results of DCT,FFT and k-wave ')

% figure; plot(t_array/1e-9, obj); title('Obj in time'); xlabel('ns')

%% reconstruction algorithm for uniform object

y0 = obj_resp_int.';
% y0 = y0/max(y0);
Niter=1000; L=sqrt(2); la = 1e-12; 
p0rec = deconviter(y0, zeros(size(y0)), psf_hat, L, Niter, la, 'fista' );

% plot i/p and o/p of reconstruction algo
figure;
subplot(1,2,1); plot(t_array/1e-6,y0); xlabel('t \mus'); title('Integrated PA resp')
subplot(1,2,2); plot(t_array/1e-6,p0rec/max(p0rec)); xlabel('t \mus'); title('Reconstruction')
suptitle('Uniform Illum')

%% reconstruction algorithm for speckle object

% setup forward kernel
A = psf_hat*ones(1,N_speckle); % forward kernel in k-space
Y = speckle_resp_int;
% Y = Y./max(Y);

% setup parameters for iterative reconstruction algorithm 
Niter=1000;
L = sqrt(2); % bound Lipschitzconstant 
la1 = 1e-13;  
la2 = 0.0001;
% alg='fista';  % elastic-net ell^l-ell^2 term
alg = 'Bfista'; % elastic-net ell^l-ell^2 term
% ... 'fista'   R = la1 * ||PREC||_1 + la2/2 * ||PREC||_2^2 ... sparse
% ... 'b-fista' R = la1 * ||PREC||_21 + la2/2 * ||PREC||_2^2 ... blocksparse 

Prec = reconiter1Dv2(Y,zeros(size(Y)), A, L, Niter, la1, la2, alg);
piter = sum (Prec,2); 

% plot i/p and o/p of reconstruction of algorithm
% figure
% subplot(1,2,1); plot(t_array/1e-6,mean(Y,2)); xlabel('t \mus'); title('Integrated PA resp')
% subplot(1,2,2); plot(t_array/1e-6,piter/max(piter)); xlabel('t \mus'); title('Reconstruction')
% suptitle('Speckle Illum')

%

figure
hold on
plot(t_array/1e-9,obj,'DisplayName','Obj(t)','LineWidth',0.75,'Color','black')
plot(t_array/1e-9,y0,'DisplayName','Uniform Resp','LineStyle','--','Color','blue','LineWidth',0.75)
plot(t_array/1e-9,mean(Y,2),'DisplayName','Speckle Resp (mean)','LineStyle','--','Color',[0 0.9 0],'LineWidth',0.75)
plot(t_array/1e-9,p0rec/max(p0rec),'DisplayName','Uniform Resp Recons','LineWidth',1.25,'Color',[0 0.5 1])
plot(t_array/1e-9,piter/max(piter),'DisplayName','Speckle Resp Recons (mean)','LineWidth',1.25,'Color',[0 0.9 0.75])
grid on; legend('show','Location','southeast')
xlabel('t(ns)')
title({['Recinstruction results'],['Speckle width = ' num2str(speckle_FWHM*kgrid.dz/1e-6/1.5) 'ns']})
xlim([0 100])
