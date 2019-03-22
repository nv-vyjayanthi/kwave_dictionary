% filter and reconstruct from the PA response 

% run this to get data
kwave_BuildFromImpulse_PAresp

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

% case1: using filteres&integrated impulse response
NFFT = 2*numel(t_array) + 1;
psf = impulse_resp_int.';
psf = psf./max(psf);
psf_dct = dct(psf);
psf_ft = abs(ifftshift(fft(psf,NFFT)));
psf_hat = psf_ft;
psf_hat = psf_hat./max(psf_hat);
psf_hat = psf_hat(round(numel(psf_ft)/2)+1:end);

% % case2: using only integrated impulse response
% NFFT = 2*numel(t_array) + 1;
% psf = abs(cumtrapz(t_array,impulse_resp)).';
% % psf(psf<0) = 0;
% psf = psf./max(psf);
% psf_dct = dct(psf);
% temp_psf = fftshift(psf);
% psf_ft = abs(ifftshift(fft(psf,NFFT)));
% psf_hat = psf_ft;
% psf_hat = psf_hat./max(psf_hat);
% psf_hat = psf_hat(round(numel(psf_ft)/2)+1:end);

freq = linspace(-1,1,numel(t_array)).*1/(2*kgrid.dt);
figure;subplot(1,2,1); plot(psf); title('PSF')
subplot(1,2,2); plot(freq/1e6,psf_hat); title('PSF FT'); xlabel('MHz')

%% reconstruction algorithm for uniform object

% figure;
% subplot(1,2,1);plot(freq(max_ind(2):max_ind(2)+kgrid.Nt-1),psf);title('PSF'); xlabel('Freq (MHz)')
% subplot(1,2,2);plot(freq(max_ind(2):max_ind(2)+kgrid.Nt-1),psf_hat);title('PSF ft'); xlabel('Freq (MHz)')

y0 = obj_resp_int.';
Niter=500; L=sqrt(2); la = 0.3; 
p0rec = deconviter(y0, zeros(size(y0)), psf_hat, L, Niter, la, 'fastnet' );

% plot i/p and o/p of reconstruction algo
figure;
subplot(1,2,1); plot(t_array/1e-6,y0); xlabel('t \mus'); title('Integrated PA resp')
subplot(1,2,2); plot(t_array/1e-6,p0rec); xlabel('t \mus'); title('Reconstruction')
suptitle('Uniform Illum')

%% reconstruction algorithm for speckle object

% setup forward kernel
A = psf_hat*ones(1,N_speckle); % forward kernel in k-space
Y = speckle_resp_int;

% setup parameters for iterative reconstruction algorithm 
Niter=900;
L = sqrt(2); % bound Lipschitzconstant 
la1 = 0.10;  
la2 = 0.0001;
% alg='fista';  % elastic-net ell^l-ell^2 term
alg = 'bfista'; % elastic-net ell^l-ell^2 term
% ... 'fista'   R = la1 * ||PREC||_1 + la2/2 * ||PREC||_2^2 ... sparse
% ... 'b-fista' R = la1 * ||PREC||_21 + la2/2 * ||PREC||_2^2 ... blocksparse 

Prec = reconiter1Dv2(Y,zeros(size(Y)), A, L, Niter, la1, la2, alg);
piter = sum (Prec,2); 

% plot i/p and o/p of reconstruction of algorithm
figure
subplot(1,2,1); plot(t_array/1e-6,mean(Y,2)); xlabel('t \mus'); title('Integrated PA resp')
subplot(1,2,2); plot(t_array/1e-6,piter/max(piter)); xlabel('t \mus'); title('Reconstruction')
suptitle('Speckle Illum')



