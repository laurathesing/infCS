%% Compressed sensing with a GS matrix in 1D (uses spgl1)

%  The example produces GS reconstruction and truncated Fourier
%  representation of a discontinuous funtion when Fourier samples are taken uniformly;
%  Two sampling masks are used: 1. multilevel random sampling map, 
%  and 2. sampling map which indexes only samples of low frequencies.

clear all; 

%% Input Parameters (change as needed)

a = 6; %number of vanishing moments; it can be in {1,2,3,4,5,6,7}
eps = 1; %require 1/eps to be a natural number; eps=1 gives Nyquist rate difference of impact of eps and S?
S = 1/2; %Fourier samples range: -S*2^R/eps, ... , 2^R*S/eps-1
R = 13; %maximum scale of wavelet coefficients

M = S*2^R/eps;

scalMN = 1/8;
M0=M*scalMN; % sampling bandwidth N
Mmax = 2*M;

j0 =4; % adapt according to vm
sub_samplingrate = 0.5;

%% Precompute the fourier transforms of the scaling function, right boundary scaling functions, and left boundary scaling functions.

disp('Precomputing Fourier transforms of the scaling functions...')
if a==1  % the case of Haar wavelets
    ft_sca= haar_phi_ft(((M:-1:-M+1)./(2^R/eps))');
else     % the case of boundary corrected wavelets
    [ft_sca_L, ft_sca, ft_sca_R] = CDJV_Setup(a, eps, S, R);
end

%% create two masks for the Fourier coefficients

% uncomment for different sampling matrix type
% p1=1;p2=1.4;
% masktemp = GetSampleMask1D(M0,ceil(0.25*M0), p1,p2, 30, 200);
% mask = zeros(Mmax,1);
% mask(Mmax/2-M0/2:Mmax/2+M0/2-1,1) = masktemp;
% % 
% sparsity = sum(mask)/(M0);
% %fprintf('Taking %.2f%% of the samples...\n',sparsity*100)
% figure('Name','Mask 1'); 
% imagesc(-mask'); colormap gray; 

[idx, ~] = sph1_rect2(M0, M0, ceil(sub_samplingrate*M0)/2, j0); % uses the subsampling matrix from cww
idx = idx';
idx_plot = zeros(Mmax,1);
idx_plot(Mmax/2+idx) =1;
idx_plot(1:Mmax/2) = fliplr(idx_plot(Mmax/2+1:Mmax)');
mask = idx_plot;

%figure; imagesc(-mask'); colormap gray;
sum(mask(:))    

samp_sum = floor(sum(mask(:))/2);
mask_centre =zeros(2*M,1);
mask_centre(M-samp_sum:M+samp_sum) = 1;
% figure('Name','Mask 2');
% imagesc(-mask_centre'); colormap gray;

    
%% Define a function and its Fourier samples

disp('Calculating the Fourier samples...')

%f = @(x) sin(pi*x) + (sin(5*pi*x)-1).*(x>0.3);%.*(x<0.58);
%f = @(x) x.^5;
f = @(x) sin(3*pi*x) + sin(5*pi*x).*(x > 0.51);
xylims=[-0.02 1.02 -1.5 2];

%omega_main = Samples(eps, M, 'sin(pi*x)', 0, 1) + Samples(eps, M, 'sin(5*pi*x)-1', 0.3, 1);
%omega_main = Samples(eps,M,'x.^5',0,1);
omega_main = Samples(eps, M, 'sin(3*pi*x)', 0, 1) + Samples(eps, M, 'sin(5*pi*x)', 0.51, 1);

omega = mask.*(omega_main);
%figure; plot(abs(omega)); title('used samples');
%figure; plot(abs(omega_main)); title('all samples');
%% Compute the Generalized function GS_handle and solve for wavelet coefficients with mask
 
disp('Computing GS reconstruction...')

if a==1
    GS_handle = @(x,mode) Haar_Op_Handle(mode, x, S, eps, R, ft_sca, mask);
else
    GS_handle = @(x,mode) CDJV_Op_Handle(mode, x, S, eps, R, a, ft_sca, ft_sca_L, ft_sca_R, mask);
end
opts = spgSetParms('verbosity',1, 'iterations', 500);
wcoeff = spg_bpdn(GS_handle,omega(:),0.0005,opts);


%% Compute image from reconstructed wavelet coefficients

disp('Computing image from reconstructed wavelet coefficients...')
p=max(12, R+1);
if a==1
    
    reconim = get1DHaarReconstruction(wcoeff, p, R);
    
else
    reconim = get1DReconstruction(wcoeff, a, p, R);
end

x_supp=(0:2^-p:1-2^-p)';

orig=f(x_supp);
figure('Position', [150, 150, 1250, 500],'Name','1D reconstruction of a discontinuous funcion from sparse uniform samples for Mask 1');
subplot(2,3,1); plot(x_supp, real(orig)); axis(xylims); title('Original function'); 
subplot(2,3,4); plot(x_supp, real(orig)); xlim([0.4,0.6]); ylim([-1.5,1]);title('Zoom'); 

subplot(2,3,2); plot(x_supp, real(reconim)); axis(xylims); title('GS reconstruction'); 
subplot(2,3,5); plot(x_supp, real(reconim)); xlim([0.4,0.6]); ylim([-1.5,1]);title('Zoom'); 




%% Calulate error of GS reconstruction
% 
error_gs= norm(orig-reconim,2)/norm(orig,2);%sqrt(sum(sum(abs(orig - reconim).^2))/2^(p));
fprintf('GS reconstruction error is %d \n',error_gs);


%% Truncated Fourier representation

disp('Computing trucated Fourier representation...')

omega_Four = mask_centre.*omega_main;

Supp=2^-p:2^-p:1;
fft_rec=zeros(size(Supp));
for j=M:-1:-M+1
    fft_rec=fft_rec+sqrt(eps)*exp(2*eps*1i*j*Supp*pi)*omega_Four(M-j+1);
end

error_fft = norm(fft_rec' - orig,2)/norm(orig,2);%sqrt(sum(abs(f((2^-p:2^-p:1))-fft_rec).^2)/2^(p));
fprintf('Truncated Fourier error is %d \n',error_fft)

subplot(2,3,3); plot(x_supp, real(fft_rec)); axis(xylims); title('Truncated Fourier series'); 
subplot(2,3,6); plot(x_supp, real(fft_rec)); xlim([0.4,0.6]); ylim([-1.5,1]);title('Zoom'); 


load('cww-master/etc/cww_defaults.mat')
figure; plot(x_supp, real(reconim),'color', cww_dflt.red, 'Linewidth', 1.5); axis(xylims); set(gca,'Fontsize',18); %title('CS reconstruction'); 
figure; plot(x_supp, real(fft_rec),'r', 'Linewidth', 1.5); axis(xylims); title('Truncated Fourier series');

