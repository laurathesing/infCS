function [error_Wave,error_Walsh,wc]=Example_handle_1D_flip(R,q,q_L,q2,vm,subsampling_rate,f)

% function to reconstruct f from a flipped sample pattern
% R,q,

load('cww-master/etc/cww_defaults.mat') % load font size, line width, etc.

%% setting

% noise = 1e-3;

dims = 1; % dimension of the signal
disp_plot = 'on'; % show plots
do_save = 0; % possibility to save figures
dest = 'plots'; % location of plots
location = 'northeast';
noise = 1e-3; % noise parameter for spgl reconstruction

if (exist(dest) ~= 7) 
    mkdir(dest);
end


M = 2^(R); 
N = 2^(R+q); % sampling bandwidth
L = 2^(R+q_L); % coefficient bandwidth
N2 = 2^(R+q+q2); % size of evaluation
wname = sprintf('db%d', vm); % wavelet type
j0 = 3;%cww_compute_j0(vm);
nbr_samples = round(N*subsampling_rate) % number of samples

bd_mode = 'bd'; % choice of boundary wavelets

%% precompute values for fast walsh transform

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(R+q_L, dims, wname, bd_mode, j0);

%% samples

samples = cww_walsh_sampling_1d(f,N);

% subsampling pattern and plot

[idx, scales] = sph1_rect2(N, M, nbr_samples, j0);
idx = idx';
idx_plot = zeros(N,1);
idx_plot(idx) =1;

%inverse sampling pattern
idx_plot = fliplr(idx_plot');
figure; imagesc(idx_plot); colormap gray;
idx = (idx_plot>0);
idx_plot2 = zeros(N,1);
idx_plot2(idx) = 1;

%plot new sampling pattern
figure; imagesc(idx_plot2'); colormap gray;

%% optimization setting and evaluation
% spgl operator
G = @(x, mode) cww_handle_1d_cs(x, mode, idx, R+q_L, R+q_L, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces);

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);

% reconstruction
wc = spg_bpdn(G, y, noise, opts_spgl1); 

%% evaluation of the result

residual = norm(y-G(wc,1));
fprintf('Computed solution with error: %g\n', residual )

% compute reconstruction from coefficients
sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

A = cww_get_scaling_matrix(R+q+q2, R+q_L, wname, bd_mode);

x = A*sc; 

% compute original signal
eps = 1e-14;
t1 = linspace(0,1-eps,N2)';

ft1 = f(t1);
ymax = max(max(ft1(:), max(x(:))));
ymin = min(min(ft1(:), min(x(:))));

figure; plot(t1,x, 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);

% compute TW reconstruction
raw_samples = zeros([N2,1]);
raw_samples(1:nbr_samples) = samples(1:nbr_samples);
walsh_approx = fastwht(raw_samples)*N2;
figure; plot(t1, walsh_approx, 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);

% error evaluation - can be changed to other norm or absolute error
error_Walsh = norm(abs(ft1 - walsh_approx),2)/norm(abs(ft1),2);
error_Wave = norm(abs(ft1 - x),2)/norm(abs(ft1),2);

fprintf('Walsh error %g Wavelet error  %g \n',error_Walsh,error_Wave)

% build setting if plots will be saved

set(gca, 'FontSize', cww_dflt.font_size);

r = 0.1;
axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])

if do_save
    fname = sprintf('CS_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

end
