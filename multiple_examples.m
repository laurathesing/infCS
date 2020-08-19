%% Chose your setting

% uncomment function that will be reconstructed

%f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
%f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + 2*x;
f = @(x) cos(2*pi*x)  + 0.2 * cos(10*pi *x); 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6)+ cos(2*pi*x).*(x <= 0.5); 
%f = @(x) cos(2*pi*x) + cos(10*pi*x).*(x >= 0.5) ;%+ x.*(x < 0.3);


q_L = 6; %wavelet coefficient bandwidth
q2 = 6;  %function evaluation discritization
vm = 4;  % wavelet vanishing moments
iter = 3; %number of examples

R = 5;  % R+q gives the sampling bandwidth

q = [1,2,3]; % length decides about number of examples and N=2^(R+q)

subsampling_rate = 2.^(-q+1); % number of samples nbr = subsampling_rate*N

%% Iteration for different examples

error_Wave = zeros(iter,1); % CS reconstruction error
error_Wal = zeros(iter,1); % TW reconstruction error
wc = zeros(iter,2^(R+q_L)); % recontructed wavelet coefficients


for i = 1:iter
    % computation of the reconstruction the method also internally plot the
    % reconstruction and sampling pattern
    [error_Wave(i),error_Wal(i),wc(i,:)] = Example_handle_1D(R,q(i),q_L,q2,vm,subsampling_rate(i),f);
    % [error_Wave(i),error_Wal(i),wc(i,:)] = Example_handle_1D_flip(R,q(i),q_L,q2,vm,subsampling_rate(i),f);
end

%% evaluation

dims=1; %dimension of signal
wname = sprintf('db%d', vm); 
j0 =3;
bd_mode = 'bd';

% plot wavelet coefficients
t = 1:2^(R+3);
figure; plot(t,abs(wc(1,1:2^(R+3))));%,t,abs(wc(2,1:2^(R+3))),t,abs(wc(3,1:2^(R+3))),'Linewidth',1.1);

% plot CS and TW error for different examples
figure; plot(R+q,error_Wave,R+q,error_Wal,'Linewidth',1.5);
legend('CS error','TW error');