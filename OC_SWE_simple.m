
% path gains
alpha_1 = 1.75;
beta_1 = 1;

% number of IRS elements
M = 1024;

% path delays
tau_1 = 1e-6;
tau_2 = 5*1e-6;

% carrier frequency, bandwidth, wavelength and distance between IRS elements
f_c = 30e9;
W = 400*1e6;
lamda_c = 3e8/f_c;
d = lamda_c/2;

% angles of arrival and departure
DoA = pi/3;
DoD = 0;

% normalized angles of arrival and departure
psi_1_TR = d*sin(DoA)/lamda_c;
psi_1_RR = d*sin(DoD)/lamda_c;

% effective angles and gains of the cascaded BS-IRS-UE channel
psi_1_C = psi_1_TR - psi_1_RR;
gamma_C = alpha_1*beta_1;
%  dont need to calculate effective delay as we will synchronize the signals
tau_C = tau_1 + tau_2;

% set the IRS phase configuration

% align to f_c
% phi = zeros(M,1);
% for i = 1:M
%     phi(i) = 4*pi*(i-1)*(psi_1_C);
% end

% number of OFDM subcarriers
N = 128;

% the number of time slots
T = 10;

for t = 1:T
    
    % set random phases for Opportunistic communications
    phi = zeros(M,1);

    % a = unifrnd(0,1);
    % for i = 1:M
    %     phi(i) = 2*pi*a*(i-1);
    % end

    % which subcarrier to transmit
    n = 64;
    for i = 1:M
        phi(i) = 4*pi*(i-1)*(psi_1_C) + 2*pi*(i-1)*(n*W/N - W/2)*(psi_1_C)/f_c;
    end


    f = linspace(f_c-W/2, f_c+W/2, 1000);

    H = zeros(size(f));
    H_averaged = zeros(size(f));
    for i = 1:length(f)
        h = 0;
        for m = 1:M
            h = h + gamma_C*exp(1i*phi(m))*exp(-1i*2*pi*(m-1)*psi_1_C*(1+f(i)/f_c));
        end 
        H(i) = h;
    end

    H_averaged = H_averaged + H;
    fprintf('Iteration %d\n', t);
end

% Plot magnitude and phase
figure;
plot(f, abs(H_averaged));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


