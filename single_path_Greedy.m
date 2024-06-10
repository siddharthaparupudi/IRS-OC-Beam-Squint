
M = 256;    % number of IRS elements
N = 256;    % number of OFDM subcarriers
T = 10;      % the number of time slots
P = 1e-3;   % Total Power at the BS (equal power allocation to all subcarriers)
No = 1e-9;  % Noise power

L1 = 1;     % number of paths in the BS-IRS channel
L2 = 1;     % number of paths in the IRS-user channel

% carrier frequency, bandwidth, wavelength and distance between IRS elements
f_c = 30e9;
W = 400*1e6;
lamda_c = 3e8/f_c;
d = lamda_c/2;

% the subcarrier frequencies
f = linspace(-W/2, W/2, N);

% pathloss exponents
pathloss_BS_IRS = 2;
pathloss_IRS_users = 4;

% the tolerance for the Beam Squint effect (90% threshold)
eps = 0.178/M;

% BS is at (500,0)
% 1st IRS element is at (0,276.725)
% 1024th IRS element is at (0,276.725+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K_set = [1,5,10,50,100,500,1000,5000,10000];


users_x = unifrnd(800,800.1,max(K_set));
users_y = unifrnd(800,800.1,max(K_set));

d_BS_IRS = sqrt((500-0)^2 + (0-276.725)^2);
d_IRS_users = zeros(max(K_set),1);
for k = 1:max(K_set)
    d_IRS_users(k) = sqrt((users_x(k)-0)^2 + (users_y(k)-276.725)^2);
end


% resolvable anglebook of the IRS
anglebook = zeros(M,1);
for i = 1:M
    anglebook(i) = -1+ 2*(i-1)/M;
end

% cascaded normalised angles of the BS-IRS-user channels
psi_C = zeros(max(K_set),1);
for k = 1:max(K_set)
    psi_C(k) = anglebook(randi([1,M]));
end

% normalised angles of the IRS-user channels in presence of beam-squint for all paths, users, frequencies
theta = bsxfun(@times,psi_C,reshape((1+f/f_c), [1,N]));

% the array response of the IRS for all paths, users, frequencies
array_response = ULA_array_2(M,max(K_set),N,theta);

% channel gains
P_alpha = 1e9;
P_beta = 1e6;

% channel gains of the BS-IRS channel
alpha =  sqrt((P_alpha*exp(-1/2))/((d_BS_IRS)^(pathloss_BS_IRS)));


% channel gains of the IRS-user channels
beta = zeros(max(K_set),1);
for k = 1:max(K_set)
    beta(k) = sqrt((P_beta*exp(-1/2))/(d_IRS_users(k))^(pathloss_IRS_users));
end

% channel gains of the cascaded BS-IRS-user channels
gamma_C = zeros(max(K_set),1);
for k = 1:max(K_set)
    gamma_C(k) = alpha*beta(k);
end
% reshape gamma_C (we see same channel gains on all subcarriers)
gamma_C = repmat(gamma_C, [1,N]);


% channel delays (will neglect later as we will apply a synchronization offset)

% max delay offset for nLOS paths
tau_offset = 1e-6;

% channel delays for BS-IRS channel
tau_TR = d_BS_IRS/(3e8);

% channel delays for IRS-user channels
tau_RR = zeros(max(K_set));
for k = 1:max(K_set)
    tau_RR(k) = d_IRS_users(k)/(3e8);
end

% channel delays of the cascaded BS-IRS-user channels
tau_C = zeros(max(K_set));
for k = 1:max(K_set)
    tau_C(k) = tau_TR + tau_RR(k);
end

rates = zeros(length(K_set),1);         % the average rate acheived
max_rates = zeros(length(K_set),1);     % the maximum rate achievable (BF on all subcarriers)

for index = 1:length(K_set)
    K = K_set(index);

    users_x_k = users_x(1:K);
    users_y_k = users_y(1:K);
    d_IRS_users_k = d_IRS_users(1:K);

    psi_C_k = psi_C(1:K);
    theta_k = theta(:,:,1:K,:);
    array_response_k = array_response(:,:,:,1:K,:);
    
    alpha_k = alpha;
    beta_k = beta(1:K);
    gamma_C_k = gamma_C(1:K);
    
    tau_TR_k = tau_TR;
    tau_RR_k = tau_RR(1:K);
    tau_C_k = tau_C(1:K);

    % the set of scheduled users
    schedule = zeros(N,T);

    % Rate on each subcarrier
    Rate = zeros(N,T);

    % total rate in each time slot
    Rate_total = zeros(T,1);

    % the average gain squared
    gain_squared = 0;

    % the average channel gain at each subcarrier
    H_averaged = zeros(N,1);

    % the jain index
    jain_index_slot = zeros(T,1);
    jain_index_slot_gain = zeros(T,1);

    for t = 1:T
        tic
        % generate random phase shifts for the IRS
        a = unifrnd(-1,1);
        phi = 2*pi*a*(0:M-1);

        % the array response vector of the IRS
        array_configuration = zeros(M,1);
        for m = 1:M
            array_configuration(m) = exp(1i*phi(m));
        end
     
        % calculate the channel conditions at each subcarrier
        inner_product = squeeze(sum(array_configuration.*array_response,1));
        H_k = gamma_C.*inner_product;

  
        % at each time slot, schedule the user with best channel conditions on each subcarrier
        for i = 1:N
            % Find the user-subcarrier pair with the highest metric
            [best_gain, best_user] = max(abs(H_k(:,i)));
            
            % Schedule the user on the subcarrier
            schedule(i, t) = best_user;

            % Calculate the rate achieved by the scheduled user on the scheduled subcarrier
            Rate(i,t) = W/N*log2(1 + (P/(No*N))*abs(best_gain)^2);
            
        end
        

        % Calculate the total rate in each time slot (sum of rates across all subcarriers)
        Rate_total(t) = sum(Rate(:,t));
        
        gain_2_slot = 0;
        gain_4_slot = 0;

        for i = 1:N
            gain_squared = gain_squared + abs(H_k(schedule(i,t),i))^2;
            gain_2_slot = gain_2_slot + abs(H_k(schedule(i,t),i))^2;
            gain_4_slot = gain_4_slot + abs(H_k(schedule(i,t),i))^4;
        end

        for i = 1:N
            H_averaged(i) = H_averaged(i) + abs(H_k(schedule(i,t),i))^2/T;
        end

        jain_index_slot(t) = sum(Rate(:,t))^2/(N*sum(Rate(:,t).^2));
        jain_index_slot_gain(t) = gain_2_slot^2/(N*gain_4_slot);
        
        toc
        fprintf('Iteration %d\n', t);
    end

    % Plot magnitude and phase
    figure;
    plot(f, H_averaged);
    xlim([min(f), max(f)]);
    ylim([0, 1.2*max(H_averaged)]);
    title('Magnitude of |H|^2');
    xlabel('Frequency (Hz)');
    ylabel('Average channel gain at each subcarrier'); 

    avg_rate = sum(Rate_total)/T;
    fprintf('Average rate: %f\n', avg_rate);

    d_IRS_UE = min(d_IRS_users);
    rates(index) = avg_rate;
    max_rates(index) = W*log2(1+(P/(No*N))*((M^2*P_alpha*P_beta*(sinc(M*eps))^2)/(exp(1)*d_BS_IRS^(pathloss_BS_IRS)*d_IRS_UE^(pathloss_IRS_users)))*((0.7498*log(K))^(1.71) + 346.474*0.5772));

    gain_squared = gain_squared/(T*N);
    fprintf('Average gain squared on each subcarrier: %f\n', gain_squared);

    jain_index = sum(jain_index_slot)/T;
    fprintf('Jain index: %f\n', jain_index);

    jain_index_gain = sum(jain_index_slot_gain)/T;
    fprintf('Jain index for channel gain: %f\n', jain_index_gain);

    jain_index_avg_channel = sum(H_averaged)^2/(N*sum(H_averaged.^2));
    fprintf('Jain index for average channel gain: %f\n', jain_index_avg_channel);

end

% plot variation of average rate with number of users
figure;
semilogx(K_set, rates,"-o");
hold on;
semilogx(K_set, max_rates,"-*");
xlim([min(K_set), max(K_set)]);
ylim([0, 1.2*max(max_rates)]);
title('Average rate vs. Number of users');
xlabel('Number of users');
ylabel('Average rate (bps)');
legend('Average rate', 'Max rate');

% function to calculate the array response vector of the IRS
% function ULA = ULA_array(M, L1, L2, theta)
%     ULA = zeros(M, L1, L2);
%     factor = -1i*2*pi*(0:M-1).';
%     for m = 1:M
%         ULA(m,:,:) = exp(factor(m).*(theta));
%     end
% end

function ULA = ULA_array_2(M,K,N,theta)
    ULA = zeros(M,K,N);
    factor = -1i*2*pi*(0:M-1).';
    for m = 1:M
        ULA(m,:,:) = exp(factor(m).*(theta));
    end
end
