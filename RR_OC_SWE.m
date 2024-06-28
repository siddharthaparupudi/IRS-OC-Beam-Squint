
M = 512;    % number of IRS elements
N = 128;    % number of OFDM subcarriers
T = 1;      % the number of time slots
P = 1e-3;   % Total Power at the BS (equal power allocation to all subcarriers)
No = 1e-9;  % Noise power

L1 = 2;     % number of paths in the BS-IRS channel
L2 = 3;     % number of paths in the IRS-user channel

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

% BS is at (500,0)
% 1st IRS element is at (0,276.725)
% 1024th IRS element is at (0,276.725+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K_set = [10];

% resolvable anglebook of the IRS
anglebook = zeros(M,1);
for i = 1:M
    anglebook(i) = -1+ 2*(i-1)/M;
end

users_x = unifrnd(800,800.1,max(K_set));
users_y = unifrnd(800,800.1,max(K_set));

d_BS_IRS = sqrt((500-0)^2 + (0-276.725)^2);
d_IRS_users = zeros(max(K_set),1);
for k = 1:max(K_set)
    d_IRS_users(k) = sqrt((users_x(k)-0)^2 + (users_y(k)-276.725)^2);
end

% cascaded normalised angles of the BS-IRS-user channels
psi_C = zeros(L1,L2,max(K_set));
for k = 1:max(K_set)
    for l1 = 1:L1
        for l2 = 1:L2
            psi_C(l1,l2,k) = anglebook(randi([1,M]));
        end
    end
end

% normalised angles of the IRS-user channels in presence of beam-squint for all paths, users, frequencies
theta = bsxfun(@times,psi_C,reshape((1+f/f_c), [1,1,1,N]));

% the array response of the IRS for all paths, users, frequencies
array_response = ULA_array(M,L1,L2,max(K_set),N,theta);

% channel gains
P_alpha = 1e9;
P_beta = 1e6;

% channel gains of the BS-IRS channel
alpha = zeros(L1,1);
for l1 = 1:L1
    alpha(l1) = sqrt((P_alpha*exprnd(1)*exp(-2*l1))/((d_BS_IRS)^(pathloss_BS_IRS)));
end

% channel gains of the IRS-user channels
beta = zeros(L2,max(K_set));
for k = 1:max(K_set)
    for l2 = 1:L2
        beta(l2,k) = sqrt((P_beta*exprnd(1)*exp(-2*l2))/(d_IRS_users(k))^(pathloss_IRS_users));
    end
end

% channel gains of the cascaded BS-IRS-user channels
gamma_C = zeros(L1,L2,max(K_set));
for k = 1:max(K_set)
    for l1 = 1:L1
        for l2 = 1:L2
            gamma_C(l1,l2,k) = alpha(l1)*beta(l2,k);
        end
    end
end
% reshape gamma_C (we see same channel gains on all subcarriers)
gamma_C = repmat(gamma_C, [1,1,1,N]);


% channel delays

% max delay offset for nLOS paths
tau_offset = 1e-6;

% channel delays for BS-IRS channel
tau_TR = zeros(L1,1);
tau_TR(1) = d_BS_IRS/3e8;
for l1 = 2:L1
    tau_TR(l1) = unifrnd(tau_TR(1), tau_TR(1) + tau_offset);
end

% channel delays for IRS-user channels
tau_RR = zeros(L2,max(K_set));
for k = 1:max(K_set)
    tau_RR(1,k) = d_IRS_users(k)/(3e8);
    for l2 = 2:L2
        tau_RR(l2,k) = unifrnd(tau_RR(1,k), tau_RR(1,k) + tau_offset);
    end
end

% channel delays of the cascaded BS-IRS-user channels
tau_C = zeros(L1,L2,max(K_set));
for k = 1:max(K_set)
    for l1 = 1:L1
        for l2 = 1:L2
            tau_C(l1,l2,k) = tau_TR(l1) + tau_RR(l2,k);
        end
    end
end


rates = zeros(length(K_set),1);

for index = 1:length(K_set)
    K = K_set(index);
        
    users_x_k = users_x(1:K);
    users_y_k = users_y(1:K);

    d_IRS_users_k = d_IRS_users(1:K);

    psi_C_k = psi_C(:,:,1:K);
    theta_k = theta(:,:,1:K,:);
    array_response_k = array_response(:,:,:,1:K,:);
    
    alpha_k = alpha;
    beta_k = beta(:,1:K);
    gamma_C_k = gamma_C(:,:,1:K,:);

    tau_TR_k = tau_TR;
    tau_RR_k = tau_RR(:,1:K);
    tau_C_k = tau_C(:,:,1:K);

    % the set of scheduled users
    schedule = zeros(T);

    % The average rate in the system
    Rate = 0;

    % the average channel gain squared
    average_gain = zeros(T,1);

    % the average channel gain at each subcarrier
    H_averaged = zeros(1,N);

    % jain index
    jain_index_slot_gain = zeros(T,1);

    for t = 1:T
        
        % the phase configuration of the IRS
        phi = zeros(M,1);

        % RR user scheduling
        user = mod(t,K)+1;
        schedule(t) = user;

        % BF phase configuration for the 1st path (i.e LoS path) for the scheduled user
        % which subcarrier to beamform on
        n = N/2; % centre subcarrier
        for i = 1:M
            phi(i) = 2*pi*(i-1)*(psi_C(1,1,user))*(1+ (-W/2 + n*W/N)/f_c);
        end

        % Beamforming IRS configuration when considering all paths (LoS + nLoS) for the scheduled user
        % which subcarrier to beamform on
        % n = N/2;
        % for m = 1:M
        %     c = 0;
        %     for l1 = 1:L1
        %         for l2 = 1:L2
        %             c = c + gamma_C(l1,l2,user)*exp(-1i*2*pi*(tau_C(l1,l2,user))*(-W/2 + n*W/N))*exp(-1i*2*pi*(m-1)*psi_C(l1,l2,user)*(1+(-W/2 + n*W/N)/f_c));
        %         end
        %     end
        %     phi(m) = - angle(c);
        % end

        % the array response vector of the IRS
        array_configuration = zeros(M,1);
        for m = 1:M
            array_configuration(m) = exp(1i*phi(m));
        end

        inner_product = squeeze(sum(array_configuration.*array_response_k,1));
        H_k = squeeze(sum(sum(gamma_C_k.*inner_product.*exp(-1i*2*pi*bsxfun(@times, tau_C_k, reshape(f, [1,1,1,N]))),1),2));


        % add the rate achieved by the scheduled user
        Rate = Rate + sum(W/N*log2(1 + (P/(N*No))*abs(H_k(user,:)).^2));

        % average the channel over time (only add the channels for the scheduled user in each time slot)
        H_averaged = H_averaged + abs(H_k(user,:)).^2./T;

        gain_2_slot = 0;
        gain_4_slot = 0;
        for i = 1:N
            gain_2_slot = gain_2_slot + abs(H_k(user,i))^2;
            gain_4_slot = gain_4_slot + abs(H_k(user,i))^4;
        end

        jain_index_slot_gain(t) = gain_2_slot^2/(N*gain_4_slot);

        centre_subcarrier_gain = abs(H_k(user,N/2))^2;
        sum_gain = 0;
        count = 0;

        for i = 1:N
            subcarrier_gain = abs(H_k(user, i))^2;
            if subcarrier_gain >= 0.9 * centre_subcarrier_gain
                sum_gain = sum_gain + subcarrier_gain;
                count = count + 1;
            end
        end

        average_gain(t) = sum_gain / count;
    
        fprintf('Iteration %d\n', t);
    end

    % Plot magnitude and phase
    figure('Color', 'w');
    plot(f, abs(H_averaged(1,:)));
    title('\textbf{Magnitude of} $\mathbf{|H|^2}$','Interpreter','latex');
    xlabel('\textbf{Frequency (Hz)}','Interpreter','latex');
    ylabel('\textbf{Average channel gain for scheduled user}','Interpreter','latex');
    saveas(gcf, 'Beam_Squint_demo_multipath_2.png');

    % The average rate
    Rate = Rate/T;
    fprintf('The average rate in the system is %f\n', Rate);

    % The average channel gain squared
    average_gain_squared = sum(average_gain)/T;
    fprintf('The average channel gain squared is on centre subcarriers %f\n', average_gain_squared);

    % The Jain's index
    jain_index_gain = sum(jain_index_slot_gain)/T;
    fprintf('The Jain''s index for the gain is %f\n', jain_index_gain);

    jain_average_index = sum(H_averaged)^2/(N*sum(H_averaged.^2));
    fprintf('The Jain''s index for the average gain is %f\n', jain_average_index);

    rates(index) = Rate;

end

figure('Color', 'w');
semilogx(K_set, rates);
title('Variation of Rate with Number of Users');
xlabel('Number of Users (K)', "Interpreter", "latex");
ylabel('Rate (bps)', "Interpreter", "latex");


function ULA = ULA_array(M,L1,L2,K,N,theta)
    ULA = zeros(M,L1,L2,K,N);
    factor = -1i*2*pi*(0:M-1).';
    for m = 1:M
        ULA(m,:,:,:,:) = exp(factor(m).*(theta));
    end
end
