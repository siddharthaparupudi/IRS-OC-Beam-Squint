rng(0);

M = 256;    % number of IRS elements
N = 512;    % number of OFDM subcarriers
T = 100;   % the number of time slots
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

% BS is at (500,0)
% 1st IRS element is at (0,276.725)
% 1024th IRS element is at (0,276.725+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K_set = [1,10,100,1000];


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


rates = zeros(length(K_set),1);

for index = 1:length(K_set)
    K = K_set(index);
        
    users_x_k = users_x(1:K);
    users_y_k = users_y(1:K);
    d_IRS_users_k = d_IRS_users(1:K);

    psi_C_k = psi_C(1:K);
    theta_k = theta(1:K,:);
    array_response_k = array_response(:,1:K,:);

    alpha_k = alpha;
    beta_k = beta(1:K);
    gamma_C_k = gamma_C(1:K,:);
    
    tau_TR_k = tau_TR;
    tau_RR_k = tau_RR(1:K);
    tau_C_k = tau_C(1:K);
    
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
        % which subcarrier to transmit
        n = N/2; % centre subcarrier
        for i = 1:M
            phi(i) = 2*pi*(i-1)*(psi_C_k(user))*(1+ (-W/2 + n*W/N)/f_c);
        end


        % the array response vector of the IRS
        array_configuration = zeros(M,1);
        for m = 1:M
            array_configuration(m) = exp(1i*phi(m));
        end
        
        H_k = zeros(K,N);
    
        % calculate the channel for the scheduled user
        for k = user
            for i = 1:N
                matrix = ULA_array(M, L1, L2, psi_C_k(k)*(1+f(i)/f_c)); 
                inner_product_2 = squeeze(sum(array_configuration .* matrix, 1));
                % H_k(k,i) = sum(sum(gamma_C(k).*inner_product_2.*exp(-1i*2*pi*f(i)*(tau_C(k)))));  
                H_k(k,i) = sum(sum(gamma_C_k(k).*inner_product_2));    
            end
        end


        % % calculate the channel conditions at each subcarrier
        % inner_product = squeeze(sum(array_configuration.*array_response,1));
        % H_k = gamma_C.*inner_product;

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
    figure("Color", 'w');
    plot(f, abs(H_averaged(1,:)));
    xlim([-W/2, W/2]);
    ylim([0, 1.2*max(H_averaged(1,:))]);
    title('\textbf{Magnitude of} $\mathbf{|H|^2}$', 'Interpreter', 'latex');
    xlabel('\textbf{Frequency (Hz)}','Interpreter', 'latex');
    ylabel('\textbf{Average channel gain for scheduled user}','Interpreter', 'latex');  
    saveas(gcf, "Beam_Squint_demo_4.png")

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

% figure;
% semilogx(K_set, rates);
% title('Variation of Rate with Number of Users');
% xlabel('Number of Users');
% ylabel('Rate');

function ULA = ULA_array(M, L1, L2, theta)
    ULA = zeros(M, L1, L2);
    factor = -1i*2*pi*(0:M-1).';
    for m = 1:M
        ULA(m,:,:) = exp(factor(m).*(theta));
    end
end

function ULA = ULA_array_2(M,K,N,theta)
    ULA = zeros(M,K,N);
    factor = -1i*2*pi*(0:M-1).';
    for m = 1:M
        ULA(m,:,:) = exp(factor(m).*(theta));
    end
end


