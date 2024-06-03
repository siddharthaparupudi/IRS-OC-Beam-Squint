
% number of IRS elements
M = 1024;

% carrier frequency, bandwidth, wavelength and distance between IRS elements
f_c = 30e9;
W = 400*1e6;
lamda_c = 3e8/f_c;
d = lamda_c/2;

% pathloss exponents
pathloss_BS_IRS = 2;
pathloss_IRS_users = 4;

% BS is at (500,0)
% 1st IRS element is at (0,276.725)
% 1024th IRS element is at (0,276.725+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K_set = [1000];

rates = zeros(length(K_set),1);

for index = 1:length(K_set)
    K = K_set(index);
        
    users_x = unifrnd(800,900,K);
    users_y = unifrnd(800,900,K);

    d_BS_IRS = sqrt((500-0)^2 + (0-276.725)^2);
    d_IRS_users = zeros(K,1);
    for k = 1:K
        d_IRS_users(k) = sqrt((users_x(k)-0)^2 + (users_y(k)-276.725)^2);
    end


    % number of paths
    % L1 = number of paths from BS to IRS
    % L2_K(k) = number of paths from IRS to kth user
    L1 = 2;
    L2_K = randi([2,4],K,1);


    % resolvable anglebook of the IRS
    anglebook = zeros(M,1);
    for i = 1:M
        anglebook(i) = -1+ 2*(i-1)/M;
    end

    % cascaded normalised angles of the BS-IRS-user channels
    psi_C = cell(K,1);
    for k = 1:K
        psi_C{k} = zeros(L1,L2_K(k));
        for l1 = 1:L1
            for l2 = 1:L2_K(k)
                psi_C{k}(l1,l2) = anglebook(randi([1,M]));
            end
        end
    end

    % channel gains

    % channel gains of the BS-IRS channel
    alpha = zeros(L1,1);
    for l1 = 1:L1
        alpha(l1) = 1e6*sqrt(exprnd(1))/(sqrt(d_BS_IRS))^(pathloss_BS_IRS)*exp(-l1/2);
    end

    % channel gains of the IRS-user channels
    beta = cell(K,1);
    for k = 1:K
        beta{k} = zeros(L2_K(k),1);
        for l2 = 1:L2_K(k)
            beta{k}(l2) = 1e3*sqrt(exprnd(1))/(sqrt(d_IRS_users(k)))^(pathloss_IRS_users)*exp(-l2/2);
        end
    end


    % channel gains of the cascaded BS-IRS-user channels
    gamma_C = cell(K,1);
    for k = 1:K
        gamma_C{k} = zeros(L1,L2_K(k));
        for l1 = 1:L1
            for l2 = 1:L2_K(k)
                gamma_C{k}(l1,l2) = alpha(l1)*beta{k}(l2);
            end
        end
    end


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
    tau_RR = cell(K,1);
    for k = 1:K
        tau_RR{k} = zeros(L2_K(k),1);
        tau_RR{k}(1) = d_IRS_users(k)/3e8;
        for l2 = 2:L2_K(k)
            tau_RR{k}(l2) = unifrnd(tau_RR{k}(1), tau_RR{k}(1) + tau_offset);
        end
    end

    % channel delays of the cascaded BS-IRS-user channels
    tau_C = cell(K,1);
    for k = 1:K
        tau_C{k} = zeros(L1,L2_K(k));
        for l1 = 1:L1
            for l2 = 1:L2_K(k)
                tau_C{k}(l1,l2) = tau_TR(l1) + tau_RR{k}(l2);
            end
        end
    end


    % number of OFDM subcarriers
    N = 64;

    % the number of time slots
    T = 500;

    % the set of scheduled users
    schedule = zeros(T);

    % The average rate in the system
    Rate = 0;

    % Total power in the system (equal power allocation to all subcarriers)
    P = 1e-3;

    % Noise power
    No = 1e-9;

    % the average channel gain squared
    average_gain = zeros(T,1);

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
        % n = N/2; % centre subcarrier
        % for i = 1:M
        %     phi(i) = 4*pi*(i-1)*(psi_C{user}(1,1)) + 2*pi*(i-1)*(n*W/N - W/2)*(psi_C{user}(1,1))/f_c;
        % end

        % Beamforming IRS configuration when considering all paths (LoS + nLoS) for the scheduled user
        for m = 1:M
            c = 0;
            for l1 = 1:L1
                for l2 = 1:L2_K(user)
                    c = c + gamma_C{user}(l1,l2)*exp(-1i*2*pi*(tau_C{user}(l1,l2))*f_c)*exp(-1i*2*pi*(m-1)*psi_C{user}(l1,l2)*(1+f_c/f_c));
                end
            end
            phi(m) = - angle(c);
        end

        % the array response vector of the IRS
        array_vector = zeros(M,1);
        for m = 1:M
            array_vector(m) = exp(1i*phi(m));
        end
        

        % calculate the channel conditions at each subcarrier
        f = linspace(f_c-W/2, f_c+W/2, N);

        H_k = zeros(K,length(f));
    
        % calculate the channel for the scheduled user
        for k = user
            for i = 1:length(f)
                matrix = ULA_array(M, L1, L2_K(k), psi_C{k}*(1+f(i)/f_c)); 
                inner_product = squeeze(sum(array_vector .* matrix, 1));
                H_k(k,i) = sum(sum(gamma_C{k}.*inner_product.*exp(-1i*2*pi*f(i)*(tau_C{k}))));      
            end
        end


        % add the rate achieved by the scheduled user
        Rate = sum(W/N*log2(1 + (P/(N*No))*abs(H_k(user,:)).^2));

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
    figure;
    plot(f, abs(H_averaged(1,:)));
    title('Magnitude of |H|^2');
    xlabel('Frequency (Hz)');
    ylabel('Average channel gain for scheduled user');



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