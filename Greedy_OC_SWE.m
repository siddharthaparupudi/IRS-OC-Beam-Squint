
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
K_set = [1,10,50,100,250,500,750,1000];

rates = zeros(length(K_set),1);
max_rates = zeros(length(K_set),1);

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
    L1 = 1;
    L2_K = randi([1],K,1);


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
    P_alpha = 1e6;
    P_beta = 1e3;

    % channel gains of the BS-IRS channel
    alpha = zeros(L1,1);
    for l1 = 1:L1
        alpha(l1) = sqrt((P_alpha*exprnd(1)*exp(-l1/2))/((d_BS_IRS)^(pathloss_BS_IRS)));
    end

    % channel gains of the IRS-user channels
    beta = cell(K,1);
    for k = 1:K
        beta{k} = zeros(L2_K(k),1);
        for l2 = 1:L2_K(k)
            beta{k}(l2) = sqrt((P_beta*exprnd(1)*exp(-l2/2))/(d_IRS_users(k))^(pathloss_IRS_users));
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
    tau_TR(1) = 0;
    for l1 = 2:L1
        tau_TR(l1) = unifrnd(tau_TR(1), tau_TR(1) + tau_offset);
    end

    % channel delays for IRS-user channels
    tau_RR = cell(K,1);
    for k = 1:K
        tau_RR{k} = zeros(L2_K(k),1);
        tau_RR{k}(1) = 0;
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
    N = 128;

    % the number of time slots
    T = 3;

    % the set of scheduled users
    schedule = zeros(N,T);

    % Rate on each subcarrier
    Rate = zeros(N,T);

    % total rate in each time slot
    Rate_total = zeros(T,1);

    % Total Power in the system (equal power allocation to all subcarriers)
    P = 1e-3;

    % Noise power
    No = 1e-9;

    % the average gain squared
    gain_squared = 0;

    H_averaged = zeros(N,1);

    % the jain index
    jain_index_slot = zeros(T,1);
    jain_index_slot_gain = zeros(T,1);

    for t = 1:T
        tic

        % the phase configuration of the IRS
        a = unifrnd(-1,1);
        phi = 4*pi*a*(0:M-1);

        % the array response vector of the IRS
        array_vector = zeros(M,1);
        for m = 1:M
            array_vector(m) = exp(1i*phi(m));
        end
     
        % calculate the channel conditions at each subcarrier
        f = linspace(f_c-W/2, f_c+W/2, N);

        H_k = zeros(K,length(f));

        % channel calculation with for loops
        % for k = 1:K
        %     for i = 1:length(f)
        %         h_k_matrix = zeros(L1,L2_K(k));
        %         for m = 1:M
        %              h_k_matrix = h_k_matrix + exp(1i*phi(m))*gamma_C{k}.*exp(-1i*2*pi*(m-1)*psi_C{k}*(1+f(i)/f_c)).*exp(-1i*2*pi*f(i)*(tau_C{k}));
        %         end
        %         H_k(k,i) = sum(sum(h_k_matrix)); 
        %     end
        % end

        % channel calculation with matrix multiplication
        for k = 1:K
            for i = 1:length(f)
                matrix = ULA_array(M, L1, L2_K(k), psi_C{k}*(1+f(i)/f_c)); 
                inner_product = squeeze(sum(array_vector .* matrix, 1));
                H_k(k,i) = sum(sum(gamma_C{k}.*inner_product.*exp(-1i*2*pi*f(i)*(tau_C{k}))));      
            end
        end

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
    title('Magnitude of |H|^2');
    xlabel('Frequency (Hz)');
    ylabel('Average channel gain at each subcarrier'); 


    avg_rate = sum(Rate_total)/T;
    fprintf('Average rate: %f\n', avg_rate);

    d_IRS_UE = min(d_IRS_users);
    rates(index) = avg_rate;
    max_rates(index) = W*log2(1+(P/(No*N))*((M^2*P_alpha*P_beta)/(exp(1)*d_BS_IRS^2*d_IRS_UE^4))*((0.7498*log(K))^(1.71) + 346.474*0.5772));

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
semilogx(K_set, rates);
hold on;
semilogx(K_set, max_rates);
title('Average rate vs. Number of users');
xlabel('Number of users');
ylabel('Average rate (bps)');
legend('Average rate', 'Max rate');

% function to calculate the array response vector of the IRS
function ULA = ULA_array(M, L1, L2, theta)
    ULA = zeros(M, L1, L2);
    factor = -1i*2*pi*(0:M-1).';
    for m = 1:M
        ULA(m,:,:) = exp(factor(m).*(theta));
    end
end

