M_set = [1,4,16,64,256,1024];    % number of IRS elements
N_set = [512];    % number of OFDM subcarriers
T = 500;      % the number of time slots
P = 1e-3;   % Total Power at the BS (equal power allocation to all subcarriers)
No = 1e-9;  % Noise power

L1 = 1;     % number of paths in the BS-IRS channel
L2 = 1;     % number of paths in the IRS-user channel

% carrier frequency, bandwidth, wavelength and distance between IRS elements
f_c = 30e9;
W = 400*1e6;
lamda_c = 3e8/f_c;
d = lamda_c/2;

% pathloss exponents
pathloss_BS_IRS = 2;
pathloss_IRS_users = 4;

% the tolerance for the Beam Squint effect (90% threshold)
epsilon = 0.1;
delta = 0;

% BS is at (500,0)
% 1st IRS element is at (0,276.725)
% 1024th IRS element is at (0,276.725+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K_set = [1000];

rates_greedy = zeros(length(K_set),length(M_set), length(N_set));         % the average rate acheived
rates_RR = zeros(length(K_set),length(M_set), length(N_set));             % the average rate acheived by the RR scheme
max_rates_greedy = zeros(length(K_set),length(M_set), length(N_set));     % the maximum rate achievable (BF on all subcarriers)

H_variation_RR = zeros(N_set(end), length(M_set));
H_variation_greedy = zeros(N_set(end), length(M_set));

for index_m = 1:length(M_set)
    for index_n = 1:length(N_set)
        M = M_set(index_m);
        N = N_set(index_n);

        % the subcarrier frequencies
        f = linspace(-W/2, W/2, N);


        users_x = unifrnd(800,900,max(K_set));
        users_y = unifrnd(800,900,max(K_set));

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
            schedule_greedy = zeros(N,T);
            schedule_RR = zeros(T,1);

            % Rate on each subcarrier
            Rate_greedy = zeros(N,T);

            % total rate in each time slot
            Rate_total_greedy = zeros(T,1);
            Rate_total_RR = zeros(T,1);

            % the average gain squared
            gain_squared_greedy = 0;
            gain_squared_RR = 0;
            gain_squared_centre_RR = 0;

            % the average channel gain at each subcarrier
            H_averaged_greedy = zeros(N,1);
            H_averaged_RR = zeros(N,1);

            % the jain index
            jain_index_slot_gain_greedy = zeros(T,1);
            jain_index_slot_gain_RR = zeros(T,1);

            for t = 1:T
                tic
                % generate random phase shifts for the IRS
                a = unifrnd(-1,1);
                phi_greedy = 2*pi*a*(0:M-1);
                % the array response vector of the IRS
                array_configuration_greedy = zeros(M,1);
                for m = 1:M
                    array_configuration_greedy(m) = exp(1i*phi_greedy(m));
                end
                
                if mod(t,K)~=0
                    user_RR = mod(t,K);
                else
                    user_RR = K;
                end
                schedule_RR(t) = user_RR;
                phi_RR = zeros(M,1);
                % BF phase configuration for the 1st path (i.e LoS path) for the scheduled user
                % which subcarrier to transmit
                n = N/2; % centre subcarrier
                for i = 1:M
                    phi_RR(i) = 2*pi*(i-1)*(psi_C_k(user_RR))*(1+ (-W/2 + n*W/N)/f_c);
                end
                % the array response vector of the IRS
                array_configuration_RR = zeros(M,1);
                for m = 1:M
                    array_configuration_RR(m) = exp(1i*phi_RR(m));
                end
            
                % calculate the channel conditions at each subcarrier
                inner_product = squeeze(sum(array_configuration_greedy.*array_response_k,1));
                H_k_greedy = gamma_C_k.*inner_product;

                % calculate the channel conditions at each subcarrier
                inner_product = squeeze(sum(array_configuration_RR.*array_response_k,1));
                H_k_RR = gamma_C_k.*inner_product;


                % at each time slot, schedule the user with best channel conditions on each subcarrier
                for i = 1:N
                    % Find the user-subcarrier pair with the highest metric
                    [best_gain, best_user] = max(abs(H_k_greedy(:,i)));
                    
                    % Schedule the user on the subcarrier
                    schedule_greedy(i, t) = best_user;

                    % Calculate the rate achieved by the scheduled user on the scheduled subcarrier
                    Rate_greedy(i,t) = W/N*log2(1 + (P/(No*N))*abs(best_gain)^2);
                end
                

                % Calculate the total rate in each time slot (sum of rates across all subcarriers)
                Rate_total_greedy(t) = sum(Rate_greedy(:,t));
                % add the rate achieved by the scheduled user
                Rate_total_RR(t) = Rate_total_RR(t) + sum(W/N*log2(1 + (P/(N*No))*abs(H_k_RR(user_RR,:)).^2));
                
                gain_2_slot_greedy = 0;
                gain_4_slot_greedy = 0;
                gain_2_slot_RR = 0;
                gain_4_slot_RR = 0;
            
                for i = 1:N
                    gain_squared_greedy = gain_squared_greedy + abs(H_k_greedy(schedule_greedy(i,t),i))^2;
                    gain_squared_RR = gain_squared_RR + abs(H_k_RR(user_RR,i))^2;
                    gain_2_slot_greedy = gain_2_slot_greedy + abs(H_k_greedy(schedule_greedy(i,t),i))^2;
                    gain_4_slot_greedy = gain_4_slot_greedy + abs(H_k_greedy(schedule_greedy(i,t),i))^4;
                    gain_2_slot_RR = gain_2_slot_RR + abs(H_k_RR(user_RR,i))^2;
                    gain_4_slot_RR = gain_4_slot_RR + abs(H_k_RR(user_RR,i))^4;
                end

                for i = 1:N
                    H_averaged_greedy(i) = H_averaged_greedy(i) + abs(H_k_greedy(schedule_greedy(i,t),i))^2/T;
                    H_averaged_RR(i) = H_averaged_RR(i) + abs(H_k_RR(user_RR,i))^2/T;
                end

                H_variation_RR(:,index_m) = H_averaged_RR;
                H_variation_greedy(:,index_m) = H_averaged_greedy;

                centre_subcarrier_gain = abs(H_k_RR(user_RR,N/2))^2;
                sum_gain = 0;
                count = 0;

                for i = 1:N
                    subcarrier_gain = abs(H_k_RR(user_RR, i))^2;
                    if subcarrier_gain >= 0.9 * centre_subcarrier_gain
                        sum_gain = sum_gain + subcarrier_gain;
                        count = count + 1;
                    end
                end

                gain_squared_centre_RR = gain_squared_centre_RR + sum_gain/(count*T); 

                jain_index_slot_gain_RR(t) = gain_2_slot_RR^2/(N*gain_4_slot_RR);
                jain_index_slot_gain_greedy(t) = gain_2_slot_greedy^2/(N*gain_4_slot_greedy);
                
                toc
                fprintf('Iteration %d\n', t);
            end

            % Plot magnitude and phase
            % figure("Color", 'w');
            % plot(f, H_averaged_RR);
            % hold on;
            % plot(f, H_averaged_greedy);
            % xlim([min(f), max(f)]);
            % ylim([0, 1.2*max(max(H_averaged_greedy), max(H_averaged_RR))]);
            % title('\textbf{Magnitude of} $\mathbf{|H|^2}$', 'Interpreter', 'latex');
            % xlabel('\textbf{Frequency (Hz)}', 'Interpreter', 'latex'); 
            % ylabel('\textbf{Average channel gain at each subcarrier}', 'Interpreter', 'latex'); 
            % legend('Average channel in RR', 'Average channel in Greedy', 'Interpreter', 'latex')

            avg_rate_greedy = sum(Rate_total_greedy)/T;
            avg_rate_RR = sum(Rate_total_RR)/T;
            fprintf('Average rate RR: %f\n', avg_rate_RR/W);
            fprintf('Average rate Greedy: %f\n', avg_rate_greedy/W); 

            d_IRS_UE = mean(d_IRS_users);
            rates_greedy(index, index_m, index_n) = (1/W)*avg_rate_greedy;
            rates_RR(index, index_m, index_n) = (1/W)*avg_rate_RR;
            max_rates_greedy(index, index_m, index_n) = (1/W)*(1-delta)*W*log2(1+(P/(No*N))*((M^2*P_alpha*P_beta*(1-epsilon))/(exp(1)*d_BS_IRS^(pathloss_BS_IRS)*d_IRS_UE^(pathloss_IRS_users)))*((0.7498*log(K))^(1.71) + pi));
         
            gain_squared_greedy = gain_squared_greedy/(T*N);
            gain_squared_RR = gain_squared_RR/(T*N);
            fprintf('Average gain squared on each subcarrier in RR: %f\n', gain_squared_RR);
            fprintf('Average gain squared on centre subcarriers in RR: %f\n', gain_squared_centre_RR);
            fprintf('Average gain squared on each subcarrier in Greedy: %f\n', gain_squared_greedy);

            jain_index_gain_RR = sum(jain_index_slot_gain_RR)/T;
            jain_index_gain_greedy = sum(jain_index_slot_gain_greedy)/T;
            fprintf('Jain index for channel gain in RR: %f\n', jain_index_gain_RR);
            fprintf('Jain index for channel gain in Greedy: %f\n', jain_index_gain_greedy);

            jain_index_avg_channel_RR = sum(H_averaged_RR)^2/(N*sum(H_averaged_RR.^2));
            jain_index_avg_channel_greedy = sum(H_averaged_greedy)^2/(N*sum(H_averaged_greedy.^2));
            fprintf('Jain index for average channel gain in RR: %f\n', jain_index_avg_channel_RR);
            fprintf('Jain index for average channel gain in Greedy: %f\n', jain_index_avg_channel_greedy);
        end
    end
end


% figure('Color', 'w');
% legendLabels = {};
% for index_m = 1:length(M_set)
%     hold on;
%     plot(f, H_variation_RR(:, index_m));
%     legendLabels{end+1} = sprintf('Average channel gain RR (M = %d, N = %d)', M_set(index_m), N_set(end));
%     hold on;
%     plot(f, H_variation_greedy(:, index_m));
%     legendLabels{end+1} = sprintf('Average channel gain Greedy (M = %d, N = %d)', M_set(index_m), N_set(end));
%     xlim([min(f), max(f)]);
%     ylim([0, 1.2*max(max(H_variation_greedy(:, index_m)), max(H_variation_RR(:, index_m)))]);
%     title('\textbf{Magnitude of} $\mathbf{|H|^2}$', 'Interpreter', 'latex');
%     xlabel('\textbf{Frequency (Hz)}', 'Interpreter', 'latex');
%     ylabel('\textbf{Average channel gain at each subcarrier}', 'Interpreter', 'latex');
% end
% legend(legendLabels, 'Interpreter', 'latex'); 


% % plot variation of average rate with number of users
% figure('Color', 'w');
% set(gca, 'XScale', 'log');
% legendLabels = {};
% for index_m = 1:length(M_set)
%     for index_n = 1:length(N_set)
%         hold on;
%         semilogx(K_set, rates_RR(:, index_m, index_n), "-x");
%         legendLabels{end+1} = sprintf('Average rate RR (M = %d, N = %d)', M_set(index_m), N_set(index_n));
%         hold on;
%         semilogx(K_set, rates_greedy(:, index_m, index_n),"-o");
%         legendLabels{end+1} = sprintf('Average rate Greedy (M = %d, N = %d)', M_set(index_m), N_set(index_n));
%         hold on;
%         semilogx(K_set, max_rates_greedy(:, index_m, index_n),"-*");
%         legendLabels{end+1} = sprintf('Max rate (M = %d, N = %d)', M_set(index_m), N_set(index_n));
%         xlim([min(K_set), max(K_set)]);
%         ylim([0, 1.2*max(max(max(max_rates_greedy(:, index_m, index_n)), max(rates_RR(:, index_m, index_n))), max(rates_greedy(:, index_m, index_n)))]);
%         title('\textbf{Average rate vs. Number of users}', 'Interpreter', 'latex');
%         xlabel('\textbf{Number of users (K)}', 'Interpreter', 'latex');
%         ylabel('\textbf{Average rate (bps/Hz)}', 'Interpreter', 'latex');
%     end
% end
% legend(legendLabels, 'Interpreter', 'latex'); 

figure('Color', 'w');
for index_n = 1:length(N_set)
    hold on;
    semilogx(log2(M_set), rates_RR(length(K_set), :, index_n), "-x");
    hold on;
    semilogx(log2(M_set), rates_greedy(length(K_set), :, index_n), "-o");
    hold on;
    semilogx(log2(M_set), max_rates_greedy(length(K_set), :, index_n), "-*");
    xlim([min(log2(M_set)), max(log2(M_set))]);
    ylim([0, 1.2*max(max(max(max_rates_greedy(length(K_set), :, index_n)), max(rates_RR(length(K_set), :, index_n))), max(rates_greedy(length(K_set), :, index_n)))]);
    xt = get(gca, 'XTick');
    set (gca, 'XTickLabel', 2.^xt);
    title('\textbf{Average rate vs. Number of IRS elements}', 'Interpreter', 'latex');
    xlabel('\textbf{Number of IRS elements (M)}', 'Interpreter', 'latex');
    ylabel('\textbf{Average rate (bps/Hz)}', 'Interpreter', 'latex');
    legend('Average rate RR', 'Average rate Greedy', 'Max rate', 'Interpreter', 'latex');
end



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
