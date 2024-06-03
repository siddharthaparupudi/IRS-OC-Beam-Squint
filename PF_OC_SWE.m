% number of IRS elements
M = 1024;

% carrier frequency, bandwidth, wavelength and distance between IRS elements
f_c = 30e9;
W = 400*1e6;
lamda_c = 3e8/f_c;
d = lamda_c/2;

% pathloss exponent/s
pathloss_BS_IRS = 2;
pathloss_IRS_users = 4;

% BS is at (500,0)
% 1st IRS element is at (0,200)
% 1024th IRS element is at (0,200+1023*d)
% users are randomly distributed in the rectangle (800,800), (800,900), (900,800), (900,900)
% K users
K = 12;
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


% angles of arrival and departure
DoAs = unifrnd(atan(200/500),atan(353.45/500),L1,1);
DoDs = cell(K,1);
for k = 1:K
    DoDs{k} = unifrnd(atan(523.725/850), atan(623.725/850), L2_K(k),1);
end

% normalised angles of arrival and departure
psi_1_TR = d*sin(DoAs)/lamda_c;
psi_2_RR = cell(K,1);
for k = 1:K
    psi_2_RR{k} = d*sin(DoDs{k})/lamda_c;
end

% cascaded normalised angles of the BS-IRS-user channels
% psi_C = cell(K,1);
% for k = 1:K
%     psi_C{k} = zeros(L1,L2_K(k));
%     for l1 = 1:L1
%         for l2 = 1:L2_K(k)
%             psi_C{k}(l1,l2) = psi_1_TR(l1) - psi_2_RR{k}(l2);
%         end
%     end
% end

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
    alpha(l1) = 1e3/(sqrt(d_BS_IRS))^(pathloss_BS_IRS);
end

% channel gains of the IRS-user channels
beta = cell(K,1);
for k = 1:K
    beta{k} = zeros(L2_K(k),1);
    for l2 = 1:L2_K(k)
        beta{k}(l2) = 1e6/(sqrt(d_IRS_users(k)))^(pathloss_IRS_users);
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
N = 128;

% the number of time slots
T = 100;

% the acceptable delay in the system to reach the long term gain
tau = 10;

% the set of scheduled users
schedule = zeros(N,T);

% the long term gains of each user
T_k = ones(K,1);


for t = 1:T
    
    % the phase configuration of the IRS
    phi = zeros(M,1);

    % set random phases for opportunistic communications at each time slot
    % a = unifrnd(-1,1);
    % for i = 1:M
    %     phi(i) = 4*pi*a*(i-1);
    % end

    % BF phase configuration for the 1st path (i.e LoS path) with 
    % the subcarrier chosen for BF
    n = 64;
    for i = 1:M
        phi(i) = 4*pi*(i-1)*(psi_C{1}(1,1)) + 2*pi*(i-1)*(n*W/N - W/2)*(psi_C{1}(1,1))/f_c;
    end

    % choose phase configuration to make h(f) real at the n'th subcarrier for 1'st user
    % user = 1;
    % for m = 1:M
    %     c = 0;
    %     for l1 = 1:L1
    %         for l2 = 1:L2_K(user)
    %             c = c + gamma_C{user}(l1,l2)*exp(-1i*2*pi*(tau_TR(l1) + tau_RR{user}(l2))*f_c)*exp(-1i*2*pi*(m-1)*psi_C{user}(l1,l2)*(1+f_c/f_c));
    %         end
    %     end
    %     phi(m) = - angle(c);
    % end

    % calculate the channel conditions at each subcarrier
    f = linspace(f_c-W/2, f_c+W/2, N);

    H_k = zeros(K,length(f));
    H_averaged = zeros(K,length(f));
    for k = 1:K
        for i = 1:length(f)
            h_k = 0;
            for m = 1:M
                for l1 = 1:L1
                    for l2 = 1:L2_K(k)
                        % due to snchronization offset applied at the receiver,
                        % we don't need the term *exp(-1i*2*pi*f_c*tau_C{k}(1,1))
                        h_k = h_k + gamma_C{k}(l1,l2)*exp(1i*phi(m))*exp(-1i*2*pi*(m-1)*psi_C{k}(l1,l2)*(1+f(i)/f_c))*exp(-1i*2*pi*f(i)*(tau_TR(l1) + tau_RR{k}(l2)));
                    end 
                end
            end
            H_k(k,i) = h_k;
        end
    end

    % at each time slot, at each subcarrier,
    % calculate the metric abs(H_k(f,t))/T_k(t) for each user
    % schedule the user with the highest metric at each subcarrier
    % If a user is the best user for more than one subcarrier,
    % then the user is scheduled for the subcarrier with the highest metric
    % If a user is the best user for more than one subcarrier with the same metric,
    % then the user is scheduled for the subcarrier with the lowest index
    % Do until all subcarriers are scheduled or we have scheduled all users

    % If the user is scheduled at subcarrier i, then the user's long term gain is updated as:
    % T_k(t+1) = (1-1/tau)T_k(t) + (1/tau)abs(H_k(f_i, t)) 
    % If the user is not scheduled at any subcarrier, then the user's long term gain is updated as:
    % T_k(t+1) = (1-1/tau)T_k(t)

    PF_metrics = zeros(K,N);
    for i = 1:N
        PF_metrics(:,i) = abs(H_k(:,i))./T_k;
    end

    while nnz(schedule(:,t)) < min(N,K)
        % Find the user-subcarrier pair with the highest metric
        [best_metric, idx] = max(PF_metrics(:));
        [best_user, best_subcarrier] = ind2sub(size(PF_metrics), idx);
        
        % Schedule the user on the subcarrier
        schedule(best_subcarrier, t) = best_user;
        fprintf('User %d is scheduled on subcarrier %d at time slot %d\n' , best_user, best_subcarrier, t);
        fprintf('%d ', schedule(:,t));
        fprintf('\n');
        % Update the user's long term gain
        T_k(best_user) = (1 - 1/tau) * T_k(best_user) + (1/tau) * abs(H_k(best_user,best_subcarrier));
        
        % Remove all other subcarriers for this user from consideration
        PF_metrics(best_user, :) = -inf;
        % Remove all other users for this subcarrier from consideration
        PF_metrics(:, best_subcarrier) = -inf;
    end

    % update the long term gains of the users that are not scheduled
    not_scheduled = setdiff(1:K, schedule(:,t));
    for i = 1:length(not_scheduled)
        T_k(not_scheduled(i)) = (1 - 1/tau) * T_k(not_scheduled(i));
    end


    H_averaged = H_averaged + abs(H_k);
   
    fprintf('Iteration %d\n', t);
end

% Plot magnitude and phase
figure;
subplot(3,2,1);
plot(f, abs(H_averaged(1,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 1');


subplot(3,2,2);
plot(f, abs(H_averaged(2,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 2');


subplot(3,2,3);
plot(f, abs(H_averaged(3,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 3');


subplot(3,2,4);
plot(f, abs(H_averaged(4,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 4');

subplot(3,2,5);
plot(f, abs(H_averaged(5,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 5');


subplot(3,2,6);
plot(f, abs(H_averaged(6,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 6');

figure;
subplot(3,2,1);
plot(f, abs(H_averaged(7,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 7');


subplot(3,2,2);
plot(f, abs(H_averaged(8,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 8');


subplot(3,2,3);
plot(f, abs(H_averaged(9,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 9');


subplot(3,2,4);
plot(f, abs(H_averaged(10,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 10');

subplot(3,2,5);
plot(f, abs(H_averaged(11,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 11');


subplot(3,2,6);
plot(f, abs(H_averaged(12,:)));
title('Magnitude of H');
xlabel('Frequency (Hz)');
ylabel('Magnitude of user 12');

T_k