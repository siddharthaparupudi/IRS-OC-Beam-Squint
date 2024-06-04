
z = linspace(0,10,1000)';
% the CCDF of Z (bessel function of the second kind)
bessel = 2*sqrt(z).*besselk(1,2*sqrt(z));

% Remove NaN values
validData = ~isnan(bessel);
z = z(validData);
bessel = bessel(validData);

% the CDF of Z
F = 1-bessel;

% Fit the CCDF of Z to an exponential function
expfit = fittype('(exp(b*x))', 'independent', 'x', 'coefficients', {'b'});
f1 = fit(z, bessel, expfit, 'StartPoint', [1]);

disp(f1)

% Plot the CCDF of Z and the fitted curve
figure;
plot(z, bessel)
hold on;
plot(f1)
xlabel('z')
ylabel('1-F(z)')
legend('Original Data','Fitted Curve')
title('CCDF of F(z) and fitted curve')

% compute the inverse function of the CDF of Z
F_inverse = @(x) interp1(F, z, x);

% Generate the number of users
users = linspace(1, 10000, 10000);

% l_K = F_inverse(1-1/users);
l_K = F_inverse(1-1./users);

% Remove NaN values
validData = ~isnan(l_K);
users = users(validData);
l_K = l_K(validData);

% Fit the inverse function to a logarithmic function
logfit = fittype('(b*log(x))^c', 'independent', 'x', 'coefficients', {'b','c'});
f2 = fit(users', l_K', logfit, 'StartPoint', [1, 1]);

disp(f2)

% Plot the inverse function
figure;
plot(users, l_K);
hold on;
plot(f2);
xlabel('Number of users (K)');
ylabel('l_K');
legend('Original Data','Fitted Curve')
title('Gain contribution due to multiuser diversity');