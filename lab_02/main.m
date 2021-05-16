% Вариант 17.

x = [14.90, 14.40, 13.56, 15.55, 13.97, 16.33, 14.37, 13.46, 15.51, 14.69,...
     13.41, 14.24, 15.65, 14.54, 13.55, 13.15, 14.32, 15.04, 13.27, 14.60,...
     13.83, 13.93, 14.11, 14.15, 15.48, 15.96, 14.46, 13.87, 13.67, 15.30, ...
     13.95, 16.08, 18.25, 14.93, 15.37, 14.38, 15.56, 13.92, 14.23, 12.80, ...
     13.16, 13.89, 14.24, 13.90, 12.82, 13.20, 13.89, 13.50, 13.44, 16.13, ...
     14.68, 15.27, 13.35, 13.62, 16.16, 16.46, 13.83, 14.13, 15.68, 15.22, ...
     12.59, 12.94, 13.09, 16.54, 14.61, 14.63, 14.17, 13.34, 16.74, 16.30, ...
     13.74, 15.02, 14.96, 15.87, 16.03, 12.87, 14.32, 14.48, 14.57, 14.43, ...
     12.61, 14.52, 15.29, 12.07, 14.58, 11.74, 14.97, 14.31, 12.94, 12.82, ...
     14.13, 14.48, 12.25, 14.39, 15.08, 12.87, 14.25, 15.12, 15.35, 12.27, ...
     14.43, 13.85, 13.16, 16.77, 14.47, 14.89, 14.95, 14.55, 12.80, 15.26, ...
     13.32, 14.92, 13.44, 13.48, 12.81, 15.01, 13.19, 14.68, 14.44, 14.89]; 

gamma = 0.9;
n = length(x);

% 1.a Точечная оценка мат. ожидания.
mu = expectation(x);
fprintf('mu = %.2f\n', mu); 
% 1.a Точечная оценка дисперсии.
sSqr = variance(x);
fprintf('S^2 = %.2f\n\n', sSqr);

% 1.b
% tinv(a, n) - квантиль уровня a распределения Стьюдента с n степенями свободы.
mu_low = mu + (sqrt(sSqr) * tinv((1 - gamma) / 2, n - 1)) / sqrt(n);
mu_high = mu + (sqrt(sSqr) * tinv((1 + gamma) / 2, n - 1)) / sqrt(n);

fprintf('mu_low = %.2f\n', mu_low);
fprintf('mu_high = %.2f\n', mu_high);

% 1.c
% chi2inv(a, n) - квантиль уровня a распределения хи квадрат с n степенями свободы.
sigma_low  = ((n - 1) * sSqr) / chi2inv((1 + gamma) / 2, n - 1);
sigma_high = ((n - 1) * sSqr) / chi2inv((1 - gamma) / 2, n - 1);

fprintf('sigma_low = %.2f\n', sigma_low);
fprintf('sigma_high = %.2f\n', sigma_high);

% 3.a
figure;
% y = mu(x_n) (n = 120)
xs = 1:1:length(x);
ys = mu * ones(n);
plot(xs,ys, 'r'); % blue

hold on
% y = mu(x_n) (n = 1..120)
ys2 = [];
for i = 1:n
    ys2(end + 1) = expectation(x(1:i));
end
plot(xs,ys2, 'y'); % yellow

% y = mu_low(x_n) (n = 1..120)
ys3 = [];
for i = 1:n
    curr_mu = expectation(x(1:i));
    curr_sSqr = variance(x(1:i));
    curr_n = i;
    ys3(end + 1) = curr_mu + (sqrt(curr_sSqr) * tinv((1 - gamma) / 2, curr_n - 1)) / sqrt(curr_n);
end
plot(xs,ys3, 'g');

% y_high = mu(x_n) (n = 1..120)
ys4 = [];
for i = 1:n
    curr_mu = expectation(x(1:i));
    curr_sSqr = variance(x(1:i));
    curr_n = i;
    ys4(end + 1) = curr_mu + (sqrt(curr_sSqr) * tinv((1 + gamma) / 2, curr_n - 1)) / sqrt(curr_n);
end
plot(xs,ys4, 'b');

% 3.b
figure;
% y = sugma(x_n) (n=120)
xs = 1:1:length(x);
ys = sSqr * ones(n);
plot(xs,ys, 'r'); % blue

hold on
% y = sugma(x_n) (n=1..120)
ys2 = [];
for i = 1:n
    ys2(end + 1) = variance(x(1:i));
end
plot(xs,ys2, 'y'); % yellow


% y = sugma_low(x_n) (n=1..120)
ys3 = [];
for i = 1:n
    curr_sSqr = variance(x(1:i));
    curr_n = i;
    ys3(end + 1) = ((curr_n - 1) * curr_sSqr) / chi2inv((1 + gamma) / 2, curr_n - 1);
end
plot(xs,ys3, 'g');

% y = sugma_high(x_n) (n=1..120)
ys4 = [];
for i = 1:n
    curr_sSqr = variance(x(1:i));
    curr_n = i;
    ys4(end + 1) = ((curr_n - 1) * curr_sSqr) / chi2inv((1 - gamma) / 2, curr_n - 1);
end
plot(xs,ys4, 'b');


% Точечная оценка дисперсии.
function sSqr = variance(x)
   sSqr = 0;
   n = length(x);
   mu = expectation(x);
    for i = 1:n
        sSqr = sSqr + (x(i) - mu)^2;
    end
    sSqr = sSqr / (n - 1);
end

% Точечная оценка мат. ожидания.
function mu = expectation(x)
    mu = sum(x) / length(x);
end