clear
% 1er test
N = 10000;
K = 1;
theta = 3.3;
param = 1.5;

Y = generer(theta,param,N,K);
histfit(Y,[],'weibull')
moy_cal = mean(Y);
var_cal = var(Y);
moy_exact = theta*gamma(1+1/param);
var_exact = theta^2*gamma(1+2/param) - moy_exact^2;

%2nd test
clear

N = 1000;
K = 500;
theta = 3.3;
param = 1.5;

Y = generer(theta,param,N,K);
plot(Y(:,1))
moy_cal = mean(Y).';
var_cal = var(Y).';
abscisses = (1:K);
% Il faut Plot la ligne qu'il faut
plot(abscisses, moy_cal, abscisses, sqrt(var_cal))

% 2.Estimateur

[alpha_est,BRC] = estimateur_mv(Y,param,N);
disp(['BRC = ', num2str(BRC)])
disp(['moyenne = ', num2str(mean(alpha_est))])
disp(['Variance = ', num2str(var(alpha_est))]) % Reste à comparer avec les valeurs théoriques
disp(['moyenne exatcte = ', num2str(theta^param)])
disp(['Variance = ', num2str(theta^param^2/N)])
% 3. Détection

% 3.1 Théorique
a0 = 0.9;
a1 = 1.5;
Y = generer(theta,param,N,K);
alpha = (0.01:0.01:0.99);
figure
for N=[10,20,50]
    hold on
    p = pi_theorique(a0,a1,2*N, alpha);
    plot(alpha, p)
end
hold off

legend('N=10','N=20','N=50')
N = 20;
figure('Name','a1')
for a1=[1.2,1.5,2]
    hold on
    p = pi_theorique(a0,a1,2*N, alpha);
    plot(alpha, p)
end
hold off

legend('a1=1.2','a1=1.5','a1=2')

%3.2 Pratique
N = 20;
K = 50000;
a0 = 0.9;
a1 = 1.5;
figure('Name','Pratique')
p = pi_estimee(a0,a1,param,N,K,alpha);
plot(alpha,p)


function Y = generer(theta,p,N,K)
    Y = zeros(N,K);
    X = rand(N,K);
    for i=1:N
        for j=1:K
            Y(i,j) = F_inv(X(i,j),theta,p);
        end
    end
end

function y = F_inv(x,theta,p)
    y = theta*(-log(1-x))^(1/p);
end

function [alpha_est,BRC] = estimateur_mv(Y,p,N)
    alpha_est = 1/N*sum(Y.^p);
    BRC = var(alpha_est); % Mettre la vraie valeur
end

function p = pi_theorique(a0, a1, L, alpha)
    lambda = a0/2*chi2inv(1-alpha,L);
    p = 1 - chi2cdf(2*lambda/a1,L);
end

function T = generer_H1(a1,param,N,K)
    theta = a1^(1/param);
    T = generer(theta,param,N,K);
end

function p = pi_estimee(a0,a1,param,N,K,alpha)
    lambda = a0/2*chi2inv(1-alpha,2*N);
    T = generer_H1(a1,param,N,K);
    S = sum(T);
    R = zeros(K,length(alpha));
    for i=1:K
        for j=1:length(alpha)
            R(i,j) = lambda(j) < S(i);
        end
    end
    p = mean(R);
end