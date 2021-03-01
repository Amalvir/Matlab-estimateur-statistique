clear
% 1er test
N = 10000;
K = 1;
theta = 3.3;
p = 1.5;

Y = generer(theta,p,N,K);
histfit(Y,[],'weibull')
moy_cal = mean(Y);
var_cal = var(Y);
moy_exact = theta*gamma(1+1/p);
var_exact = theta^2*gamma(1+2/p) - moy_exact^2;

%2nd test
clear

N = 1000;
K = 500;
theta = 3.3;
p = 1.5;

Y = generer(theta,p,N,K);
plot(Y(:,1))
moy_cal = mean(Y).';
var_cal = var(Y).';
abscisses = (1:K);
% Il faut Plot la ligne qu'il faut
plot(abscisses, moy_cal, abscisses, sqrt(var_cal))

% 2.Estimateur

[alpha_est,BRC] = estimateur_mv(Y,p,N);
disp(['BRC = ', num2str(BRC)])
disp(['moyenne = ', num2str(mean(alpha_est))])
disp(['Variance = ', num2str(var(alpha_est))]) % Reste à comparer avec les valeurs théoriques
disp(['moyenne exatcte = ', num2str(theta^p)])
disp(['Variance = ', num2str(theta^p^2/N)])
% 3. Détection

% 3.1 Théorique
a0 = 0.9;
a1 = 1.5;
Y = generer(theta,p,N,K);
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

p = pi_estimee(a0,a1,N,K,alpha);
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

function T = generer_H1(a1,N,K)
    theta = a1^(1/p);
    T = generer(theta,p,N,K);
end

function p = pi_estimee(a0,a1,N,K,alpha)
    lambda = a0/2*chi2inv(1-alpha,2*N);
    T = generer_H1(a1,N,K);
    R = sum(T) > lambda;
    p = mean(R);
end