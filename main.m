%Dans l'énoncé, il y a 2 utilisations differentents de la variable p.
%'param' est donc le paramère de la loi et 'p' les differents pis calculés.

close all
clear
% 1er test

N = 10000;
K = 1;
theta = 3.3;
param = 1.5;

Y = generer(theta,param,N,K);
figure('Name','1. Histogramme')
histfit(Y,[],'weibull')
title('Histogramme de la loi de Weibull et sa densité')
xlabel('Vitesse du vent (m/s)')
ylabel('Densité')
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
figure('Name',"1. Réalisation d'un signal test")

subplot(2,1,1)
plot(Y(:,1))
title('Une réalisation de y')
xlabel('N éléments')
ylabel('Vitesse du vent (m/s)')
moy_cal = mean(Y).';
var_cal = var(Y).';
abscisses = (1:K);

subplot(2,1,2)
plot(abscisses, moy_cal, abscisses, sqrt(var_cal))
title('Moyennes et écart-types des y')
xlabel('K réalisations')
ylabel('(m/s)')
legend('Moyennes','Écart-types')
% 2.Estimateur

[alpha_est,BRC] = estimateur_mv(Y,theta,param,N);
disp(['BRC = ', num2str(BRC)])
disp(['moyenne = ', num2str(mean(alpha_est))])
disp(['Variance = ', num2str(var(alpha_est))])
disp(['moyenne exacte = ', num2str(theta^param)])
disp(['Variance exacte = ', num2str(theta^param^2/N)])

% 3. Détection
% 3.1 Théorique

a0 = 0.9;
a1 = 1.5;
Y = generer(theta,param,N,K);
alpha = (0.01:0.01:0.99);
figure('Name','3.1 pi_théorique avec N varie')
for N=[10,20,50]
    hold on
    p = pi_theorique(a0,a1,2*N, alpha);
    plot(alpha, p)
end
hold off
title('Puissance théorique du test avec différentes valeurs de N')
xlabel('\alpha')
ylabel('\pi_{théorique}')
legend('N=10','N=20','N=50')
N = 20;
figure('Name','3.1 pi_théorique avec a1 varie')
for a1=[1.2,1.5,2]
    hold on
    p = pi_theorique(a0,a1,2*N, alpha);
    plot(alpha, p)
end
hold off
title('Puissance théorique du test avec différentes valeurs de a1')
xlabel('\alpha')
ylabel('\pi_{théorique}')
legend('a1=1.2','a1=1.5','a1=2')

%3.2 Pratique
N = 20;
K = 50000;
a0 = 0.9;
a1 = 1.5;
figure('Name','3.2 Pratique')
p_estim = pi_estimee(a0,a1,param,N,K,alpha);
p_th = pi_theorique(a0,a1,2*N,alpha);
plot(alpha, p_estim, alpha, p_th)
title('Comparaison entre puissance estimée et théorique')
xlabel('\alpha')
ylabel('\pi')
legend('Estimée','Théorique')


function Y = generer(theta,param,N,K)
    % Renvoie Y de taille N x K. K réalisations de N éléments de loi de
    % Weibull W(theta,param).
    
    X = rand(N,K);
    Y = weibullinv(X,theta,param);
end

function y = weibullinv(x,theta,param)
    % Fonction de répartion inverse de la loi de Weibull W(theta,param) en
    % fonction de x réel ou tenseur.
    
    y = theta.*(-log(1-x)).^(1/param);
end

function [alpha_est,BRC] = estimateur_mv(Y,theta,param,N)
    % Renvoie l'estimateur alpha_est des K réalisations de N éléments ainsi
    % que la BCR des estimateurs de a.
    
    alpha_est = 1/N*sum(Y.^param);
    BRC = (theta^param)^2/N; % Mettre la vraie valeur
end

function p = pi_theorique(a0, a1, L, alpha)
    % Renvoie la puissance théorique p du test pour alpha en fonction de
    % a0,a1 et L.
    
    lambda = a0/2*chi2inv(1-alpha,L);
    p = 1 - chi2cdf(2*lambda/a1,L);
end

function T = generer_H1(a1,param,N,K)
    % Génère K réalisations de longueur N associés à l'hypothèse H1
    
    theta = a1^(1/param);
    T = generer(theta,param,N,K);
end

function p = pi_estimee(a0,a1,param,N,K,alpha)
    % % Renvoie la puissance estimée p du test pour alpha en fonction de
    % a0, a1, param, N et K.
    
    lambda = a0/2*chi2inv(1-alpha,2*N);
    Y = generer_H1(a1,param,N,K);
    T = sum(Y);
    R = zeros(K,length(alpha));
    for i=1:K
        for j=1:length(alpha)
            R(i,j) = T(i) > lambda(j);
        end
    end
    p = mean(R); % Moyenne des K séries
end