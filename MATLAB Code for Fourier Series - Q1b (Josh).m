% Fourier sine series 
% Series terms have form B_k sin(lambda_k x)
% Dan Coombs, Sep 17 2020
% Adapted by Joshua Santo, Oct. 19, 2021
% for Group Assignment 2, Q1(b)

clear all; 
close all;

L=1; % set length of domain
pts=11; % set number of points to be plotted
M=100; % number of terms to be used in the sum
t=1; %timestamp (do this for 0, 0.05, 0.1, 0.2, and 1)

B = @(k) (2/((k^3)*(pi^3)))*(1-(-1)^k); 
%set B_k terms as function of k
lambda = @(k) pi * k / L; 
% set lambda_k terms similarly

x=(0:L/pts:1); % set up x values
u=zeros(1,length(x)); % zero vector of u values

% this loop sums the series from k=1 to k=M
for k=1:M
    u = u + B(k) * sin(lambda(k)*x)*exp(-(k^2)*(pi^2)*t); 
end

figure(1);


%finally plot y against x:
plot(x,u,'LineWidth',2.0)
xlabel('Length')
ylabel('Concentration')
title('Q1(b) - Fourier Series Solution Concentration Distribution at t=1')

