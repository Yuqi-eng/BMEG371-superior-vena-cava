% Fourier sine series 
% Series terms have form B_k sin(lambda_k x)
% Dan Coombs, Sep 17 2020
% Adapted by Joshua Santo, Oct. 19, 2021
% for Group Assignment 2, Q1(b)

clear all; 
close all;

L=1; % set length of domain
pts=10; % set number of points to be plotted
M=100; % number of terms to be used in the sum
t=0; %timestamp (do this for 0, 0.05, 0.1, 0.2, and 1)

%% Iteration 1, t=0
B = @(k) (4/((k^3)*(pi^3)))*(1-(-1)^k); %set B_k terms as function of k
lambda = @(k) pi * k / L; % set lambda_k terms similarly

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
title('Q1(b) - Fourier Series Solution Concentration Distribution')
hold on

%% Iteration 2, t=0.05
t_2=0.05;
B_2 = @(kb) (4/((kb^3)*(pi^3)))*(1-(-1)^kb); %set B_l terms as function of l
lambda_2 = @(kb) pi * kb / L; % set lambda_l terms similarly

x_2=(0:L/pts:1); % set up x values
u_2=zeros(1,length(x_2)); % zero vector of u values
for kb=1:M
    u_2 = u_2 + B_2(kb) * sin(lambda_2(kb)*x_2)*exp(-(kb^2)*(pi^2)*t_2); 
end
plot(x_2,u_2,'LineWidth',2.0)

%% Iteration 3, t=0.1
t_3=0.1;
B_3 = @(kc) (4/((kc^3)*(pi^3)))*(1-(-1)^kc); %set B_l terms as function of l
lambda_3 = @(kc) pi * kc / L; % set lambda_l terms similarly

x_3=(0:L/pts:1); % set up x values
u_3=zeros(1,length(x_3)); % zero vector of u values
for kc=1:M
    u_3 = u_3 + B_3(kc) * sin(lambda_3(kc)*x_3)*exp(-(kc^2)*(pi^2)*t_3); 
end
plot(x_3,u_3,'LineWidth',2.0)

%% Iteration 4, t=0.2
t_4=0.2;
B_4 = @(kd) (4/((kd^3)*(pi^3)))*(1-(-1)^kd); %set B_l terms as function of l
lambda_4 = @(kd) pi * kd / L; % set lambda_l terms similarly

x_4=(0:L/pts:1); % set up x values
u_4=zeros(1,length(x_4)); % zero vector of u values
for kd=1:M
    u_4 = u_4 + B_4(kd) * sin(lambda_4(kd)*x_4)*exp(-(kd^2)*(pi^2)*t_4); 
end
plot(x_4,u_4,'LineWidth',2.0)

%% Iteration 5, t=1
t_5=1;
B_5 = @(ke) (4/((ke^3)*(pi^3)))*(1-(-1)^ke); %set B_l terms as function of l
lambda_5 = @(ke) pi * ke / L; % set lambda_l terms similarly

x_5=(0:L/pts:1); % set up x values
u_5=zeros(1,length(x_5)); % zero vector of u values
for ke=1:M
    u_5 = u_5 + B_5(ke) * sin(lambda_5(ke)*x_5)*exp(-(ke^2)*(pi^2)*t_5); 
end
plot(x_5,u_5,'LineWidth',2.0)

%% Finalizing Graph
legend('t=0','t=0.05', 't=0.1', 't=0.2', 't=1')
hold off