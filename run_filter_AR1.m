close all
clear all
clc

clear all;
randn('state',1)

%% generate data
periods=1000;
rho=0.95;
sigma_bar=log(0.01);
rho_sigma=0.9;
eta_sigma=0.3;

siggma=zeros(periods,1);
x=zeros(periods,1);
siggma(1)=sigma_bar;
epsil=randn(periods,2);
for jj=2:periods
    siggma(jj)=(1-rho_sigma)*sigma_bar+rho_sigma*siggma(jj-1)+eta_sigma*epsil(jj,1);
    x(jj)=rho*x(jj-1)+exp(siggma(jj))*epsil(jj,2);
end 
true_par=[rho_sigma,rho,eta_sigma,sigma_bar];
save data_simulated x siggma rho sigma_bar rho_sigma eta_sigma true_par


%% settings

num_sim_filter          = 200;         % Number of particles for a particle filter                            
num_sim_smoother        = 200;         % Number of particles for a particle filter                            
dimy= 1;          % Number of variables in the set of observables (y)
burnin = 500;
MH_draws = 2; 

% *************************LOAD DATA **************************************
load data_simulated x rho sigma_bar rho_sigma eta_sigma

observable_series=x;
observable_series=observable_series(~isnan(observable_series),:);

% transform to unbounded
[rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans]=par_transform_AR1(rho_sigma,rho,eta_sigma,sigma_bar);
npar=4;

%% initialize random numbers
PFstreamrandn = RandStream('mt19937ar','Seed',2);
% The normally distributed shocks to the state equation for the particle filter
shocks      = randn(PFstreamrandn,num_sim_filter,periods);

% The random numbers (from a uniform(0,1) distribution) needed for the
% resampling routine in the particle filter
PFstreamrand = RandStream('mt19937ar','Seed',2);
randnr             = rand(PFstreamrand,periods,num_sim_filter);


%initialize matrices
draws=zeros(MH_draws+burnin,npar);
likelihood=zeros(MH_draws+burnin,1);


draws(1,:)=[rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans];

%create save name
tempdate=deblank(num2str(fix(clock)));
tempdate=[tempdate(1:4),'_',strtrim(tempdate(7:10)),'_',strtrim(tempdate(13:17)),'_',strtrim(tempdate(20:24)),'_',strtrim(tempdate(27:30))];
fname=['simulation_',tempdate];


%% MCMC

scale_mh = 2;
accept=0;

inv_Hessian = 1e-3*eye(4);
Sigma_chol = cholcov(inv_Hessian)';
old_posterior = PF_caller(draws(1,:),observable_series,num_sim_filter,num_sim_smoother,shocks,randnr,0);
likelihood(1,1)=old_posterior;

proposal_draws = scale_mh*Sigma_chol*randn(4,MH_draws+burnin);
for ii=2:burnin+MH_draws
    xhatstar = draws(ii-1,:)+proposal_draws(:,ii)';
    new_posterior = PF_caller(xhatstar,observable_series,num_sim_filter,num_sim_smoother,shocks,randnr,0);
    likelihood(ii,1)=new_posterior;
    accprob=exp(new_posterior-old_posterior);
    if rand(1,1)<=accprob
        draws(ii,:)=xhatstar;
        old_posterior = new_posterior;
        accept=accept+1;
    else
        draws(ii,:)=draws(ii-1,:);
    end
    if mod(ii,100)==0
        ratio=accept/ii;
        disp(['Acceptance Rate:' num2str(ratio)])
        eval(['save ',fname]);
    end
end

[rho_sigma_vec,rho_1_vec,eta_sigma_vec,sigma_bar_vec]=par_retransform_AR1(draws(1:ii-1,:));
figure
subplot(4,1,1)
plot(rho_sigma_vec)
title('\rho_\sigma')
axis tight
subplot(4,1,2)
plot(rho_1_vec)
title('\rho')
axis tight
subplot(4,1,3)
plot(eta_sigma_vec)
title('\eta_\sigma')
axis tight
subplot(4,1,4)
plot(sigma_bar_vec)
title('\bar \sigma')
axis tight

eval(['save ',fname]);

%% Run Smoother
num_sim_smoother = 200;

% make sure that parameters are transformed
PF_caller([rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans],observable_series,num_sim_filter,num_sim_smoother,shocks,randnr,1,'Stoch_vol_AR_smoother');

load Stoch_vol_AR_smoother vola_resids level_resids x_smoother_store xhat_smoother x_resampled_store

figure
subplot(3,1,1)
plot(1:periods,xhat_smoother,'r*','LineWidth',0.5)
title('Volatility')
hold on
plot(1:periods,siggma,'b-','LineWidth',0.5)
plot(1:periods,squeeze(median(x_resampled_store,2)),'c-.','LineWidth',0.5)
legend('Smoothed','True','Filtered','Location','NorthWest')
fprintf('Correlation States: %4.3f\n',corr(xhat_smoother',siggma))

subplot(3,1,2)
plot(1:periods,median(vola_resids,2)/eta_sigma,'r*','LineWidth',0.5)
hold on
plot(1:periods,eta_sigma*epsil(:,1),'b-','LineWidth',0.5)
title('Vola Residuals')
fprintf('Correlation State Resids: %4.3f\n',corr(median(vola_resids,2),eta_sigma*epsil(:,1)))
fprintf('Std State Resids: %4.3f\n',std(vola_resids(:)/eta_sigma))

subplot(3,1,3)
plot(1:periods,median(level_resids,2),'r*','LineWidth',0.5)
hold on
plot(1:periods,epsil(:,2),'b-','LineWidth',0.5)
title('Level Residuals')
fprintf('Correlation Level Resids: %4.3f\n',corr(median(level_resids,2),epsil(:,2)))
fprintf('Std Level Resids: %4.3f\n',std(level_resids(:)))

