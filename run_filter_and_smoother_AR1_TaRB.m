% Generate Data from AR1-stochastic volatility model, run MCMC using particle filter on it. Then use smoother.
% 
% Notes: 
% - This code does not do mode-finding and computing Hessian at the mode, because the likelihood from the particle filter is not differentiable. 
%   Instead, a diagonal matrix is used. This implies loss of efficiency but any positive definite matrix will still yield draws from the ergodic distribution 
%   (see Chib/Greenberg (1995))
% - To perform mode-finding using the CMA-ES algorithm, uncomment the code in lines 101 following and 
%   download the cmaes.m from the provided link.
% - Using an insufficient number of particles will yield problems with the acceptance rate, see Pitt et al. (2012):
%   "On some properties of Markov chain Monte Carlo simulation methods based on the particle filter", Journal of Econometrics, 171, 134-151
% - In contrast to SMM approaches, the simulations should not be conducted with fixed random numbers, unless the number of particles is really large. See 
%   Thomas Flury and Neil Shephard (2011): "Bayesian Inference based only on simulated likelihood", Econometric Theory, 27, 933–956
% 
% Copyright (C) 2013-2015 Benjamin Born + Johannes Pfeifer
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% See <http://www.gnu.org/licenses/>.

close all
clear all
clc

randn('state',1)
addpath('Tools')
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

num_sim_filter          = 1000;         % Number of particles for a particle filter                            
num_sim_smoother        = 1000;         % Number of particles for a particle filter                            
dimy= 1;          % Number of variables in the set of observables (y)
burnin = 500;
MH_draws = 5000; 

% *************************LOAD DATA **************************************
load data_simulated x rho sigma_bar rho_sigma eta_sigma

observable_series=x;
observable_series=observable_series(~isnan(observable_series),:);

% transform to unbounded support
[rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans]=par_transform_AR1(rho_sigma,rho,eta_sigma,sigma_bar);
npar=4;

%% initialize random numbers

%initialize matrices
draws=zeros(MH_draws+burnin,npar);
likelihood=zeros(MH_draws+burnin,1);

%create save name depending on clock time
tempdate=deblank(num2str(fix(clock)));
tempdate=[tempdate(1:4),'_',strtrim(tempdate(7:10)),'_',strtrim(tempdate(13:17)),'_',strtrim(tempdate(20:24)),'_',strtrim(tempdate(27:30))];
fname=['simulation_',tempdate];


%% MCMC

scale_mh = 1e-8; %scale factor for proposal density to get acceptance rate of 23-40 percent
accept=0;
new_block_probability=0.3;

crit=1e-4;
nit=100;

%set startig value, starts at true values and makes mode-finding redundant
draws(1,:)=[rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans];

% set jumping matrix for proposal density
% uses simple identity matrix instead of Hessian at the mode
inv_Hessian = 1e-3*eye(4); 
Sigma_chol = cholcov(inv_Hessian)';

last_posterior = PF_caller(draws(1,:),observable_series,num_sim_filter,num_sim_smoother,0);
likelihood(1,1)=last_posterior;

proposal_draws = scale_mh*Sigma_chol*randn(4,MH_draws+burnin);
blocked_draws_counter=0;
accepted_draws_counter=0;
for draw_iter=2:burnin+MH_draws
    %% randomize indices for blocking in this iteration
    indices=randperm(npar)';
    blocks=[1; (1+cumsum((rand(length(indices)-1,1)>(1-new_block_probability))))];
    nblocks=blocks(end,1); %get number of blocks this iteration
    current_draw=draws(draw_iter-1,:)'; %get starting point for current draw for updating
    for block_iter=1:nblocks
        blocked_draws_counter=blocked_draws_counter+1;
        nxopt=length(indices(blocks==block_iter,1)); %get size of current block
        par_start_current_block=current_draw(indices(blocks==block_iter,1));
        [fval,xopt_current_block,~,hessian_mat,~,~,exitflag] = ...
            csminwel('PF_caller_TaRB',par_start_current_block,1e-5*eye(nxopt),[],crit,nit,...
                current_draw,indices(blocks==block_iter,1),observable_series,num_sim_filter,num_sim_smoother,0); %inputs for objective

        if any(any(isnan(hessian_mat))) || any(any(isinf(hessian_mat)))
            inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
        else
            inverse_hessian_mat=inv(hessian_mat); %get inverse Hessian
            if any(any((isnan(inverse_hessian_mat)))) || any(any((isinf(inverse_hessian_mat))))
                inverse_hessian_mat=eye(nxopt)*1e-4; %use diagonal
            end
        end
        [proposal_covariance_Cholesky_decomposition_upper,negeigenvalues]=chol(inverse_hessian_mat);
        %if not positive definite, use generalized Cholesky of Eskow/Schnabel
        if negeigenvalues~=0
            proposal_covariance_Cholesky_decomposition_upper=chol_SE(inverse_hessian_mat,0);
        end
        proposal_covariance_Cholesky_decomposition_upper=proposal_covariance_Cholesky_decomposition_upper*scale_mh;
        %get proposal draw

        tempdraw=current_draw;
        r = randn(1,nxopt) * proposal_covariance_Cholesky_decomposition_upper;
        proposed_par=xopt_current_block+r';
        tempdraw(indices(blocks==block_iter,1))=proposed_par;

        % check whether draw is valid and compute posterior
        try
            logpost = - (PF_caller_TaRB(tempdraw(indices(blocks==block_iter,1)),current_draw,indices(blocks==block_iter,1),observable_series,num_sim_filter,num_sim_smoother,0)); %inputs for objective
        catch
            logpost = -inf;
        end
        
        if (logpost > -inf)
            %get ratio of proposal densities, required because proposal depends
            %on current mode via Hessian and is thus not symmetric anymore
            proposal_density_proposed_move_forward=multivariate_normal_pdf(proposed_par',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,nxopt);
            proposal_density_proposed_move_backward=multivariate_normal_pdf(par_start_current_block',xopt_current_block',proposal_covariance_Cholesky_decomposition_upper,nxopt);
            accprob=logpost-last_posterior+ log(proposal_density_proposed_move_backward)-log(proposal_density_proposed_move_forward); %Formula (6), Chib/Ramamurthy
            exp(accprob)
            if  (log(rand) < accprob)
                current_draw(indices(blocks==block_iter,1))=proposed_par;
                last_posterior=logpost;
                accepted_draws_counter =accepted_draws_counter +1;
            else %no updating
                %do nothing, keep old value
            end
        end
    end
    accepted=accepted_draws_counter/blocked_draws_counter;
    draws(draw_iter,:) = current_draw;
    logpost = last_posterior; %make sure not a temporary draw is returned;
    if mod(draw_iter,10)==0
        disp(['Acceptance Rate:' num2str(accepted)])
        eval(['save ',fname]);
    end
end
% retransform parameters and plot them
[rho_sigma_vec,rho_1_vec,eta_sigma_vec,sigma_bar_vec]=par_retransform_AR1(draws(burnin:draw_iter,:));
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

% comment the following four line out to obtain smoothed estimates at the
% true values

rho_sigma_trans=mean(draws(burnin:draw_iter,1));
rho_trans=mean(draws(burnin:draw_iter,2));
eta_sigma_trans=mean(draws(burnin:draw_iter,3));
sigma_bar_trans=mean(draws(burnin:draw_iter,4));

% make sure to use parameters as they are the input arguments!
PF_caller([rho_sigma_trans,rho_trans,eta_sigma_trans,sigma_bar_trans],observable_series,num_sim_filter,num_sim_smoother,1,'Stoch_vol_AR_smoother');

% load smoother results
load Stoch_vol_AR_smoother vola_resids level_resids x_smoother_store xhat_smoother x_resampled_store

%plot them
figure
subplot(3,1,1)
plot(1:periods,xhat_smoother,'r--','LineWidth',0.5)
title('Volatility')
hold on
plot(1:periods,siggma,'b-','LineWidth',0.5)
plot(1:periods,squeeze(median(x_resampled_store,2)),'c-.','LineWidth',0.5)
legend('Smoothed','True','Filtered','Location','NorthWest')
fprintf('Correlation States: %4.3f\n',corr(xhat_smoother',siggma))

subplot(3,1,2)
plot(1:periods,median(vola_resids,2)/eta_sigma,'r--','LineWidth',0.5)
hold on
plot(1:periods,eta_sigma*epsil(:,1),'b-','LineWidth',0.5)
title('Vola Residuals')
fprintf('Correlation State Resids: %4.3f\n',corr(median(vola_resids,2),eta_sigma*epsil(:,1)))
fprintf('Std State Resids: %4.3f\n',std(vola_resids(:)/eta_sigma))

subplot(3,1,3)
plot(1:periods,median(level_resids,2),'r--','LineWidth',0.5)
hold on
plot(1:periods,epsil(:,2),'b-','LineWidth',0.5)
title('Level Residuals')
fprintf('Correlation Level Resids: %4.3f\n',corr(median(level_resids,2),epsil(:,2)))
fprintf('Std Level Resids: %4.3f\n',std(level_resids(:)))


