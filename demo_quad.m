% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to test different restart procedures for f(z) = (1-exp(-s*sqrt(z)))/z.
%
% It demonstrates how to compute the action of a general matrix 
% function of the form \int_{\inf}^0 g(t)/(t-z) dt with funm_quad.
%
% To this end, the integrand has to be given as a function handle of the 
% form param.function = @(z,t) g(t)./(t-z), which is then passed on to
% (a minor modification of) Matlab's INTEGRAL function.
%
% Using this functionality is only recommended for Hermitian matrices A, 
% as a diagonalization of the reduced Hessenberg matrix is performed in the 
% current implementation, which can cause instabilities for non-normal 
% matrices. This is not an inherent limitation of the restart algorithm; 
% it is only necessary because INTEGRAL does not support matrix function 
% integrands.

%% Define f, A, and b.
s = 1e-3;
f = @(x) (1-exp(-s*sqrt(x)))/x;
N = 100;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N); 
b = kron(e,e);
b = b/norm(b);

load exact_solutions

%% Choose parameters for FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
param.restart_length = 50;                                      % each restart cycle consists of 50 Arnoldi iterations
param.max_restarts = 14;                                        % perform at most 20 restart cycles
param.function = @(z,t) 1./(t-z) .* sin(s*sqrt(-t))./(pi*t);    % integrand of integral representation of function f
param.tol = 1e-14;                                              % tolerance for quadrature rule
param.hermitian = 1;                                            % set 0 if A is not Hermitian
param.V_full = 0;                                               % set 1 if you need Krylov basis
param.H_full = 0;                                               % do not store all Hessenberg matrices
param.exact = exact_poisson_quad;                               % Exact solution. If not known set to []
param.stopping_accuracy = 1e-16;                                % stopping accuracy
param.inner_product = @(a,b) b'*a;                              % use standard euclidean inner product
param.thick = [];                                               % no implicit deflation is performed
param.min_decay = 0.95;                                         % we desire linear error reduction of rate < .95 
param.waitbar = 1;                                              % show waitbar 
param.reorth_number = 0;                                        % #reorthogonalizations
param.truncation_length = inf;                                  % truncation length for Arnoldi
param.verbose = 2;                                              % print more detailed information about progress of algorithm

%% Compute f(A)*b by FUNM_QUAD without implicit deflation
tic
[f1,out1] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out1.err,'g--s','Color',[0,.85,0])
hold on
for k = 2:length(out1.err)
    text(k+0.1,2*out1.err(k),num2str(out1.num_quadpoints(k)),'Color',[0 .85 0],'FontSize',16,'Rotation',45);
end

%% Choose parameters for FUNM_QUAD algorithm (this time with implicit deflation)
param.thick = @thick_quad;                                      % thick-restart function
param.number_thick = 5;                                         % Number of target eigenvalues for implicit deflation    


%% Compute f(A)*b by FUNM_QUAD with implicit deflation
tic
[f2,out2] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out2.err,'m--+')
legend('without deflation','with deflation')
xlabel('cycle'); ylabel('absolute 2-norm error')
axis([0,15,1e-16,1e-4])
for k = 1:length(out2.err)
    text(k+0.1,2*out2.err(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end
hold off

