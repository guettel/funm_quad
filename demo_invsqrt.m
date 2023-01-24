% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to compute the inverse square root of the 100x100 finite difference 
% discretization of the two-dimensional Laplace operator both with and 
% without implicit deflation.

%%
% Initialize 2D Laplacian
N = 100;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
s = eigs(A,1,'SM');
A = A/s;
b = kron(e,e);
b = b/norm(b);

% Load exact solution
load exact_solutions

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
param.function = 'invSqrt';
param.restart_length = 50;              % each restart cycle consists of 50 Arnoldi iterations
param.max_restarts = 20;                % perform at most 20 restart cycles
param.tol = 1e-14;                      % tolerance for quadrature rule
param.transformation_parameter = 1;     % parameter for the integral transformation
param.hermitian = 1;                    % set 0 if A is not Hermitian
param.V_full = 0;                       % set 1 if you need Krylov basis
param.H_full = 0;                       % do not store all Hessenberg matrices
param.exact = exact_poisson_invsqrt;    % Exact solution. If not known set to []
param.stopping_accuracy = 1e-16;        % stopping accuracy
param.inner_product = @(a,b) b'*a;      % use standard euclidean inner product
param.thick = [];                       % no implicit deflation is performed
param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95 
param.waitbar = 1;                      % show waitbar 
param.reorth_number = 0;                % #reorthogonalizations
param.truncation_length = inf;          % truncation length for Arnoldi
param.verbose = 1;                      % print information about progress of algorithm

% Compute A^(-1/2)*b using FUNM_QUAD without implicit deflation
tic
[f1,out1] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out1.err,'g-.s','Color',[0,.85,0])
hold on
for k = 2:length(out1.err)
    text(k+0.1,3*out1.err(k),num2str(out1.num_quadpoints(k)),'Color',[0 .85 0],'FontSize',16,'Rotation',45);
end


%% Now with implicit deflation, adapt the FUNM_QUAD parameters
param.thick = @thick_quad;              % Thick restart function for implicit deflation
param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation

% compute A^(-1/2)*b using FUNM_QUAD with implicit deflation
tic
[f2,out2] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out2.err,'m--+')
for k = 2:length(out2.err)
    text(k+0.1,2*out2.err(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end

axis([0,20,1e-14,1])
legend('without deflation','implicit deflation')
xlabel('cycle'); ylabel('absolute 2-norm error')
hold off

