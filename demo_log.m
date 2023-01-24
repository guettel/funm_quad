% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to compute log(I+A)/A*b for a 40x40 finite difference discretization of 
% the two-dimensional Laplace operator both with and without implicit
% deflation.
%%

f = @(x) log(1+x)/x;
N = 40;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
b = kron(e,e);
b = b/norm(b);

load exact_solutions

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
param.function = 'log';
param.restart_length = 20;              % Each restart cycle consists of 20 Arnoldi iterations
param.max_restarts = 20;                % We perform at most 20 restart cycles
param.tol = 1e-14;                      % Tolerance for quadrature rule
param.hermitian = 1;                    % set 0 if A is not Hermitian
param.V_full = 0;                       % set 1 if you need Krylov basis
param.H_full = 0;                       % do not store all Hessenberg matrices
param.exact = exact_poisson_log;        % Exact solution. If not known set to []
param.stopping_accuracy = 1e-16;        % stopping accuracy
param.inner_product = @(a,b) b'*a;      % Use standard euclidean inner product
param.thick = [];                       % No implicit deflation is performed
param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95 
param.waitbar = 1;                      % show waitbar 
param.reorth_number = 0;                % #reorthogonalizations
param.truncation_length = inf;          % truncation length for Arnoldi
param.verbose = 1;                      % print information about progress of algorithm

% compute log(I+A)/A*b using FUNM_QUAD without implicit deflation
tic
[f1,out1] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out1.err,'g--+')
hold on
for k = 2:length(out1.err)
    text(k+0.1,2*out1.err(k),num2str(out1.num_quadpoints(k)),'Color',[0 1 0],'FontSize',16,'Rotation',45);
end

%% adapt parameters for FUNM_QUAD algorithm (with with implicit deflation)
param.thick = @thick_quad;              % Thick restart function for implicit deflation
param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation

% compute log(I+A)/A*b using FUNM_QUAD with implicit deflation
tic
[f2,out2] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out2.err,'m--+')
for k = 2:length(out2.err)
    text(k+0.1,2*out2.err(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end

legend('without deflation','implicit deflation')
xlabel('cycle'); ylabel('absolute 2-norm error')
hold off
