% This demonstrates the behavior of quadrature-based restart algorithm 
% FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% for different restart lengths m = 10,...,50. The test is to compute 
% the inverse square root of a 100x100 finite difference discretization 
% of the Laplace operator.

%%
% Initialize 2D Laplacian
N = 100;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
s = eigs(A,1,'SM');
A = A/s;
b = kron(e,e);
b = b/norm(b);

load exact_solutions

addpath('funm_quad')
%%
% Run FUNM_QUAD algorithm 5 times for different restart lengths

param.function = 'invSqrt';
param.tol = 1e-14;                      % tolerance for quadrature rule
param.transformation_parameter = 1;     % parameter for the integral transformation
param.hermitian = 1;                    % set 0 if A is not Hermitian
param.V_full = 0;                       % set 1 if you need Krylov basis
param.H_full = 0;                       % do not store all Hessenberg matrices
param.exact = exact_poisson_invsqrt;    % exact solution. If not known set to []
param.stopping_accuracy = 1e-16;        % stopping accuracy
param.inner_product = @(a,b) b'*a;      % use standard euclidean inner product
param.thick = [];                       % no implicit deflation is performed
param.min_decay = .95;                  % we desire linear error reduction of rate < .95
param.waitbar = 1;                      % show waitbar
param.reorth_number = 0;                % #reorthogonalizations
param.truncation_length = inf;          % truncation length for Arnoldi
param.verbose = 1;                      % print information about progress of algorithm


for m = 10:10:50,
    
    % adapt restart length and max_restarts
    param.restart_length = m;               % each restart cycle consists of m Arnoldi iterations
    param.max_restarts = 1000/m;            % choose number of restarts so that overall at most 1000 matrix vector products are performed
    
    % compute A^(-1/2)*b by quadrature based restart algorithm
    tic
    [f1,out1] = funm_quad(A,b,param);
    toc
    x = 1:length(out1.err);
    x = m*x;
    % plot convergence curves
    if m == 10
        semilogy(x,out1.err,'k-+')
    end
    if m == 20
        semilogy(x,out1.err,'b-o')
    end
    if m == 30
        semilogy(x,out1.err,'g-s')
    end
    if m == 40
        semilogy(x,out1.err,'r-*')
    end
    if m == 50
       semilogy(x,out1.err,'m-d') 
    end
    hold all
end

legend('m = 10','m = 20','m = 30','m = 40','m = 50')
xlabel('number of matrix-vector products')
ylabel('absolute 2-norm error')

