% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to compute the exponential function of the 500x500 finite difference 
% discretization of a two-dimensional convection-diffusion problem for 
% convection parameters 0 (symmetric), 100 and 200 (non-symmetric) to 
% demonstrate the different behavior concerning speed of convergence 
% and required number of quadrature nodes for symmetric and non-symmetric 
% problems.


%% Build discretization matrix for 2D convection-diffusion problem 
N = 500;
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + 0*D1;

% choose time step s = 2*1e-3
s = 2*1e-3;
A = -s*A;

% choose right-hand side as normalized vector of all ones
b = ones(N^2,1); b = b/norm(b);

% load exact solutions
load exact_solutions

%% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 70;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-13;                  % tolerance for quadrature rule
param.hermitian = 1;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = exact_convdiff_1;     % exact solution. If not known set to []
param.stopping_accuracy = 1e-14;    % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95 
param.waitbar = 1;                  % show waitbar 
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi 
param.verbose = 1;                  % print information about progress of algorithm

%% compute exp(A)b by quadrature-based restart algorithm
tic
[f1,out1] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out1.err,'g--+')
for k = 2:length(out1.err)
    text(k+0.1,2*out1.err(k),num2str(out1.num_quadpoints(k)),'Color',[0 1 0],'FontSize',16,'Rotation',45);
end
hold all


%% Now increase the convection parameter to 100
A = D2 + 100*D1;
A = -s*A;
param.exact = exact_convdiff_2;     % update exact solution
param.hermitian = 0;                % the matrix A is no longer Hermitian!

% compute exp(A)b by quadrature based restart algorithm
tic
[f2,out2] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out2.err,'m-*')
for k = 2:length(out2.err)
    text(k+0.1,2*out2.err(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end

%% Further increase the convection parameter even further to 200
A = D2 + 200*D1;
A = -s*A;
param.exact = exact_convdiff_3;     % update exact solution

% compute exp(A)b by quadrature-based restart algorithm
tic
[f3,out3] = funm_quad(A,b,param);
toc

% plot convergence curve and number of quadrature points
semilogy(out3.err,'b-.x')
for k = 2:length(out3.err)
    text(k+0.1,2*out3.err(k),num2str(out3.num_quadpoints(k)),'Color',[0 0 1],'FontSize',16,'Rotation',45);
end

legend('\nu = 0','\nu = 100','\nu = 200')
xlabel('cycle'); ylabel('absolute 2-norm error')
axis([0,10.9,1e-15,10])

