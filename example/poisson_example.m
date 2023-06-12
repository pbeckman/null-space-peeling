%% Compute HSS factorization from matvecs

% number of vertical points in grid (size of schur complement matrix)
n   = 1000;
% number of horizontal points in grid
m   = 51; if mod(m,2) == 0; m = m+1; end
% rank of block row
r   = 15;
% number of random matvecs
s   = 3*r;
% level of factorization hierarchy
lvl = floor(log2(n)-log2(r));

fprintf("\nn: %i, r: %i, s: %i, level: %i\n", n, r, s, lvl)

% 5-point finite difference negative Laplacian on a rectangle 
% with homogeneous Dirichlet boundary conditions
ey = ones(n,1);
Ty = spdiags([-ey 2*ey -ey], -1:1, n, n);

ex = ones(fix(m/2),1);
Tx = spdiags([-ex 2*ex -ex], -1:1, fix(m/2), fix(m/2));

Cii = kron(speye(fix(m/2)),Ty) + kron(Tx,speye(n));

C13 = sparse(n*fix(m/2), n);
C13(end-n+1:end,:) = -speye(n);

C23 = sparse(n*fix(m/2), n);
C23(1:n,:) = -speye(n);

C33 = spdiags([-ey 4*ey -ey], -1:1, n, n);

% compute Schur complement using sparse matvecs and linear solves
p = dissect(Cii);
R = chol(Cii(p,p));
fwd = @(v) C33*v - C13(p,:)'*(R\(R'\(C13(p,:)*v))) - C23(p,:)'*(R\(R'\(C23(p,:)*v)));
adj = fwd;

% form index tree
tree = IndexTree(n, lvl);

% compute factorization
tic;
A = HSSMatrix(fwd, adj, tree, r, s);
t = toc;
fprintf("Factorization time: %.2e s\n", t);

%% Compute error (this is slow for large n because of the matvec speed)

% compute dense version of operator and compressed matrix for error computations
tic; K = fwd(eye(n)); t = toc;
fprintf("Operator matvec time: %.2e s\n", t/n);
tic; HK = A*eye(n); t = toc;
fprintf("HSS matvec time:      %.2e s\n", t/n);

fprintf("Dense matrix: %.2f MB\n", whos('K').bytes * 9.53674e-7)
fprintf("HSS   matrix: %.2f MB\n", whos('A').bytes * 9.53674e-7)
fprintf("Frobenius error: %.3e\n", norm(K - HK,'fro'))

%% Plot dense and compressed matrices and the log error

figure(1)
clf

colormap(parula) 

subplot('Position', [0.02, 0.1, 0.3, 0.8]);
imagesc(log10(abs(K)));
axis square
colorbar
clim([-16,0])
title('$\log_{10}|K|$','Interpreter','latex','FontSize',24)

subplot('Position', [0.35, 0.1, 0.3, 0.8]);
imagesc(log10(abs(HK)));
axis square
colorbar
clim([-16,0])
title(strcat('$\log_{10}|\tilde{K}|$', sprintf(' (rank %i, level %i)', r, lvl)),'Interpreter','latex','FontSize',24)

ax = subplot('Position', [0.68, 0.1, 0.3, 0.8]);
imagesc(log10(abs(K - HK)));
axis square
colormap(ax, turbo)
colorbar
title('$\log_{10}|K - \tilde{K}|$','Interpreter','latex','FontSize',24)