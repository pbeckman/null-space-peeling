%% Compute HSS factorization from matvecs

% size of matrix
n   = 1000;
% rank of block row
r   = 10;
% number of random matvecs
s   = 3*r;
% level of factorization hierarchy
lvl = floor(log2(n)-log2(r));

fprintf("\nn: %i, r: %i, s: %i, level: %i\n", n, r, s, lvl)

% fast periodic convolution with exp(-2*pi*c*abs(x))
f   = [0:floor(n/2) -ceil(n/2)+1:-1];
x   = linspace(0, 1, n);
c   = 1;
G   = n / pi * c ./ (f.^2 + c^2);
fwd = @(v) ifft(G' .* fft(v, [], 1), [], 1);
adj = fwd;

% form index tree
tree = IndexTree(n, lvl);

% compute factorization
tic;
A = HSSMatrix(fwd, adj, tree, r, s);
t = toc;
fprintf("Factorization time: %.2e s\n", t);

%% Compute error

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
imagesc(K);
axis square
colorbar
title('$K$','Interpreter','latex','FontSize',24)

subplot('Position', [0.35, 0.1, 0.3, 0.8]);
imagesc(HK);
axis square
colorbar
title(strcat('$\tilde{K}$', sprintf(' (rank %i, level %i)', r, lvl)),'Interpreter','latex','FontSize',24)

ax = subplot('Position', [0.68, 0.1, 0.3, 0.8]);
imagesc(log10(abs(K - HK)));
axis square
colormap(ax, turbo)
colorbar
title('$\log_{10}|K - \tilde{K}|$','Interpreter','latex','FontSize',24)

