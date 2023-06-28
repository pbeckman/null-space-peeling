%% Compute row butterfly factorization from matvecs

% size of matrix
n   = 45;
% rank of block row
r   = 10;
% number of random matvecs
s   = r + 10;
% level of factorization hierarchy
lvl = floor(log2(n)-log2(r));
% construct index tree
tree = IndexTree(n, lvl);

fprintf("\nn: %i, r: %i, s: %i, level: %i\n", n, r, s, lvl)

% choose an example from below
example = "random";

if example == "low-rank"
    % symmetric low rank kernel M*M'
    r_true = 10;
    M = flip(sort(rand(n, r_true)));
    fwd = @(v) M*(M'*v);
    adj = fwd;
elseif example == "gaussian"
    % fast periodic convolution with Gaussian kernel
    % numerically low-rank depending on c
    f   = [0:floor(n/2) -ceil(n/2)+1:-1];
    x   = linspace(0, 1, n);
    c   = 2;
    G   = n * exp(-f.^2 / c.^2);
    fwd = @(v) ifft(G' .* fft(v, [], 1), [], 1);
    adj = fwd;
elseif example == "DFT"
    % apply the DFT matrix using the FFT
    fwd = @(v) 1/sqrt(n) *  fft(v, [], 1);
    adj = @(v)   sqrt(n) * ifft(v, [], 1);
elseif example == "random"
    % random butterfly matrix
    r_true = 10;
    B = RBFMatrix();
    B.tree = tree;
    B.V = BDMatrix(cellfun(@(idx) randn(length(range(idx)),r_true), B.tree.idx{lvl+1}, 'UniformOutput', false));
    B.F = BDMatrix(cellfun(@(idx) randn(length(range(idx)),r_true), B.tree.idx{lvl+1}, 'UniformOutput', false));
    for l=1:lvl
        W_blocks = cell(2^(l-1),1);
        for b=1:2^(l-1)
            W_blocks{b} = BlockMatrix({ ...
                BDMatrix(arrayfun(@(~) randn(r_true, 2*r_true), 1:2^(lvl-l), 'UniformOutput', false)); ...
                BDMatrix(arrayfun(@(~) randn(r_true, 2*r_true), 1:2^(lvl-l), 'UniformOutput', false))  ...
                });
        end
        B.W{l} = BDMatrix(W_blocks);
    end
    fwd = @(v) B*v;
    adj = @(v) B'*v;
else
    error("Example '%s' not found - please run a valid example.", example)
end



% % reordered DFT matrix
% B = RBFMatrix();
% B.tree = tree;
% B.V = BDMatrix(arrayfun(@(~) 1, 1:n, 'UniformOutput', false));
% B.F = BDMatrix(arrayfun(@(~) 1, 1:n, 'UniformOutput', false));
% for l=1:lvl
%     W_blocks = cell(2^(l-1),1);
%     for b=1:2^(l-1)
%         W_blocks{b} = BlockMatrix({ ...
%             BDMatrix(arrayfun(@(j) [1  exp(-1i*j/2^(lvl-l+1))], 0:2^(lvl-l)-1, 'UniformOutput', false)); ...
%             BDMatrix(arrayfun(@(j) [1 -exp(-1i*j/2^(lvl-l+1))], 0:2^(lvl-l)-1, 'UniformOutput', false)); ...
%             });
%     end
%     B.W{l} = BDMatrix(W_blocks);
% end
% fwd = @(v) B*v;
% adj = @(v) B'*v;

%% Compute factorization

tic;
verbose = false;
A = RBFMatrix(fwd, adj, tree, r, s, verbose);
t = toc;
fprintf("\nFactorization time: %.2e s\n", t);

%% Compute error

% number of matvec samples to compute error
ns = 100;
M = randn(n, ns);
fprintf("\nRelative matvec error: %.3e\n", norm(fwd(M) - A*M, 'fro') / norm(M) / ns)

% compute dense version of operator and compressed matrix for error computations
tic; K = fwd(eye(n)); t = toc;
fprintf("Operator matvec time:  %.2e s\n", t/n);
tic; BK = A*eye(n); t = toc;
fprintf("Butterfly matvec time: %.2e s\n", t/n);

fprintf("\nRelative Frobenius error: %.3e\n", norm(K - BK, 'fro') / norm(K, 'fro'))
fprintf("Dense matrix:     %.2f MB\n", whos('K').bytes * 9.53674e-7)
fprintf("Butterfly matrix: %.2f MB\n", whos('A').bytes * 9.53674e-7)


%% Plot dense matrix, compressed matrix, and log relative error

figure(1)
clf

colormap(parula)

subplot('Position', [0.02, 0.1, 0.3, 0.8]);
imagesc(real(K));
axis square
colorbar
title('$K$','Interpreter','latex','FontSize',24)

subplot('Position', [0.35, 0.1, 0.3, 0.8]);
imagesc(real(BK));
axis square
clim([min(real(K), [], "all"), max(real(K), [], "all")])
colorbar
title(strcat('$\tilde{K}$', sprintf(' (rank %i, level %i)', r, lvl)),'Interpreter','latex','FontSize',24)

ax = subplot('Position', [0.68, 0.1, 0.3, 0.8]);
imagesc(log10(abs(K - BK) ./ abs(K)));
axis square
colormap(ax, turbo)
colorbar
title('$\log_{10}|K - \tilde{K}| / |K|$','Interpreter','latex','FontSize',24)