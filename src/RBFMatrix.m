classdef RBFMatrix
    % RBFMatrix Row Butterfly Matrix
    % Row butterfly factorization A = F W_L ... W_1 V^* with constructor 
    % based on "Butterfly factorization via randomized matrix-vector
    % multiplications" by Liu et al. 2021.
    %
    % F and V are each block diagonal (implemented using BDMatrix).
    %
    % Each W factor is block diagonal, where each block is itself a 2x1
    % block matrix consisting of two block diagonal matrices. This nested 
    % class based structure is a little unwieldy, but is designed to 
    % make everything easy to visualize, and to allow most methods to 
    % operate at the matrix / factor level rather than the block level.
    %
    % Note that this is an O(n^2) implementation requiring O(n) matvecs
    % because for simplicity we do not use a hybrid factorization.
    properties
        tree 
        F
        W
        V
    end
    methods
        function A = RBFMatrix(fwd, adj, tree, r, s)
            if nargin == 0; A.tree = []; A.F = []; A.W = {}; A.V = []; return; end

            A.tree = tree;

            % number of points n
            i0 = tree.idx{1}; b0 = i0{1}; n = b0(2);

            % compute row bases using adjoint solves
            A.V = compute_row_bases(A, adj, r, s);

            %% compute factors W
            % initialize factors
            A.W = cell(tree.lvl,1);

            for l=1:tree.lvl
                % box indices in complementary downward pass
                idxc = tree.idx{l+1};

                % initialize block diagonal of factor
                W_blocks   = cell(2^(l-1),1);
                top_blocks = cell(2^(tree.lvl-l),1);
                bot_blocks = cell(2^(tree.lvl-l),1);

                for c=1:2^l
                    % sample sparse standard Gaussian random matrix
                    G = zeros(n,s);
                    G(range(idxc{c}),:) = randn(length(range(idxc{c})),s);
                    % compute adjoint solves of random matrix
                    Z = adj(G);
                    % apply row bases
                    Z = A.V'*Z;
                    for m=1:l-1
                        % apply factors up the tree
                        a = ancestor(tree, m, c, l);
                        Z = A.W{m}.D{fix((a-1)/2)+1}.M{rem(a-1,2)+1} * Z;
                    end
                    for b=1:2^(tree.lvl-l)
                        % extract row sample for block
                        Zb = Z(((b-1)*2*r+1):(b*2*r),:);

                        % compute transfer matrix from sample
                        if mod(c, 2) == 1
                            top_blocks{b} = col_basis(Zb, r)';
                        else
                            bot_blocks{b} = col_basis(Zb, r)';
                        end
                    end
                    if mod(c, 2) == 0
                        % write block
                        W_blocks{c/2} = BlockMatrix({BDMatrix(top_blocks); BDMatrix(bot_blocks)});
                    end
                end

                % write factor
                A.W{l} = BDMatrix(W_blocks);
            end

            % compute column bases using forward solves
            A.F = compute_column_bases(A, fwd, r, s);
        end
        function A = RBFMatrix_by_inversion(A, fwd, adj, tree, r, s, rtol)
            if nargin == 0; A.tree = []; A.F = []; A.W = {}; A.V = []; return; end

            A.tree = tree;

            % number of points n
            i0 = tree.idx{1}; b0 = i0{1}; n = b0(2);

            % compute row bases using adjoint solves
            A.V = compute_row_bases(A, adj, r, s);

            %% compute factors W
            % initialize factors
            A.W = cell(tree.lvl,1);

            for l=1:tree.lvl
                % box indices in complementary downward pass
                idxc = tree.idx{l+1};

                % initialize block diagonal of factor
                W_blocks   = cell(2^(l-1),1);
                top_blocks = cell(2^(tree.lvl-l),1);
                bot_blocks = cell(2^(tree.lvl-l),1);

                % sample sparse standard Gaussian random matrix
                G1 = zeros(n,s);
                G2 = zeros(n,s);
                for i=1:2^(l-1)
                    G1(range(idxc{2*i-1}),:) = randn(length(range(idxc{2*i-1})),s);
                    G2(range(idxc{2*i}),  :) = randn(length(range(idxc{2*i})),  s);
                end
                % compute adjoint solves of random matrix
                Z1 = adj(G1);
                Z2 = adj(G2);
                % apply row bases
                Z1 = A.V'*Z1;
                Z2 = A.V'*Z2;

                for c=1:2^(l-1)
                    for m=1:l-1
                        % apply inverse factors up the tree
                        Z1 = invert_factor(A, m, Z1, rtol);
                        Z2 = invert_factor(A, m, Z2, rtol);
                    end
                    for b=1:2^(tree.lvl-l)
                        % extract row sample for block
                        Z1b = Z1(((b-1)*2*r+1):(b*2*r),:);
                        Z2b = Z1(((b-1)*2*r+1):(b*2*r),:);

                        % compute transfer matrix from sample
                        top_blocks{b} = col_basis(Z1b, r)';
                        bot_blocks{b} = col_basis(Z2b, r)';
                    end
                    % write block
                    W_blocks{c} = BlockMatrix({BDMatrix(top_blocks); BDMatrix(bot_blocks)});
                end

                % write factor
                A.W{l} = BDMatrix(W_blocks);
            end

            % compute column bases using forward solves
            A.F = compute_column_bases(A, fwd, r, s);
        end
        function V = compute_row_bases(A, adj, r, s)
            %% compute row bases V

            % number of points n
            i0 = A.tree.idx{1}; b0 = i0{1}; n = b0(2);

            % sample standard Gaussian random matrix
            G = randn(n,s);
            % compute adjoint solves of random matrix
            Z = adj(G);

            % leaf indices
            idxl = A.tree.idx{A.tree.lvl+1};

            row_bases = cell(2^A.tree.lvl, 1);
            for b=1:2^A.tree.lvl
                % extract row sample for block
                Zb = Z(range(idxl{b}),:);

                % compute row basis
                row_bases{b} = col_basis(Zb, r);
            end

            % write row bases
            V = BDMatrix(row_bases);
        end
        function F = compute_column_bases(A, fwd, r, s)
            %% compute column bases F assuming V and W factors have already been computed

            % number of points n
            i0 = A.tree.idx{1}; b0 = i0{1}; n = b0(2);

            col_bases = cell(2^A.tree.lvl, 1);
            % sample standard Gaussian random matrix
            G = randn(n,s);
            % compute forward solves of random matrix
            Z = fwd(G);
            % apply row bases (this time to G because it's a forward solve)
            G = A.V'*G;

            % leaf indices
            idxl = A.tree.idx{A.tree.lvl+1};

            for m=1:A.tree.lvl
                % apply factors up the tree
                G = A.W{m}*G;
            end
            for b=1:2^A.tree.lvl
                % extract row sample for block
                Zb = Z(range(idxl{b}),:);
                Gb = G(((b-1)*r+1):(b*r),:);

                % compute column basis from sample by solving least squares problem 
                col_bases{b} = Zb * pinv(Gb);
            end

            % write column bases
            F = BDMatrix(col_bases);
        end
        function s = size(A, dim)
            iendd = A.tree.idx{1}; iend = iendd{1};
            if nargin > 1
                if dim <= 2
                    s = iend(2);
                else
                    error("dimension must be either 1 or 2")
                end
            else
                s = [iend(2), iend(2)];
            end
        end
        function At = ctranspose(A)
            At = RBFMatrix();
            At.tree = A.tree;
            At.V = A.F;
            At.F = A.V;
            for l=1:A.tree.lvl
                At.W{A.tree.lvl-l+1} = A.W{l}';
            end
        end
        function B = mtimes(A, X)
            B = A.V' * X;
            for l=1:A.tree.lvl
                B = A.W{l} * B;
            end
            B = A.F * B;
        end
        function M = dense(A)
            M = A * eye(size(A));
        end
        function X = invert_factor(A, l, B, rtol)
            X = zeros(size(B));
            r = size(A.W{l}.D{1}.M{1}.D{1},1);
            i0 = 1;
            for b=1:length(A.W{l}.D)
                for c=1:length(A.W{l}.D{b}.M{1}.D)
                    [U,S,V] = svd([A.W{l}.D{b}.M{1}.D{c}; A.W{l}.D{b}.M{2}.D{c}]);
                    % rtol-rank of block
                    rb = sum(diag(S) / max(diag(S)) >= rtol);

                    i1 = 2^(A.tree.lvl-l+1)*r*(b-1) + r*(c-1) + 1;
                    i2 = i1 + 2^(A.tree.lvl-l)*r;
                    perm_inds = [i1:i1+r-1, i2:i2+r-1];

                    X(i0:i0+2*r-1,:) = V(:,1:rb)*(S(1:rb,1:rb)\(U(:,1:rb)'*B(perm_inds,:)));
                    i0 = i0 + 2*r;
                end
            end
        end
        function BF = to_yang_format(A, r)
            for b=1:2^A.tree.lvl
                BF.level(1).blocks(1,b).mat            = A.V.D{b}';
                BF.level(A.tree.lvl+2).blocks(b,1).mat = A.F.D{b};
            end
            for l=1:A.tree.lvl
                for i=1:2^l
                    for j=1:2^(A.tree.lvl - l)
                        block = A.W{l}.D{fix((i-1)/2)+1}.M{rem(i-1,2)+1}.D{j};
                        BF.level(l+1).blocks(i,2*j-1).mat = block(:,1:r);
                        BF.level(l+1).blocks(i,2*j).mat   = block(:,r+1:end);
%                         D{fix((a-1)/2)+1}.M{rem(a-1,2)+1}
                    end
                end
            end
        end
    end
end