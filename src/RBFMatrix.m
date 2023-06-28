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
    properties
        tree 
        F
        W
        V
    end
    methods
        function A = RBFMatrix(fwd, adj, tree, r, s, verbose)
            if nargin == 0; A.tree = []; A.F = []; A.W = {}; A.V = []; return; end
            if nargin < 6; verbose = false; end

            A.tree = tree;

            % number of points n
            i0 = tree.idx{1}; b0 = i0{1}; n = b0(2);

            %% compute row bases V
            if verbose; fprintf("\n--- Level 0 ---\n"); end
            % sample standard Gaussian random matrix
            G = randn(n,s);
            % compute adjoint solves of random matrix
            Z = adj(G);

            % leaf indices
            idxl = tree.idx{tree.lvl+1};

            row_bases = cell(2^tree.lvl, 1);
            for b=1:2^tree.lvl
                % extract row sample for block
                Zb = Z(range(idxl{b}),:);
                if verbose; fprintf("row box %i of level %i, column box %i of level %i, extracting indices %i-%i\n", 1, 0, b, tree.lvl, idxl{b}); end

                % compute row basis
                row_bases{b} = col_basis(Zb, r);
            end

            % write row bases
            A.V = BDMatrix(row_bases);

            %% compute factors W
            % initialize factors
            A.W = cell(tree.lvl,1);

            for l=1:tree.lvl
                if verbose; fprintf("\n--- Level %i ---\n", l); end
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
                        if verbose; fprintf("\trow box %i of level %i, applying row box %i of level %i\n", c, l, a, m); end
                    end
                    for b=1:2^(tree.lvl-l)
                        % extract row sample for block
                        Zb = Z(((b-1)*2*r+1):(b*2*r),:);
                        if verbose; fprintf("row box %i of level %i, column box %i of level %i, extracting indices %i-%i\n", c, l, b, tree.lvl-l+1, (b-1)*2*r+1, b*2*r); end

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

            %% compute column bases F
            if verbose; fprintf("\n--- Level %i ---\n", tree.lvl+1); end
            col_bases = cell(2^tree.lvl, 1);
            % sample standard Gaussian random matrix
            G = randn(n,s);
            % compute forward solves of random matrix
            Z = fwd(G);
            % apply row bases (this time to G because it's a forward solve)
            G = A.V'*G;
            for b=1:2^tree.lvl
                % make a copy so we can apply row bases to G at each level
                Gb = G;
                for m=1:tree.lvl
                    % apply factors
                    a = ancestor(tree, m, b, tree.lvl);
                    Gb = A.W{m}.D{fix((a-1)/2)+1}.M{rem(a-1,2)+1} * Gb;
                    if verbose; fprintf("\trow box %i of level %i, applying row box %i of level %i\n", b, tree.lvl+1, a, m); end
                end
                
                % extract row sample for block
                Zb = Z(range(idxl{b}),:);
                if verbose; fprintf("row box %i of level %i, column box %i of level %i, extracting indices %i-%i\n", b, tree.lvl+1, b, 1, idxl{b}); end

                % compute column basis from sample by solving least squares problem 
                col_bases{b} = Zb * pinv(Gb);
            end

            % write column bases
            A.F = BDMatrix(col_bases);
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
    end
end