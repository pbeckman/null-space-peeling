classdef HSSMatrix
    % HSSMatrix Hierarchically Semi-Separable Matrix
    % HSS factorization in "telescoping" form 
    % A = A_L = U_L A_{L-1} V_L^* + D_L with constructor  based on 
    % "Linear-complexity black-box randomized compression of hierarchically 
    % block separable matrices" by Levitt and Martinsson 2022.
    %
    % Each of the U, V, and D factors are block diagonal (implemented 
    % using BDMatrix).
    properties
        tree
        U
        V
        D
    end
    methods
        function A = HSSMatrix(fwd, adj, tree, r, s)
            A.tree = tree;

            % number of points n
            i0 = tree.idx{1}; b0 = i0{1}; n = b0(2);

            % sample standard Gaussian random matrices
            W = randn(n,s);
            S = randn(n,s);
            % compute forward and adjoint solves of random matrices
            Y = fwd(W);
            Z = adj(S);

            % initialize factors
            A.U = cell(tree.lvl+1,1); A.V = cell(tree.lvl+1,1); A.D = cell(tree.lvl+1,1);

            % box indices at leaf level
            idxl = tree.idx{tree.lvl+1};

            for l=tree.lvl:-1:0
                % initialize blocks at current level
                UB = cell(2^l,1); VB = cell(2^l,1); DB = cell(2^l,1);

                for b=1:2^l
                    % box indices
                    idxb = idxl{b};
                    
                    % extract column basis
                    Wb = W(range(idxb),:);
                    Yb = Y(range(idxb),:);
                    
                    % extract row basis
                    Sb = S(range(idxb),:);
                    Zb = Z(range(idxb),:);

                    if l > 0
                        % compute column null space and basis
                        Pb    = null_basis(Wb, r);
                        UB{b} = col_basis(Yb*Pb, r);

                        % compute row null space and basis
                        Qb    = null_basis(Sb, r);
                        VB{b} = col_basis(Zb*Qb, r);

                        % compute discrepancy
                        DB{b} = (Yb - UB{b}*(UB{b}'*Yb)) * pinv(Wb) + UB{b}*(UB{b}'*((Zb - VB{b}*(VB{b}'*Zb)) * pinv(Sb))');
                    else
                        % compute discrepancy
                        DB{b} = Yb * pinv(Wb);
                    end
                end

                % write block diagonal factors
                A.U{l+1} = BDMatrix(UB); A.V{l+1} = BDMatrix(VB); A.D{l+1} = BDMatrix(DB);

                if l > 0
                    % project onto columnspaces
                    Y = A.U{l+1}'*(Y - A.D{l+1}*W);
                    W = A.V{l+1}'*W;

                    % project onto rowspaces
                    Z = A.V{l+1}'*(Z - A.D{l+1}'*S);
                    S = A.U{l+1}'*S;

                    % box indices at next level
                    idxl = combine_idx(A.U{l+1}.idx2);
                end
            end
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
        function B = mtimes(A, X)
            Xs = cell(A.tree.lvl, 1);

            Xs{end} = A.V{end}'*X;
            for l=A.tree.lvl-1:-1:1
                Xs{l} = A.V{l+1}'*Xs{l+1};
            end

            B = A.D{1}*Xs{1};
            for l=1:A.tree.lvl-1
                B = A.U{l+1}*B + A.D{l+1}*Xs{l+1};
            end

            B = A.U{end}*B + A.D{end}*X;
        end
        function M = dense(A)
            M = A * eye(size(A, 2));
        end
    end
end