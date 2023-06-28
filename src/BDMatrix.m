classdef BDMatrix
    % BDMatrix Block Diagonal Matrix
    % Simple implementation of block diagonal matrices.
    properties
        idx1
        idx2
        D
    end
    methods
        function A = BDMatrix(blocks)
            A.D = blocks;
            A.idx1 = cell(length(blocks),1);
            A.idx2 = cell(length(blocks),1);

            i1 = 1; i2 = 1;
            for b=1:length(blocks)
                A.idx1{b} = [i1, i1+size(blocks{b},1)-1];
                A.idx2{b} = [i2, i2+size(blocks{b},2)-1];
                i1 = i1 + size(blocks{b},1);
                i2 = i2 + size(blocks{b},2);
            end
        end
        function s = size(A, dim)
            iend = A.idx1{end};
            jend = A.idx2{end};
            if nargin > 1
                if dim==1
                    s = iend(2);
                elseif dim==2
                    s = jend(2);
                else
                    error("dimension must be either 1 or 2")
                end
            else
                s = [iend(2), jend(2)];
            end
        end
        function At = ctranspose(A)
            Dt = cell(length(A.D),1);
            for b=1:length(A.D)
                Dt{b} = A.D{b}';
            end
            At = BDMatrix(Dt);
        end
        function B = mtimes(A, X)
            m = size(A, 1);
            k = size(X, 2);
            B = zeros(m, k);
            for b=1:length(A.D)
                B(range(A.idx1{b}),:) = A.D{b}*X(range(A.idx2{b}),:);
            end
        end
        function S = sparse(A)
            % TODO: fix this to use the fast version (theres an indexing problem)
%             if all(cellfun(@(block) isa(block, 'double'), A.D))
%                 nnz = sum(cellfun(@numel, A.D));
%                 I = zeros(nnz,1); J = zeros(nnz,1); V = zeros(nnz,1); 
%                 i0 = 1;
%                 for b=1:length(A.D)
%                     nnzb = numel(A.D{b});
%                     [Ib, Jb] = meshgrid(range(A.idx1{b}), range(A.idx2{b}));
%                     I(i0:(i0+nnzb-1)) = reshape(Ib, 1, []);
%                     J(i0:(i0+nnzb-1)) = reshape(Jb, 1, []);
%                     V(i0:(i0+nnzb-1)) = reshape(A.D{b}, 1, []);
%                     i0 = i0 + nnzb;
%                 end
%                 S = sparse(I, J, V);
%             else 
                S = sparse(size(A,1), size(A,2));
                for b=1:length(A.D)
                    S(range(A.idx1{b}), range(A.idx2{b})) = sparse(A.D{b});
                end
%             end
        end
        function M = dense(A)
            M = A * eye(size(A));
        end
    end
end