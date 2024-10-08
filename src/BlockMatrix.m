classdef BlockMatrix
    % BlockMatrix Block Matrix
    % Simple implementation of block matrices.
    properties
        idx1
        idx2
        M
    end
    methods
        function A = BlockMatrix(blocks)
            A.M = blocks;
            A.idx1 = cell(size(blocks,1),1);
            A.idx2 = cell(size(blocks,2),1);

            i0 = 1;
            for i=1:size(blocks,1)
                heights    = cellfun(@(b) size(b,1), blocks(i,:));
                row_height = max(heights);
                if min(heights(heights > 0)) ~= row_height
                    error("Blocks in row %i have incompatible heights - cannot construct block matrix.", i)
                end
                A.idx1{i} = [i0, i0+row_height-1];
                i0 = i0 + row_height;
            end

            j0 = 1;
            for j=1:size(blocks,2)
                widths    = cellfun(@(b) size(b,2), blocks(:,j));
                row_width = max(widths);
                if min(widths(widths > 0)) ~= row_width
                    error("Blocks in column %i have incompatible wdiths - cannot construct block matrix.", j)
                end
                A.idx2{j} = [j0, j0+row_width-1];
                j0 = j0 + row_width;
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
            Mt = cell(flip(size(A.M)));
            for i=1:size(A.M, 1)
                for j=1:size(A.M, 2)
                    Mt{j,i} = A.M{i,j}';
                end
            end
            At = BlockMatrix(Mt);
        end
        function B = mtimes(A, X)
            m = size(A, 1);
            k = size(X, 2);
            B = zeros(m, k);
            for i=1:size(A.M, 1)
                for j=1:size(A.M, 2)
                    if all(size(A.M{i,j}) > 0)
                        B(range(A.idx1{i}),:) = B(range(A.idx1{i}),:) + A.M{i,j}*X(range(A.idx2{j}),:);
                    end
                end
            end
        end
        function S = sparse(A) % TODO: fix this to use the fast version
            S = sparse(size(A,1), size(A,2));
            for i=1:size(A.M, 1)
                for j=1:size(A.M, 2)
                    if all(size(A.M{i,j}) > 0)
                        S(range(A.idx1{i}), range(A.idx2{j})) = sparse(A.M{i,j});
                    end
                end
            end
%             C = cellfun(@sparse, A.M, 'UniformOutput', false);
%             S = cat(2, cat(1, C{:,:}));
        end
        function M = dense(A)
            M = A * eye(size(A, 2));
        end
    end
end