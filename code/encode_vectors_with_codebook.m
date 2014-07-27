function [code_mat,dots,cosines,dists] = encode_vectors_with_codebook(vectors,codebook,method,param)

[dim,n_vecs]    = size(vectors);
[dim2,k]        = size(codebook);
if dim ~= dim2
    error('vectors and codebook must have the same first dimension');
end

vector_norms    = sqrt(sum(vectors.^2));
codebook_norms  = sqrt(sum(codebook.^2));

dots            = codebook' * vectors;
cosines         = dots ./ (codebook_norms' * vector_norms);
dists           = (codebook_norms'.^2*ones(1,n_vecs)) + (ones(k,1)*vector_norms.^2) - dots;

switch method
    case 'cosine'
        absCosines  = abs(cosines);
        small_val   = absCosines < param;
        absCosines(small_val)   = 0;
        absCosines(~small_val)  = absCosines(~small_val) - param;
        code_mat    = sign(cosines) .* absCosines;
    case 'vq'
        code_mat    = zeros(k,n_vecs);
        for vi = 1:n_vecs
            [vals,inds]     = sort(dists(:,vi),'ascend');
            mark_inds       = inds(1:param);
            code_mat(mark_inds,vi)  = 1;
        end
end

end