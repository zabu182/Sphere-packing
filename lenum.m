function [v_short, v_idx] = lenum(mtx_M) 
%% LENUM Find the shortest vector of a lattice using exhaustive enumeration
% Exact answer using enumeration (slow for big mtx_M).
% 
% Inputs: 
%   mtx_M = real matrix, columns define a lattice. Must not have any zero
%   singular values. If it does, either mtx_M does not define a discrete
%   lattice or it contains redundant columns.
% Outputs:
%   v_short = shortest vector 
%   v_idx = index of v_short with respect to the basis passed in.
%           v_short = mtx_M*v_idx

% Nearly verbatim implementation of `ENUM` from:
%{
@article{schnorr1994lattice,
  title={Lattice basis reduction: Improved practical algorithms and solving subset sum problems},
  author={Schnorr, Claus-Peter and Euchner, Martin},
  journal={Mathematical programming},
  volume={66},
  number={1-3},
  pages={181--199},
  year={1994},
  publisher={Springer}
}
%}

n_sp_dim = size(mtx_M,1);
n_l_dim = size(mtx_M,2);

% Initialize

% Gram-Schmidt coefficients
mtx_Mstar = zeros(size(mtx_M));
mtx_mu = zeros(n_l_dim);
v_c = zeros(n_l_dim,1);
for idx=1:n_l_dim
    
    for jidx=1:(idx-1)
        mtx_mu(idx,jidx) = mtx_M(:,idx)'*mtx_Mstar(:,jidx)/v_c(jidx);
    end
    
    mtx_muscale = repmat(mtx_mu(idx,1:(idx-1)), n_sp_dim,1);
    mtx_Mstar(:,idx) = mtx_M(:,idx) - sum(mtx_muscale.*mtx_Mstar(:,1:(idx-1)),2);
    v_c(idx) = norm(mtx_Mstar(:,idx))^2;
end
for idx=1:n_l_dim
    mtx_mu(idx,idx) = 1;
end

% Main Loop
% 1.
d_cbar = v_c(1); 
v_utilde = zeros(n_l_dim+1,1); v_utilde(1) = 1; 
v_u = zeros(n_l_dim+1,1); v_u(1) = 1;
v_y = zeros(n_l_dim+1,1); v_y(1) = 0;
v_Ddelta = zeros(n_l_dim+1,1); v_Ddelta(1) = 0;
v_delta = zeros(n_l_dim+1,1); v_delta(1) = 1;
idx_s = 1;
idx_t = 1;

v_v = zeros(n_l_dim+1,1);

v_ctilde = zeros(n_l_dim+1,1);
for idx_i=2:(n_l_dim+1)
    v_ctilde(idx_i) = 0;
    v_u(idx_i) = 0;
    v_utilde(idx_i) = 0;
    v_y(idx_i) = 0;
    v_Ddelta(idx_i) = 0;
    v_delta(idx_i) = 1;
end

% 2.
while(idx_t <= n_l_dim)
    v_ctilde(idx_t) = v_ctilde(idx_t+1) + (v_y(idx_t)+v_utilde(idx_t))^2*v_c(idx_t);
    if( v_ctilde(idx_t)<d_cbar )
        if(idx_t > 1)
            idx_t = idx_t-1;
            v_y(idx_t) = sum(v_utilde((idx_t+1):idx_s).* ...
                mtx_mu((idx_t+1):idx_s,idx_t));
            v_utilde(idx_t) = round(-v_y(idx_t));
            v_v(idx_t) = round(-v_y(idx_t));
            v_Ddelta(idx_t) = 0;
            if(v_utilde(idx_t) > -v_y(idx_t))
                v_delta(idx_t) = -1;
            else
                v_delta(idx_t) = 1;
            end
        else
            d_cbar = v_ctilde(1);
            v_u = v_utilde;
        end
    else
        idx_t = idx_t+1;
        idx_s = max(idx_s, idx_t);
        if(idx_t < idx_s)
            v_Ddelta(idx_t) = -v_Ddelta(idx_t);
        end
        if(v_Ddelta(idx_t)*v_delta(idx_t) >= 0) 
            v_Ddelta(idx_t) = v_Ddelta(idx_t) + v_delta(idx_t);
        end
        v_utilde(idx_t) = v_v(idx_t) + v_Ddelta(idx_t);
    end
end

v_idx = v_u(1:(end-1));
v_short = mtx_M*v_idx;
%v_test = min(vecnorm(mtx_M))-norm(v_short);
    
return;
end