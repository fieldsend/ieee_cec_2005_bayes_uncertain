function p = p_dominates(u,v,u_var,v_var,n_u,n_v)

% p = p_dominates(u,v,u_var,v_var,n_u,n_v)
%
% Probability of domination calculation, assumes minimisation
% Returns probability that u dominates v
%
% u = mean_objective vector
% v = mean_objective vector
% u_var = variance of objective vectors associated with u
% v_var = variance of objective vectors associated with u
% n_u = number of samples used to calculate u and u_var
% n_v = number of samples used to calculate v and v_var
%
% Jonathan Fieldsend, University of Exeter

m = (v-u)./(((v_var)/n_v+(u_var)/n_u).^0.5); %vector of m terms
p = 0.5*(1+erf(m./sqrt(2)));
p = prod(p);
end
