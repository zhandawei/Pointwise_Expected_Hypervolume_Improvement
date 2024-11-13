function y = Infill_pEHVI(x,kriging_obj,f)
[num_f, m] = size(f);
num_x = size(x,1);
u = zeros(num_x,m);
s = zeros(num_x,m);
for ii = 1:m
    [u(:, ii),s(:, ii)] = Kriging_Predictor(x, kriging_obj{ii});
end
r = max(f,[],1) + 1;
y = zeros(num_x,1);
for ii = 1:num_x
    % if the solution is predicted as a dominated solution, 
    % we assign its metric as a negative value 
    if sum(sum(f <= (u(ii,:)),2) == m) > 0
        y(ii) = - max(prod(1 + max(u(ii,:)-f,0),2) - 1);
    else
        % else we calculate its pEHVI value
        u_matrix = repelem(u(ii,:),num_f,1);
        s_matrix = repelem(s(ii,:),num_f,1);
        Phi1 = (r-u_matrix).*gausscdf((r-u_matrix)./s_matrix) + s_matrix.*gausspdf((r-u_matrix)./s_matrix);
        Phi2 = (f-u_matrix).*gausscdf((f-u_matrix)./s_matrix) + s_matrix.*gausspdf((f-u_matrix)./s_matrix);
        ehvi = prod(Phi1,2) - prod(Phi1-Phi2,2);
        y(ii) = min(ehvi,[],1);
    end
end
end

function res=gausspdf(x)
res=1/sqrt(2*pi)*exp(-x.^2/2);
end

function y=gausscdf(x)
 y=0.5*(1+erf(x/sqrt(2)));
end


