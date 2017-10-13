function [A,Ao,Av,samples,s_var,An,Ap,H,Ho,Hv,Hn] = one_plus_one_uncertain_bes(sdc,evaluations,num_reevaluations,cost_function,domain_function,l,num_obj,std_mut,func_arg,alpha,eta,mutate_all_values)

% [A,Ao,Av,samples,s_var,An,Ap] = 
% one_plus_one_uncertain_bes(sdc,evaluations,num_reevaluations,cost_function,domain_function,l,num_obj,std_mut,func_arg,alpha,eta,mutate_all_values)
%
%
% ARGUMENTS
% sdc = summation used in archiving check (1) or at element level (0)
% evaluations = number of total function evaluations to run optimiser for
% num_reevaluations = number of function evaluations per solution
% cost_function = objective function, must accept as arguments a real
%   vector of length l and a integer denoting the number of objectives, must
%   return a vector of 'number of objective' elements. All parameters are
%   on the range [0,1]
% l = number of parameters
% num_obj = number of objectives
% std_mut = standard deviation of Guassian mutation used
% func_arg = meta argument to objective function
% alpha = level of alpha for probabilistic dominance check 
% eta = discount rate for history when there is temporal noise
%   suspected (optional, set at 1.0, i.e. no discount, if argument
%   not provided) 
% mutate_all_values = 1 to mutate all values, 0 to mutate just one decision
%   variable at random (optional, set at 0, i.e. mutate just one variable 
%   if not passed)  
% RETURNS
% A = Archive solutions
% Ao = Corresponding archive objective evaluations
% Av = Corresponding varaince estimates
% samples = structure holding algorithm state, sampled every 500 iterations   
% s_var = matrix or corresponding sample variances  
% An = num reevaluations per archive member
% Ap =  highest level of probabilistic domination experienced by archive member
% H = all visited design vectors    
%
% Jonathan Fieldsend, University of Exeter, 2005,2009,2014
    
func_arg.dom_check=0; %move to point if accepted into archive
if nargin ==10
    fprintf('Eta unspecified, so set at default 1.0 (assumes stationary variance over time/space)\n');
    fprintf('Schedule unspecified, so assumed fixed at num_reevaluations value throughout optimisation\n');
    eta=1.0;
    mutate_all_values = 0;
else
    fprintf('eta of %f used',eta);
end


sample_rate=500;
sample_index=1;
s_var_index=1;
samples=cell(floor(evaluations/sample_rate),2);
H = zeros(ceil(evaluations/num_reevaluations),l);
Ho = zeros(ceil(evaluations/num_reevaluations),num_obj);
Hv = zeros(ceil(evaluations/num_reevaluations),num_obj);
Hn = zeros(ceil(evaluations/num_reevaluations),1);
s_var=zeros(evaluations/sample_rate,num_obj);
% Create random indiviual (Uniform) and evaluate
c = generate_legal_initial_solution(domain_function,l,func_arg);
a = 0.0; %uniformative prior on a
b = 0.0; %uniformative prior on b
[c_mean_objectives,c_est_var,a,b]=evaluate_f(cost_function,c,num_obj,num_reevaluations,func_arg,a,b,eta);

% In the case of the first iteration, solutions with domination count 0
% are the estimated Pareto front
A = c; %pareto set
Ao = c_mean_objectives; %pareto front
Av = c_est_var; %estimated variance of archive members
An = num_reevaluations; %num revelauations per archive member
Ap = 0; % highest level of probabilistic domination experienced by archive member
c_num = num_reevaluations;

H(1,:)= c; %history of visited locations
Ho(1,:)= c_mean_objectives; %history of estimated locations
Hn(1) = num_reevaluations;
Hv(1,:) = c_est_var;
k = 1; % history index counter


% initalise evaluation counter and print counter
num_evaluations = num_reevaluations;
[samples,s_var,sample_index,s_var_index] = update_statistics(samples,s_var,sample_index,s_var_index,num_evaluations,sample_rate,Ao,A,c_est_var);
% Iterate until number of evalutions reaches maximium
counter = num_evaluations;
while (num_evaluations<evaluations)
    
    % iterate
    m = perturb(domain_function,c,l,std_mut,func_arg,mutate_all_values);
    [m_mean_objectives,m_est_var,a,b]=evaluate_f(cost_function,m,num_obj,num_reevaluations,func_arg,a,b,eta);
    m_num = num_reevaluations;
    %update history
    k = k+1;
    H(k,:) = m;
    Ho(k,:) = m_mean_objectives;
    Hv(k,:) = m_est_var;
    Hn(k) = num_reevaluations;
    if (eta==1) %if noise is stationary always use most recent noise variance estimate
        Av_const = repmat(m_est_var,size(Av,1),1); % copy variance matrix to include most recent estimate of var
        if (p_dominates(c_mean_objectives,m_mean_objectives,m_est_var,m_est_var,c_num,m_num)<=alpha) %if c probably dominates m greater than alpha level then discard
            [c,c_mean_objectives,c_est_var,c_num,A,Ao,Av,An,Ap]=update_archive(sdc,c,c_mean_objectives,m_est_var,c_num,m,m_mean_objectives,m_est_var,m_num,A,Ao,Av_const,An,Ap,num_obj,alpha,func_arg);
        end
    else
        if (p_dominates(c_mean_objectives,m_mean_objectives,c_est_var,m_est_var,c_num,m_num)<=alpha) %if c probably dominates m greater than alpha level then discard
            [c,c_mean_objectives,c_est_var,c_num,A,Ao,Av,An,Ap]=update_archive(sdc,c,c_mean_objectives,c_est_var,c_num,m,m_mean_objectives,m_est_var,m_num,A,Ao,Av,An,Ap,num_obj,alpha,func_arg);
        end
    end
    if (func_arg.A_plus_1 == 1) %if (|A|+1) rather than (1+1) model used
        [c,c_mean_objectives,c_est_var,c_num] = A_plus_one_sample(A,Ao,Av,An,Ap,func_arg);
    end
    counter = counter +num_reevaluations;
    num_evaluations = num_evaluations+num_reevaluations;
    if  counter>=500
        fprintf('Evaluations %d, Archive size %d, alpha %f\n', num_evaluations, size(A,1), alpha);
        if (eta==1)
            fprintf('Eta equals 1')
        end
        counter = 0;
        [samples,s_var,sample_index,s_var_index] = update_statistics(samples,s_var,sample_index,s_var_index,num_evaluations,sample_rate,Ao,A,m_est_var);
    end
end

end

%--------------------------------------------------------------------------
function x = generate_legal_initial_solution(domain_function, l, func_arg)

%generate an initial legal solution
illegal =1;
while (illegal)
    x = rand(1,l);
    x = x.*func_arg.range+func_arg.lwb; %put on data range
    if( feval(domain_function,x,func_arg) ) %check if legal
        illegal = 0;
    end
end
end

%--------------------------------------------------------------------------
function [samples,s_var,sample_index,s_var_index] = update_statistics(samples,s_var,sample_index,s_var_index,num_evaluations,sample_rate,Ao,A,est_var)

% keep track of statistics
if (num_evaluations >= (sample_rate*sample_index))
    samples{sample_index,1}=Ao;
    samples{sample_index,2}=A;
    s_var(sample_index,:)=est_var;
    sample_index=sample_index+1; %increment counter
end
end

%--------------------------------------------------------------------------
function [c,c_objectives,c_var,c_num,A,Ao,Av,An,Ap]=update_archive(summation_domed_check, c,c_objectives,c_var,c_num,m,m_objectives,m_var,m_num,A,Ao,Av,An,Ap,num_obj,alpha,func_arg)

p_domed=0;
max_level =0;

if summation_domed_check % domination compared by summing probibility of domination across archive members
    sum_level=0;
    for i=1:size(A,1)
        level = p_dominates(Ao(i,:),m_objectives,Av(i,:),m_var,An(i),m_num);
        sum_level = sum_level +level;
        if (sum_level>alpha); %reject if sum of archive p_doms at alpha level
            p_domed=1;
            break;
        end
    end
    if (p_domed==0)
        
        Ao=[Ao; m_objectives];
        A=[A; m];
        Av=[Av; m_var];
        An=[An; m_num];
        Ap=[Ap; sum_level];
        I=[];
        for i=size(A,1)-1:-1:1
            sum_level =0;
            for j = size(A,1):-1:1 %compare every archive member to every other member
                if (i~=j) %don't compare to oneself
                    level = p_dominates(Ao(j,:),Ao(i,:),Av(j,:),Av(i,:),An(j),An(i));
                else
                    level =0;
                end
                sum_level = sum_level+level;
                if (sum_level > alpha)
                    %mark to remove any archive members dominated at the alpha level
                    I = [I i];
                    break;
                end
            end
            if (sum_level > Ap(i))
                Ap(i) =level; % update maximum domination level experienced by archive member
            end
        end
        Ao(I,:) = [];
        A(I,:) = [];
        Av(I,:) = [];
        An(I) = [];
        Ap(I) = [];
           
        c=m;
        c_var=m_var;
        c_objectives=m_objectives;
        c_num = m_num;
    end
else %domination compared at an individual rather than a set level
    for i=1:size(A,1)
        level = p_dominates(Ao(i,:),m_objectives,Av(i,:),m_var,An(i),m_num);
        if level> max_level
            max_level =level;
        end
        if (level>alpha); %sum_level or individual level?
            p_domed=1;
            break;
        end
    end
    if (p_domed==0)
        % not probabilistically dominated -- however as we are using a
        % (1+1)--ES, and checking individually rather than summing across
        % archive members, we actually only need store the leading edge of
        % solutions, which means checking if the entrant is non-dominated
        % by the archive
        domed =0;
        for i=1:size(A,1)
            if (dominates(Ao(i,:),m_objectives))
                domed =1;
                break
            end
        end
        if (domed==0) % not dominated, so remove and leading edge members
            % and insert into archive
            to_remove=[];
            for i=1:size(A,1)
                if (dominates(m_objectives, Ao(i,:)))
                    to_remove = [to_remove i];
                end
            end
            % remove domed
            Ao(to_remove,:)=[];
            A(to_remove,:)=[];
            Av(to_remove,:)=[];
            An(to_remove) = [];
            % add into archive
            Ao=[Ao; m_objectives];
            A=[A; m];
            Av=[Av; m_var];
            An=[An; m_num];
            Ap = [Ap; max_level];
        end
        % set new location to be perturbed location which isn't p_domed by
        % archive
        c=m;
        c_var=m_var;
        c_objectives=m_objectives;
        c_num = m_num;
        %end
    end
end
end

%--------------------------------------------------------------------------
function p=m_dominates(Ao,x)

% compare all members of Ao to x for dominance, assumes one copy of x
% resides in Ao
[n,m] = size(Ao);
p=zeros(n,1);
for i=1:m
    p = p+ Ao(:,i)<=x(i);
end

p = (p==m);
p = (sum(p) > 1); % larger than one to discount itself at the end as will be compared
end

%--------------------------------------------------------------------------
function d=dominates(x,xp)

% compare x to xp for dominance

[~,m] = size(x);
wdom = sum(x<=xp);
dom = sum(x<xp);

d = (wdom==m) & (dom>0);
end

%--------------------------------------------------------------------------
function [c_mean_objectives,c_estimated_var,a,b]=evaluate_f(cost_function,c,num_obj,num_reevaluations,func_arg,a,b,eta)

% repeatedly evaluate the solution 'c'
% calculate the mean objective value and estimate the variance

results=zeros(num_reevaluations,num_obj); %preallocate matrix for efficiency
for i=1:num_reevaluations
    results(i,:)=feval(cost_function,c,num_obj,func_arg);
end

c_mean_objectives=mean(results);
b=eta*b+0.5*sum((results-repmat(c_mean_objectives,num_reevaluations,1)).^2); %bprime=eta*b+0.5(n-1)S^2
a=eta*a+0.5*(num_reevaluations-1);  %aprime = eta*a+0.5(n-1)
c_estimated_var=b./(a-1);
end
%--------------------------------------------------------------------------
function c = perturb(domain_function, c,l,std_mut,func_arg,mutate_all_values)

if (mutate_all_values==1)
    p=randn(size(c))*std_mut+c;
    while (feval(domain_function,p,func_arg)==0)    %ensure in valid range
        p=randn(size(c))*std_mut+c; %mutate with additive Gaussian noise
    end
    c=p;
else
    %Perturbs a single parameter of a solution vector
    I=randperm(l);
    i=I(1); %select a decision variable at random
    p = c;
    p(i) = randn()*std_mut+c(i); %mutate with additive Gaussian noise
    while (feval(domain_function,p,func_arg)==0)    %ensure in valid range
        p(i)=randn()*std_mut+c(i); %mutate with additive Gaussian noise
    end
    c=p;
end

end
%--------------------------------------------------------------------------
function [c,c_mean_objectives,c_est_var,c_num] = A_plus_one_sample(A,Ao,Av,An,Ap,func_arg)

% pseudo-random sampling from archive

% select dimension at random
[n,d]=size(Ao);

dimension_index = randperm(d);
dimension_index = dimension_index(1);

selected = 0;
while selected==0
    mx = max(Ao(:,dimension_index));
    mn = min(Ao(:,dimension_index));
    if (mx ==mn) % no range on data,
        % take first point if n ==1
        if n==1
            index = 1;
        else
            I=1:n;
            %multiple points in the bin, so roulette wheel selection based
            %on the probability of dominance
            
            [~, II] = min(Ap(I));
            index = I(II);
        end
        selected=1;
        if (func_arg.dom_check==1)
            if (m_dominates(Ao,Ao(index,:))==1)
                selected=0;
            end
        end
    else
        r = rand()*(mx-mn)+mn; %uniform point in range of data
        dd = abs(Ao(:,dimension_index)-r)/(mx-mn); %normalise distances to (0,1) range
        I = find (dd<=0.05); % get indices of thouse points within 5% either side of the sample
        
        if (length(I)==1) %only one point in the bin, so chose that
            index = I(1);
            selected=1;
            if (func_arg.dom_check==1)
                if (m_dominates(Ao,Ao(index,:))==1)
                    selected=0;
                end
            end
        elseif (length(I)>1)
            %multiple points in the bin, so roulette wheel selection based
            %on the probability of dominance
            
            [~, II] = min(Ap(I));
            index = I(II);
            selected=1;
            if (func_arg.dom_check==1)
                if (m_dominates(Ao,Ao(index,:))==1)
                    selected=0;
                end
            end
        end
    end
end

c = A(index,:);
c_mean_objectives = Ao(index,:);
c_est_var = Av(index,:);
c_num = An(index,:);

end
