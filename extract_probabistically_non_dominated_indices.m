function [indices] = extract_probabistically_non_dominated_indices(Ao,Ho,An,Hn,Av,Hv,alpha)

p_non_dom = zeros(size(Ho,1),1); % if probabilistically non-dominated

for i=1:size(Ho,1)
    p_domed = 0;
    for k=1:size(Ao,1)
        p = p_dominates(Ao(k,:),Ho(i,:),Av(k,:),Hv(i,:),An(k),Hn(i));
        if (p>alpha)
            p_domed =1;
            break;
        end
    end
    if (p_domed==0)
        p_non_dom(i)=1;
    end
end
indices = find(p_non_dom==1);