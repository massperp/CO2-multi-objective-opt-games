function keep = is_nondominated(F)

% Check Pareto dominance (minimization)
%   keep(i) is true if no row j of F dominates row i, i.e. there is no j with
%   F(j,:) <= F(i,:) in every objective and F(j,:) < F(i,:) in at least one.
%   Exact duplicates are kept only once.

    n = size(F, 1);
    keep = true(n, 1);
    for i = 1:n
        dominated_by = all(F <= F(i,:), 2) & any(F < F(i,:), 2);
        keep(i) = ~any(dominated_by);
    end
    [~, first] = unique(F, 'rows', 'stable');
    unique_mask = false(n, 1);
    unique_mask(first) = true;
    keep = keep & unique_mask;
end
