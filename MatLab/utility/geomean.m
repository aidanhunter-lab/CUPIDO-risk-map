function G = geomean(x)
% Geometric mean of x containing non-negative values -- accounts for zeros
if any(x < 0)
    warning('Negative values not allowed!')
    G = nan;
    return
end
w = x ~= 0; % index non-zero components
if isvector(x)
    if all(w == 0)
        G = 0;
        return
    end
    xp = x(w); % positive components
    f = sum(w) / length(x); % proportion [0,1] of positive components
    G = f * exp(sum(log(xp)) / length(xp));
else
    % If x is an array then find means along first dimension
    allZeros = all(x == 0);
    xs = size(x);
    y = x(:,:);
    w = w(:,:);
    ws = sum(w);
    f = ws ./ xs(1);
    y(~w) = nan;
    y = log(y);
    y = num2cell(y, 1);
    G = f .* exp(cell2mat(cellfun(@(z) sum(z, 'omitnan'), y, ...
        'UniformOutput', false)) ./ ws);
    G = reshape(G, [1, xs(2:end)]);
    G(allZeros) = 0;
end
