function I = nearestIdx(aRef, aTest)

if size(aTest,2) > 1
    aTest = aTest.';
end

if size(aRef,2) == 1
    aRef = aRef.';
end

d = nan(numel(aTest), 1);
idx = nan(numel(aTest), 1);

edges = [-Inf, mean([aRef(2:end); aRef(1:end-1)]), +Inf];
I = discretize(aTest, edges);

end