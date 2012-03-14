function [PolicyInit Multipliers]=GetInitialApproxPolicy(xTarget,xStore,PolicyRulesStore,MultipliersStore)

%distance=(sum(((xStore-repmat(xTarget,length(xStore),1)).^2),1)).^.5;
distance=sum((xStore-repmat(xTarget,length(xStore),1)).^2,2);

[~, ref_id]=min(distance);
PolicyInit=PolicyRulesStore(ref_id,:);

end

