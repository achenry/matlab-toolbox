function [nVariation nPermutation Permutation] = CreatePermutationMatrix(Variation)
[nVariation nColumn]            = size(Variation); 

if nColumn~=2
    error('Wrong Size!')
end

for iVariation                  = 1:nVariation;
    VariationDepth(iVariation)  = length(Variation{iVariation,2}); 
end    
nPermutation                    = prod(VariationDepth);

Permutation     = [];
for iVariation  = 1:nVariation;        
	Permutation	= [reshapeRows(Permutation,VariationDepth(iVariation)) ....
        repmat([1:VariationDepth(iVariation) ]',prod(VariationDepth(1:iVariation-1)),1)];    
end    