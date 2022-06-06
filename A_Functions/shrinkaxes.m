function shrinkaxes(scalefactor)
if(nargin<1)
   scalefactor = .1;
end
pos = get(gca, 'Position');
pos(2) = pos(2)+scalefactor*pos(4);
pos(4) = (1-scalefactor)*pos(4);
set(gca, 'Position', pos)