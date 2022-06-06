function setLineStyle(h,Color,marker,markerSize,lineSty,lineWidth,varargin)
if size(varargin) > 0
    if varargin{1} == 1
        outlineMarker = 1;
    else
        outlineMarker = 0;
    end
else
    outlineMarker = 0;
end


h.Color         = Color;
h.Marker        = marker;
h.MarkerSize    = markerSize;
h.LineStyle     = lineSty;
h.LineWidth     = lineWidth;

if ~strcmp(marker,'.')
    h.MarkerFaceColor = Color;
end

if outlineMarker
    h.MarkerEdgeColor = 'k';
end