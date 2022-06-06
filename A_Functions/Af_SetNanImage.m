function Af_SetNanImage(h)
% set the color of nans to black

nanclr = [0,0,0];

a = h.CData;
cm = colormap(gcf);

% find minimum and maximum
amin=min(a(:));
amax=max(a(:));
% size of colormap
n = size(cm,1);
% color step
dmap=(amax-amin)/n;

% standard imagesc
% him = imagesc(a);
% add nan color to colormap
colormap([nanclr; cm]);
% changing color limits
caxis([amin-5*dmap amax]);
% place a colorbar
hcb = colorbar;
% change Y limit for colorbar to avoid showing NaN color
ylim(hcb,[amin amax])

if nargout > 0
    h = him;
end