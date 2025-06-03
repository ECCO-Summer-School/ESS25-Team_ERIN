function [figure1] = Plot_SDSL_Trend(SDSLtrend, Lcg, roma, titlestr)
% plot sterodynamic sea-level trend on a global Robinson projection
%   - SDSLtrend: (n x 1) linear trends in mm/yr
%   - Lcg: (n x 2) longitudes and latitudes of grid
%   - roma: (256 x 3) colormap matrix
%   - titlestr: figure title

% if nargin < 4
%     titlestr = 'SDSL [1992–2017]';
% end

% 1x1 world grid
lon = (-179.5:179.5)';
lat = (-89.5:89.5)';
[latg, long] = meshgrid(lat, lon);
L1 = [long(:), latg(:)];

idx = knnsearch(Lcg, L1);
TrendField = SDSLtrend(idx);

for i = 1:numel(idx)
    [d1km(i,1), ~] = lldistkm(L1(i,:), Lcg(idx(i),:));
end
TrendField(d1km > 120) = NaN;
TrendField = reshape(TrendField, size(latg));

figure1 = figure('Color','w');
axes1 = axes('Parent',figure1, 'Position',[0.1 0.1 0.8 0.8]);
m_proj('robinson','lon',[-180 180]);
m_pcolor(long, latg, TrendField); shading flat;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('XTick',[], 'YTick', []);
caxis([-3, 3]);
colormap(roma);
cb = colorbar('FontSize', 10, 'Location', 'eastoutside');
cb.Ticks = -3:1:3;
cb.TickLabels = {'<-3', '-2', '-1', '0', '1', '2', '>3'};
set(cb.Label, 'String', 'Trend [mm yr^-^1]', 'rotation', 90)
% title(titlestr, 'FontSize', 12, 'FontWeight', 'bold');

end

