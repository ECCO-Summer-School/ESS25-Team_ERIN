% steric contribution
clc; clear all; close all;
%% load data from Dangendorf et al. (2024) GMSL reconstruction
cd ("C:\Users\erobson\OneDrive - Tulane University\phd_work\code");

addpath('D24_data\'); % data
% toolboxes
addpath('m_map\');
addpath('wafo-master\');
addpath('cdt\');

load KSSLfin.mat; % reconstruction
load KSSL_SEfin.mat; % associated standard error
load regTG_coords.mat; % tide gauge data

clearvars -except SDSL tt LALT regionSL regTG_coords

% sterodynamic SL component 
s = find(tt>=1992 & tt <=2017); % match same timespan as ECCO
t = tt(s);
SDSL = SDSL(:,s)*1000; % mm


%% compute linear trend at each grid point
SDSLtrend = NaN(size(SDSL, 1), 1);
x = [t(:) ones(length(t), 1)];
for i = 1:size(SDSL, 1)
    y = SDSL(i, :)';
    if sum(~isnan(y)) >=20 % require 20 years of data 
        m = x(~isnan(y), :) \ y(~isnan(y));
        SDSLtrend(i) = m(1); % mm yr-1
    end
end

% plot D24 sdsl trends (1992-2017)
load roma.mat
[figure1] = Plot_SDSL_Trend(SDSLtrend, LALT, roma, 'Dangendorf et al., 2024 SDSL trend [1992-2017]');
%% 
cd ('C:\Users\erobson\OneDrive - Tulane University\ECCO\mat_files')
%% ECCO v4r4 output arrays
% GRID = v4 grid
% % grid variables 
% RAC = surface area (m^2) of each grid cell
% MSK = land/sea mask
% XC = longitude of grid cell
% YC = latutude of grid cell

% SSHDYN = monthly ocean-dynamic SL (1992-2017)

load('GRID.mat');
load('SSHDYN.mat');

% figure('color','white')
% pcolor(std(SSHDYN,0,3)'); shading flat; colorbar
% caxis([0 0.1]);
% title('Standard deviation of monthly dynamic sea level 1992-2017 [m]')
% xlabel('Longitude'); ylabel('Latitude')
%% calculate steric component from ECCO
% SSHDYN - OBPNOPAB
% OBPNOPAB = manometric sea level - IBE
% SSHDYN = absolute sea level (would compare with altimtery)

load("OBPNOPAB.mat");
stericSL = SSHDYN - OBPNOPAB; % ecco ssl
stericSL = stericSL*1000; % mm
save('stericSL.mat', 'stericSL');

% ECCO = monthly 
% convert monthly to annual data
yrs = 1992:2017;
nyrs = numel(yrs);
stericSL_annual = NaN(size(stericSL, 1), size(stericSL, 2), numel(yrs));

for ii = 1:nyrs
    m_idx = (1:12) + (ii-1)*12;
    stericSL_annual(:, :, ii) = nanmean(stericSL(:, :, m_idx), 3);
end
%% flatten ecco grid

ecco_lon = XC(:);
ecco_lat = YC(:);
mask = MSK(:) == 1;  % ocean = 1
ecco_lon = ecco_lon(mask);
ecco_lat = ecco_lat(mask);

% regrid ECCO onto D24
ECCO_newgrid = NaN(size(LALT, 1), numel(yrs));

for yr = 1:numel(yrs)
    SSL = stericSL_annual(:,:,yr);
    SSL_flat = SSL(:);
    SSL_flat = SSL_flat(mask);
    Fyr = scatteredInterpolant(ecco_lon, ecco_lat, SSL_flat, 'linear', 'none');
    ECCO_newgrid(:, yr) = Fyr(LALT(:, 1), LALT(:, 2));
end
%% linear trend and plot ECCO onto D24 grid

stericSL_trend = NaN(size(ECCO_newgrid,1), 1);
x = [yrs(:), ones(numel(yrs),1)];

for i = 1:size(ECCO_newgrid, 1)
    y = ECCO_newgrid(i, :)';
    if sum(~isnan(y)) >= 20
        m = x(~isnan(y), :) \ y(~isnan(y));
        stericSL_trend(i) = m(1);
    end
end

cd ("C:\Users\erobson\OneDrive - Tulane University\phd_work\code");
[figure2] = Plot_SDSL_Trend(stericSL_trend, LALT, roma, 'ECCO SDSL trend [1992–2017]');
%% 

% % mask land points (for 1 degree grid from Chris [360 x 360])
% stericSL(repmat(MSK==0, 1, 1, size(stericSL,3))) = NaN; % MSK = 1 is ocean, MSK = 0 is land
% % plot this
% figure('color','white')
% pcolor(MSK'); shading flat; colorbar
% % title('ECCO Grid Mask')
%% temporal correlation

ECCO = NaN(size(SDSL));

for yr = 1:26
    SSL = stericSL_annual(:, :, yr);
    SSL_flat = SSL(:);
    SSL_flat = SSL_flat(mask); % keep ocean-only values
    Fyr = scatteredInterpolant(ecco_lon, ecco_lat, SSL_flat, 'linear', 'none');
    ECCO(:, yr) = Fyr(LALT(:, 1), LALT(:, 2));
end

dtcc = NaN(size(SDSL, 1), 1);

for i = 1:size(SDSL, 1)
    temp1 = SDSL(i, :);
    temp2 = ECCO(i, :);
    valid = ~isnan(temp1) & ~isnan(temp2);
        if sum(valid) >=10
            temp1_dm = temp1(valid) - mean(temp1(valid));
            temp2_dm = temp2(valid) - mean(temp2(valid));
            [dtcc(i)] = corr(detrend(temp1_dm'), detrend(temp2_dm'));
        end
end
%% map temporal correlation

lon = (-179.5:179.5)';
lat = (-89.5:89.5)';
[latg, long] = meshgrid(lat, lon);
L1 = [long(:), latg(:)];

idx = knnsearch(LALT, L1);
dtcc_field = dtcc(idx);

% mask points > 120 km away
d1km = NaN(numel(idx), 1);
for i = 1:numel(idx)
    [d1km(i), ~] = lldistkm(L1(i, :), LALT(idx(i), :));
end
dtcc_field(d1km > 120) = NaN;
dtcc_field = reshape(dtcc_field, size(latg));

figure('color','w');
m_proj('robinson','lon',[-180 180]);
m_pcolor(long, latg, dtcc_field); shading flat;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('XTick',[],'YTick',[]);
caxis([-1 1]);
colormap(roma);
cb = colorbar('FontSize', 10, 'Location', 'eastoutside');
set(cb, 'Ticks', [-1 0 1]);
set(cb.Label, 'String', 'Correlation Coefficient', 'rotation', 90);
%% pattern correlation
valid = ~isnan(SDSLtrend) & ~isnan(stericSL_trend);
cctrend = corr(SDSLtrend(valid), stericSL_trend(valid));
%% trend difference
difference = stericSL_trend - SDSLtrend; % ECCO - D24
cd ('C:\Users\erobson\OneDrive - Tulane University\ECCO\mat_files')

% figure
lon = (-179.5:179.5)';
lat = (-89.5:89.5)';
[latg, long] = meshgrid(lat, lon);
L1 = [long(:), latg(:)];
% nearest neighbour match ECCO to D24 grid
idx = knnsearch(LALT, L1);
TrendField = difference(idx);
% distance from matched grid point to actual data location
d1km = NaN(numel(idx),1);
for i = 1:numel(idx)
    [d1km(i), ~] = lldistkm(L1(i,:), LALT(idx(i),:));
end

TrendField(d1km > 120) = NaN;
TrendField = reshape(TrendField, size(latg));

figure('color','w');
axes1 = axes('position',[0.1 0.1 0.8 0.8]);
m_proj('robinson','lon',[-180 180]);
m_pcolor(long, latg, TrendField); shading flat;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('XTick',[], 'YTick', []);
caxis([-5, 5]);
colormap(roma);
cb = colorbar('FontSize', 10, 'Location', 'eastoutside');
set(cb, 'Ticks', [-5, 0, 5]);
set(cb, 'TickLabels', {'-5', '0', '5'});
set(cb.Label, 'String', 'Trend difference [mm yr^-^1]', 'rotation', 90);
%% RMSE

rmse_grid = NaN(size(SDSLtrend));

valid = ~isnan(SDSLtrend) & ~isnan(stericSL_trend);
rmse_grid(valid) = sqrt((SDSLtrend(valid) - stericSL_trend(valid)).^2);

% map
lon = (-179.5:179.5)';
lat = (-89.5:89.5)';
[latg, long] = meshgrid(lat, lon);
L1 = [long(:), latg(:)];

% nearest-neighbor match to LALT points
idx = knnsearch(LALT, L1);
rmse_field = rmse_grid(idx);
d1km = NaN(numel(idx),1);
for i = 1:numel(idx)
    [d1km(i), ~] = lldistkm(L1(i,:), LALT(idx(i),:));
end
rmse_field(d1km > 120) = NaN;
rmse_field = reshape(rmse_field, size(latg));


figure('color','w');
axes1 = axes('position',[0.1 0.1 0.8 0.8]);
m_proj('robinson','lon',[-180 180]);
m_pcolor(long, latg, rmse_field); shading flat;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('XTick',[], 'YTick', []);
colormap(roma);
cb = colorbar('FontSize', 10, 'Location', 'eastoutside');
caxis([0 10]);
set(cb, 'Ticks', 0:2:10);
set(cb.Label, 'String', 'RMSE [mm yr^{-1}]', 'rotation', 90);
%% divide map by region(basins from [Thompson and Merrifield, 2014)

regions = unique(regionSL);
regcoords = struct();

regn =  {
    'South Atlantic',            % 1
    'Indo-Pacific',              % 2
    'Northwest Pacific',         % 3
    'Subpolar North Atlantic',   % 4
    'Subtropical North Atlantic',% 5
    'East Pacific',              % 6
    'Arctic',                    % 7
    'Southern Ocean'             % 8
};

for i = 1:length(regions)
    cr = regions(i);
    idx = (regionSL == cr);

    regcoords(i).region = cr;
    regcoords(i).name = regn{cr};
    regcoords(i).lon = LALT(idx, 1);
    regcoords(i).lat = LALT(idx, 2);
end

regcoords = regcoords(ismember([regcoords.region], 1:6));

% plot regions on map
figure;
m_proj('robinson','lon',[-180 180],'lat',[-90 90]);

hold on
lgdhandles = gobjects(length(regcoords),1);
for i = 1:length(regcoords)
    lon = regcoords(i).lon;
    lat = regcoords(i).lat;

    [x, y] = m_ll2xy(lon, lat);

    lgdhandles(i) = scatter(x, y, 2, 'filled');
end

m_coast('patch',[1 1 1],'edgecolor',[0 0 0]);
m_grid('XTick',[],'YTick',[],'Fontsize',8);
legend(lgdhandles, {regcoords.name}, 'Location','eastoutside', 'Orientation','vertical', 'Box', 'off');
% add corresponding tide gauge locations
hold on;

colors = lines(6);
for i = 1:length(regTG_coords)
    data = regTG_coords(i).data;
    ntg = length(data);

    lon = zeros(ntg, 1);
    lat = zeros(ntg, 1);

    for j = 1:ntg
        lon(j) = data(j).lon;
        lat(j) = data(j).lat;
    end

    h(i) = m_plot(lon, lat, 'o', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 5);
end 
%% pattern correlation for each region

ccregion = struct();

for i = 1:length(regions)-2
    r = regions(i);
    idx = regionSL == r;

    temp1 = nanmean(SDSL(idx, :), 1);
    temp2 = nanmean(ECCO(idx, :), 1);

    valid = ~isnan(temp1) & ~isnan(temp2);
    temp1_dm = temp1(valid) - mean(temp1(valid));
    temp2_dm = temp2(valid) - mean(temp2(valid));

    regcc = corr(detrend(temp1_dm'), detrend(temp2_dm'));
    
    ccregion(i).name = regn{r};
    ccregion(i).r = regcc;

end
%% spatial-averaged stericSL by region

for i = 1:length(regions) - 2
    r = regions(i);
    idx = regionSL == r;

    % weight
    w = cosd(LALT(idx, 2));
    w = w / sum(w);
    temp3 = sum(SDSL(idx, :) .* w, 1);
    temp4 = sum(ECCO(idx, :) .* w, 1);

    valid = ~isnan(temp3) & ~isnan(temp4);
    dttemp3 = detrend(temp3(valid) - mean(temp3(valid)));
    dttemp4 = detrend(temp4(valid) - mean(temp4(valid)));

    [w_cc] = corr(dttemp3', dttemp4');

    ccregion(i).name = regn{r};
    ccregion(i).wr = w_cc;
end

%% compare SDSL at tide gauges with steric height from ECCO at similar locations

L = [];

for i = 1:length(regTG_coords)
    lon = [regTG_coords(i).data.lon];
    lat = [regTG_coords(i).data.lat];
    L = [L; lon(:), lat(:)];
end

idx = knnsearch(LALT, L); % nearest point in grid to TG
tg_SDSL = SDSL(idx, :)';
tg_ECCO = ECCO(idx, :)';
regionTG = regionSL(idx);

tgtrend_SDSL = NaN(size(tg_SDSL, 2), 3);
tgtrend_ECCO = NaN(size(tg_ECCO, 2), 3);
tg_cc = NaN(size(tg_SDSL, 2), 1);

t = tt(tt >=1992 & tt <= 2017);

for ii = 1:size(tg_SDSL, 2)
    temp5 = tg_SDSL(:, ii);
    temp6 = tg_ECCO(:, ii);
    valid = ~isnan(temp5) & ~isnan(temp6);

    if sum(valid) > 15
        tgtrend_SDSL(ii, :) = trend_new(1, t(valid), temp5(valid), 'N');
        tgtrend_ECCO(ii, :) = trend_new(1, t(valid), temp6(valid), 'N');

        tg_cc(ii) = corr(temp5(valid), temp6(valid));
    end
    tg_diff = tgtrend_ECCO(:,1) - tgtrend_SDSL(:,1);
end
%% 
% correlations
figure('Color', 'w');
m_proj('robinson', 'lon', [-180 180], 'lat', [-80 88]);
hold on;
m_coast('patch', [1 1 1], 'edgecolor', 'k');
valid_cc = ~isnan(tg_cc);
[x, y] = m_ll2xy(L(valid_cc, 1), L(valid_cc, 2));
scatter(x, y, 30, tg_cc(valid_cc), 'filled');
m_grid('fontsize', 10, 'XTick', [], 'YTick', []);
colormap(roma);
caxis([-1 1]);
cb = colorbar;
cb.Label.String = 'Correlation Coefficient';

% trend differences
figure('Color', 'w');
m_proj('robinson', 'lon', [-180 180], 'lat', [-80 88]);
hold on;
m_coast('patch', [1 1 1], 'edgecolor', 'k');
valid_diff = ~isnan(tg_diff);
[x, y] = m_ll2xy(L(valid_diff, 1), L(valid_diff, 2));
scatter(x, y, 30, tg_diff(valid_diff), 'filled');
m_grid('fontsize', 10, 'XTick', [], 'YTick', []);
colormap(roma);
caxis([-3 3]);
cb = colorbar;
cb.Label.String = '[mm yr^-^1]';















    







