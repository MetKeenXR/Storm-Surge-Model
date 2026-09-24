clear all; close all; clc;

disp('================================================================');
disp('   HONG KONG ADCIRC MESH - SMOOTH ELLIPTICAL REFINEMENT         ');
disp('   Resolution: 0.08-0.6km (HK Coast) → 18km (Open Ocean)       ');
disp('================================================================');
disp(' ');

%% ========================================================================
%  STEP 1: INITIALIZE OCEANMESH2D
% ========================================================================
disp('[1/7] Initializing OceanMesh2D...');

% Update these paths for your computer
oceanmesh_path = '/Users/beekenshrestha/Downloads/UST_Doc/HK_Surge/OceanMesh2D/';
run([oceanmesh_path,'setup_oceanmesh2d.m']);

disp('✅ OceanMesh2D initialized successfully.');

%% ========================================================================
%  STEP 2: SETUP PATHS AND PARAMETERS
% ========================================================================
disp('[2/7] Setting up paths and parameters...');

base_dir   = '/Users/beekenshrestha/Downloads/UST_Doc/HK_Surge/';
gebco_file = [base_dir,...
    'GEBCO_01_Jul_2026_caeca1e77c26/gebco_2026_n30.0_s10.0_w105.0_e135.0.nc'];

output_dir = [base_dir,'output/'];

% Create output directory
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Domain definition
bbox = [110.0 120.0; 18.0 24.0];
hk_lon = 114.17; hk_lat = 22.32;

% Resolution parameters - OPTIMIZED
max_el = 18000;     % 18km maximum at open ocean
min_el = 80;        % 80m minimum near Hong Kong
coarse_base = 3000; % 3km base (BETTER than 5km)
grade = 0.05;       % VERY LOW = VERY SMOOTH transition
fs = 2.5;           % Feature size for coastline

disp('✅ Parameters set:');
fprintf('   Domain: %.1f°E - %.1f°E, %.1f°N - %.1f°N\n', bbox(1,1), bbox(1,2), bbox(2,1), bbox(2,2));
fprintf('   HK Resolution: %.0fm (HK Coast) → %.0fm (Open Ocean)\n', min_el, max_el);
fprintf('   Base resolution: %.0fm (3km for other areas)\n', coarse_base);
fprintf('   Grade: %.2f (VERY SMOOTH transition)\n', grade);

%% ========================================================================
%  STEP 3: CREATE GEODATA OBJECT
% ========================================================================
disp('[3/9] Creating geodata object...');

% Load coastline shapefile
coastline_shp = 'GSHHS_f_L1';

try
    gdat = geodata('shp', coastline_shp, ...
                   'dem', gebco_file, ...
                   'bbox', bbox, ...
                   'h0', coarse_base);
    disp('✅ Geodata object created successfully.');
catch ME
    fprintf('⚠️ Warning: %s\n', ME.message);
    coastline_shp = 'GSHHS_l_L1';
    gdat = geodata('shp', coastline_shp, ...
                   'dem', gebco_file, ...
                   'bbox', bbox, ...
                   'h0', coarse_base);
    disp('✅ Geodata object created with fallback coastline.');
end

%% ========================================================================
%  STEP 4: CREATE SMOOTH ELLIPTICAL RESOLUTION FUNCTION
%  VERY GRADUAL transition with stretched ellipses - EXTENDED ZONES
% ========================================================================
disp('[4/9] Creating smooth elliptical resolution function...');

% Define the resolution function with VERY SMOOTH transition
function h = smooth_elliptical_resolution(x, y, hk_lon, hk_lat, min_el, max_el, coarse_base)
    % Calculate normalized distance in elliptical coordinates
    % Stretched ellipses aligned with coastline
    
    % Ellipse parameters - stretching for coastline alignment
    a = 1.5;  % Major axis (along coastline - NW-SE)
    b = 0.7;  % Minor axis (perpendicular to coastline)
    
    % Rotate ellipse to align with coastline (45 degrees for HK)
    theta = 45 * pi / 180;
    
    % Calculate coordinates relative to Hong Kong
    dx = x - hk_lon;
    dy = y - hk_lat;
    
    % Rotate coordinates to align with coastline
    x_rot = dx * cos(theta) + dy * sin(theta);
    y_rot = -dx * sin(theta) + dy * cos(theta);
    
    % Calculate elliptical distance (in km)
    dist_ellipse = sqrt((x_rot / a).^2 + (y_rot / b).^2) * 111;
    
    % VERY GRADUAL transition with EXTENDED ZONES
    % Start with base resolution
    h = coarse_base * ones(size(x));
    
    % Hong Kong refinement zone (0-250km) - EXTENDED
    idx_hk = dist_ellipse <= 250;
    dist_hk = dist_ellipse(idx_hk);
    
    % Zone 1: 0-20km - 80m to 400m (Hong Kong core)
    idx1 = dist_hk <= 20;
    h1 = 80 + (400-80) * (dist_hk(idx1) / 20);
    
    % Zone 2: 20-50km - 400m to 1km (Inner transition)
    idx2 = dist_hk > 20 & dist_hk <= 50;
    h2 = 400 + (1000-400) * ((dist_hk(idx2)-20) / 30);
    
    % Zone 3: 50-100km - 1km to 2.5km (Outer transition)
    idx3 = dist_hk > 50 & dist_hk <= 100;
    h3 = 1000 + (2500-1000) * ((dist_hk(idx3)-50) / 50);
    
    % Zone 4: 100-160km - 2.5km to 5km (Approach to base)
    idx4 = dist_hk > 100 & dist_hk <= 160;
    h4 = 2500 + (5000-2500) * ((dist_hk(idx4)-100) / 60);
    
    % Zone 5: 160-250km - 5km to 8km (Merging with base)
    idx5 = dist_hk > 160 & dist_hk <= 250;
    h5 = 5000 + (8000-5000) * ((dist_hk(idx5)-160) / 90);
    
    % Combine Hong Kong zones
    h_hk = zeros(size(dist_hk));
    h_hk(idx1) = h1;
    h_hk(idx2) = h2;
    h_hk(idx3) = h3;
    h_hk(idx4) = h4;
    h_hk(idx5) = h5;
    
    % Apply to the full array
    h(idx_hk) = h_hk;
    
    % Apply Gaussian smoothing for EXTREMELY SMOOTH transitions
    h = smoothdata(h, 'gaussian', 5);
    
    % Ensure bounds
    h = max(h, min_el);
    h = min(h, max_el);
end

% Create a fine grid for the resolution function
lon_vec = linspace(bbox(1,1), bbox(1,2), 500);
lat_vec = linspace(bbox(2,1), bbox(2,2), 500);
[LON, LAT] = meshgrid(lon_vec, lat_vec);

% Calculate resolution on the grid
RES = smooth_elliptical_resolution(LON, LAT, hk_lon, hk_lat, min_el, max_el, coarse_base);

% Create an interpolation function
F_res = scatteredInterpolant(LON(:), LAT(:), RES(:), 'linear', 'nearest');

disp('✅ Smooth elliptical resolution function created.');
disp('   VERY GRADUAL transition with EXTENDED ZONES:');
disp('   - 0-20km: 0.08-0.4 km (HK core)');
disp('   - 20-50km: 0.4-1 km (Inner transition)');
disp('   - 50-100km: 1-2.5 km (Outer transition)');
disp('   - 100-160km: 2.5-5 km (Approach to base)');
disp('   - 160-250km: 5-8 km (Merging with base)');
disp('   - Other areas: 3km (base)');
disp('   - Open ocean: 8-18km');

%% ========================================================================
%  STEP 5: CREATE EDGEFX WITH SMOOTH ELLIPTICAL REFINEMENT
% ========================================================================
disp('[5/9] Creating edgefx with smooth elliptical refinement...');

try
    fh = edgefx('geodata', gdat, ...
                'ef', @(x,y) F_res(x,y), ...
                'fs', fs, ...
                'wl', 20, ...
                'max_el', max_el, ...
                'g', grade, ...
                'dt', 1, ...
                'h0', coarse_base);
    disp('✅ Edgefx created with smooth elliptical refinement.');
catch
    disp('⚠️ Using standard edgefx...');
    fh = edgefx('geodata', gdat, ...
                'fs', fs, ...
                'wl', 20, ...
                'max_el', max_el, ...
                'g', grade, ...
                'dt', 1, ...
                'h0', coarse_base);
    disp('✅ Edgefx created with standard gradient.');
end

%% ========================================================================
%  STEP 6: GENERATE MESH
% ========================================================================
disp('[6/9] Generating mesh (5-15 minutes)...');
tic;

mshopts = meshgen('ef', fh, ...
                  'bou', gdat, ...
                  'plot_on', 1, ...
                  'nscreen', 14, ...
                  'cleanup', 0);

mshopts = mshopts.build;

m = mshopts.grd;
nodes = m.p;
elements = m.t;

NP = size(nodes, 1);
NE = size(elements, 1);

mesh_time = toc;
disp(['✅ Mesh generated in ', num2str(mesh_time/60), ' minutes']);
disp(['   Nodes: ', num2str(NP), ', Elements: ', num2str(NE)]);

%% ========================================================================
%  STEP 7: GENTLE MESH CLEANING
% ========================================================================
disp('[7/9] Gentle mesh cleaning...');

try
    m = clean(m, 'mqa', 0.30, 'ds', 2, 'con', 9);
    disp('   ✅ Cleaning completed');
catch
    disp('   ⚠️ Cleaning skipped');
end

nodes = m.p;
elements = m.t;
NP = size(nodes, 1);
NE = size(elements, 1);

%% ========================================================================
%  STEP 8: ASSIGN BOUNDARY CONDITIONS AND BATHYMETRY
% ========================================================================
disp('[8/9] Assigning boundary conditions and bathymetry...');

m = make_bc(m, 'auto', gdat);
m = interp(m, gebco_file);

% Extract depths
if isfield(m, 'd')
    depths = m.d;
elseif isfield(m, 'z')
    depths = m.z;
elseif isfield(m, 'b')
    depths = m.b;
else
    depths = zeros(NP, 1);
end

if isempty(depths)
    depths = m.b;
end

%% ========================================================================
%  STEP 9: WRITE FORT.14 FILE
% ========================================================================
disp('[9/9] Writing fort.14 file...');

output_fort14 = [output_dir, 'fort.14'];

try
    write(m, output_fort14, 'f14');
    disp(['✅ fort.14 saved to: ', output_fort14]);
catch
    write_standard_fort14(output_fort14, nodes, elements, depths);
    disp(['✅ fort.14 saved to: ', output_fort14]);
end

write_standard_fort14([output_dir, 'fort.14_standard'], nodes, elements, depths);
write_node_file([output_dir, 'nodes.xyz'], nodes, depths);

%% ========================================================================
%  STEP 10: COMPREHENSIVE VISUALIZATION
% ========================================================================
disp('[10/10] Generating visualization...');

% Calculate element sizes
elem_sizes = zeros(NE, 1);
for i = 1:NE
    coord = nodes(elements(i,:), :);
    dist1 = sqrt((coord(1,1)-coord(2,1))^2 + (coord(1,2)-coord(2,2))^2);
    dist2 = sqrt((coord(2,1)-coord(3,1))^2 + (coord(2,2)-coord(3,2))^2);
    dist3 = sqrt((coord(3,1)-coord(1,1))^2 + (coord(3,2)-coord(1,2))^2);
    elem_sizes(i) = mean([dist1, dist2, dist3]);
end
elem_sizes_km = elem_sizes / 1000;

% Calculate elliptical distance for each element
dist_ellipse_elem = zeros(NE, 1);
theta = 45 * pi / 180;
a = 1.5; b = 0.7;
for i = 1:NE
    centroid_lon = mean(nodes(elements(i,:), 1));
    centroid_lat = mean(nodes(elements(i,:), 2));
    dx = centroid_lon - hk_lon;
    dy = centroid_lat - hk_lat;
    x_rot = dx * cos(theta) + dy * sin(theta);
    y_rot = -dx * sin(theta) + dy * cos(theta);
    dist_ellipse_elem(i) = sqrt((x_rot / a).^2 + (y_rot / b).^2) * 111;
end

% Calculate aspect ratios
aspect_ratios = calculate_aspect_ratios(nodes, elements);

% Create comprehensive figure
figure('Color', 'w', 'Name', 'Hong Kong ADCIRC Mesh - Optimized', ...
       'Position', [50, 50, 1900, 1200]);

% Panel 1: Full triangular mesh
subplot(2,3,1);
triplot(elements, nodes(:,1), nodes(:,2), 'b-', 'LineWidth', 0.08);
hold on;
plot(hk_lon, hk_lat, 'r*', 'MarkerSize', 20, 'LineWidth', 3);
title('Optimized Smooth Elliptical Mesh');
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
axis equal;
grid on;
xlim([110 120]);
ylim([18 26]);
text(110.5, 25.5, '18km', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12);
text(113.8, 22.5, '0.08-0.6km', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 12);
text(111.5, 24.5, 'Base: 3km', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 10);

% Panel 2: Zoom to Hong Kong
subplot(2,3,2);
zoom_bbox = [113.0 115.5; 21.5 23.5];
idx_nodes = nodes(:,1) >= 113.0 & nodes(:,1) <= 115.5 & ...
            nodes(:,2) >= 21.5 & nodes(:,2) <= 23.5;
idx_elements = all(ismember(elements, find(idx_nodes)), 2);
triplot(elements(idx_elements,:), nodes(:,1), nodes(:,2), 'b-', 'LineWidth', 0.5);
hold on;
plot(hk_lon, hk_lat, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
title('Zoom: Hong Kong Coast');
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
axis equal;
grid on;

% Panel 3: Resolution vs Elliptical Distance
subplot(2,3,3);
scatter(dist_ellipse_elem, elem_sizes_km, 3, 'b', 'filled', 'MarkerEdgeAlpha', 0.1);
xlabel('Elliptical Distance from HK (km)');
ylabel('Element Size (km)');
title('VERY SMOOTH Resolution Gradient');
grid on;
hold on;
xlim([0 500]);
ylim([0 20]);

% Add ideal smooth curve with multi-zone
dist_range = 0:500;
ideal_size = zeros(size(dist_range));
for i = 1:length(dist_range)
    d = dist_range(i);
    if d <= 20
        ideal_size(i) = 80 + (400-80) * (d / 20);
    elseif d <= 50
        ideal_size(i) = 400 + (1000-400) * ((d-20) / 30);
    elseif d <= 100
        ideal_size(i) = 1000 + (2500-1000) * ((d-50) / 50);
    elseif d <= 160
        ideal_size(i) = 2500 + (5000-2500) * ((d-100) / 60);
    elseif d <= 250
        ideal_size(i) = 5000 + (8000-5000) * ((d-160) / 90);
    else
        ideal_size(i) = 8000 + (18000-8000) * (1 - exp(-0.015 * (d-250)));
    end
end
plot(dist_range, ideal_size/1000, 'r--', 'LineWidth', 2);

% Add zone boundaries
xline([20, 50, 100, 160, 250], '--', 'Color', [0.5 0.5 0.5]);
text(10, 18, 'HK', 'Color', 'r', 'FontWeight', 'bold');
text(35, 18, 'Inner', 'Color', 'm', 'FontWeight', 'bold');
text(75, 18, 'Outer', 'Color', 'g', 'FontWeight', 'bold');
text(130, 18, 'Approach', 'Color', 'c', 'FontWeight', 'bold');
text(205, 18, 'Merge', 'Color', 'b', 'FontWeight', 'bold');
legend('Elements', 'Ideal Smooth Curve', 'Location', 'northwest');

% Panel 4: Element size distribution
subplot(2,3,4);
histogram(elem_sizes_km, 50, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Element Size (km)');
ylabel('Number of Elements');
title('Element Size Distribution');
grid on;
xlim([0 20]);
hold on;
xline([0.08, 0.6], 'r--', 'LineWidth', 2);
xline([18], 'g--', 'LineWidth', 2);
legend('Elements', 'HK Target', 'Open Ocean');

% Panel 5: Spatial resolution map with stretched ellipses
subplot(2,3,5);
trisurf(elements, nodes(:,1), nodes(:,2), elem_sizes_km, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);
view(2);
hold on;
plot(hk_lon, hk_lat, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
colorbar;
colormap(jet);
caxis([0 20]);
title('Spatial Resolution Map (km)');
xlabel('Longitude (°E)');
ylabel('Latitude (°N)');
axis equal;
grid on;
xlim([110 120]);
ylim([18 26]);

% Add stretched ellipses
theta_plot = linspace(0, 2*pi, 100);
a_plot = 1.5; b_plot = 0.7;
for radius = [20, 50, 100, 160, 250]
    x_ellipse = radius/111 * a_plot * cos(theta_plot);
    y_ellipse = radius/111 * b_plot * sin(theta_plot);
    x_rot_ellipse = x_ellipse * cos(theta) - y_ellipse * sin(theta);
    y_rot_ellipse = x_ellipse * sin(theta) + y_ellipse * cos(theta);
    x_ellipse_final = hk_lon + x_rot_ellipse;
    y_ellipse_final = hk_lat + y_rot_ellipse;
    plot(x_ellipse_final, y_ellipse_final, 'w--', 'LineWidth', 1.5);
end

% Label the ellipses
text(114.5, 22.8, '20km', 'Color', 'r', 'FontSize', 10);
text(115.0, 23.3, '50km', 'Color', 'm', 'FontSize', 10);
text(115.5, 24.0, '100km', 'Color', 'g', 'FontSize', 10);
text(116.0, 24.7, '160km', 'Color', 'c', 'FontSize', 10);
text(116.8, 25.5, '250km', 'Color', 'b', 'FontSize', 10);

% Panel 6: Triangle quality
subplot(2,3,6);
histogram(aspect_ratios, 50, 'FaceColor', 'g', 'EdgeColor', 'k');
xlabel('Aspect Ratio');
ylabel('Number of Elements');
title('Triangle Quality');
grid on;
xlim([1 5]);
xline(2, 'r--', 'LineWidth', 2);
legend('Elements', 'Good Quality');

% Save figure
saveas(gcf, [output_dir, 'hong_kong_mesh_optimized.png']);
disp(['✅ Visualization saved to: ', output_dir, 'hong_kong_mesh_optimized.png']);

%% ========================================================================
%  STATISTICS SUMMARY
% ========================================================================
disp(' ');
disp('================================================================');
disp('📊 MESH STATISTICS SUMMARY');
disp('================================================================');
fprintf('Total Nodes: %d\n', NP);
fprintf('Total Elements: %d\n', NE);
fprintf('Average Element Size: %.2f km\n', mean(elem_sizes_km));
fprintf('Minimum Element Size: %.2f km (%.0f m)\n', min(elem_sizes_km), min(elem_sizes));
fprintf('Maximum Element Size: %.2f km (%.0f m)\n', max(elem_sizes_km), max(elem_sizes));
fprintf('Median Element Size: %.2f km\n', median(elem_sizes_km));

% Resolution zones
zone1 = sum(elem_sizes_km >= 0.08 & elem_sizes_km <= 0.4);
zone2 = sum(elem_sizes_km > 0.4 & elem_sizes_km <= 1);
zone3 = sum(elem_sizes_km > 1 & elem_sizes_km <= 2.5);
zone4 = sum(elem_sizes_km > 2.5 & elem_sizes_km <= 5);
zone5 = sum(elem_sizes_km > 5 & elem_sizes_km <= 8);
zone6 = sum(elem_sizes_km > 8 & elem_sizes_km <= 18);

fprintf('\n📈 RESOLUTION ZONES:\n');
fprintf('HK core (0.08-0.4km): %d elements (%.1f%%)\n', zone1, zone1/NE*100);
fprintf('Inner (0.4-1km): %d elements (%.1f%%)\n', zone2, zone2/NE*100);
fprintf('Outer (1-2.5km): %d elements (%.1f%%)\n', zone3, zone3/NE*100);
fprintf('Approach (2.5-5km): %d elements (%.1f%%)\n', zone4, zone4/NE*100);
fprintf('Merge (5-8km): %d elements (%.1f%%)\n', zone5, zone5/NE*100);
fprintf('Open ocean (8-18km): %d elements (%.1f%%)\n', zone6, zone6/NE*100);

% Triangle quality
avg_aspect = mean(aspect_ratios);
max_aspect = max(aspect_ratios);
fprintf('\n📐 TRIANGLE QUALITY:\n');
fprintf('Average Aspect Ratio: %.2f\n', avg_aspect);
fprintf('Maximum Aspect Ratio: %.2f\n', max_aspect);

fprintf('\n📁 FILES GENERATED:\n');
fprintf('   ✅ fort.14 (Full mesh for ADCIRC)\n');
fprintf('   ✅ fort.14_standard (Standard format)\n');
fprintf('   ✅ nodes.xyz (Node coordinates)\n');
fprintf('   ✅ hong_kong_mesh_optimized.png (Visualization)\n');

disp('================================================================');
disp('✅ COMPLETE PIPELINE FINISHED!');
disp(['📁 Output directory: ', output_dir]);
disp('================================================================');

%% ========================================================================
%  SUPPORTING FUNCTIONS
% ========================================================================

function write_standard_fort14(filename, nodes, elements, depths)
    NE = size(elements, 1);
    NP = size(nodes, 1);
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    fprintf(fid, 'Hong Kong ADCIRC Mesh: Optimized smooth elliptical refinement\n');
    fprintf(fid, '%d %d\n', NE, NP);
    
    for i = 1:NP
        fprintf(fid, '%8d %12.6f %12.6f %12.3f\n', ...
            i, nodes(i,1), nodes(i,2), depths(i));
    end
    
    for i = 1:NE
        fprintf(fid, '%8d %8d %8d %8d\n', ...
            i, elements(i,1), elements(i,2), elements(i,3));
    end
    
    fclose(fid);
end

function write_node_file(filename, nodes, depths)
    fid = fopen(filename, 'w');
    if fid == -1
        return;
    end
    
    for i = 1:size(nodes,1)
        fprintf(fid, '%12.6f %12.6f %12.3f\n', ...
            nodes(i,1), nodes(i,2), depths(i));
    end
    
    fclose(fid);
end

function aspect_ratios = calculate_aspect_ratios(nodes, elements)
    NE = size(elements, 1);
    aspect_ratios = zeros(NE, 1);
    
    for i = 1:NE
        coord = nodes(elements(i,:), :);
        area = polyarea(coord(:,1), coord(:,2));
        sides = [sqrt((coord(1,1)-coord(2,1))^2 + (coord(1,2)-coord(2,2))^2);
                 sqrt((coord(2,1)-coord(3,1))^2 + (coord(2,2)-coord(3,2))^2);
                 sqrt((coord(3,1)-coord(1,1))^2 + (coord(3,2)-coord(1,2))^2)];
        if area > 0
            s = mean(sides);
            inradius = area / s;
            circumradius = (sides(1)*sides(2)*sides(3)) / (4*area);
            aspect_ratios(i) = circumradius / inradius;
        else
            aspect_ratios(i) = 1;
        end
    end
end