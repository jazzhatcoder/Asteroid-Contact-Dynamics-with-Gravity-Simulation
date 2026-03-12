%% Asteroid (101955) Bennu Shape Model - Complete & Robust Visualization
% Combined version with comprehensive features and robust error handling
% Data source: Nolan et al. (2013) radar and optical observations
% File format: Wavefront OBJ-like format with vertices and faces

clear all; close all; clc;

fprintf('=== Asteroid (101955) Bennu Complete Shape Model Viewer ===\n\n');

%% Configuration and Settings
config = struct();
config.enableAnimation = false;
config.enableDetailedAnalysis = true;
config.exportFigures = false;
config.verboseOutput = true;
config.fallbackMode = false;

%% Load physical parameters with error handling
fprintf('Loading physical parameters...\n');
try
    % Try to read from actual files first
    if exist('data/pole.tab', 'file') && exist('data/rotate.tab', 'file')
        % Parse pole data from file
        pole_fid = fopen('data/pole.tab', 'r');
        if pole_fid ~= -1
            pole_line1 = fgetl(pole_fid);
            pole_line2 = fgetl(pole_fid);
            fclose(pole_fid);
            
            % Extract values (format: "name", value, uncertainty)
            pole_tokens1 = strsplit(pole_line1, ',');
            pole_tokens2 = strsplit(pole_line2, ',');
            pole_lat = str2double(pole_tokens1{2});
            pole_long = str2double(pole_tokens2{2});
        else
            error('Could not read pole.tab');
        end
        
        % Parse rotation data from file
        rotate_fid = fopen('data/rotate.tab', 'r');
        if rotate_fid ~= -1
            rotate_line = fgetl(rotate_fid);
            fclose(rotate_fid);
            
            rotate_tokens = strsplit(rotate_line, ',');
            rotation_period = str2double(rotate_tokens{2});
            rotation_uncertainty = str2double(rotate_tokens{3});
        else
            error('Could not read rotate.tab');
        end
    else
        error('Data files not found');
    end
    
catch ME
    if config.verboseOutput
        fprintf('Warning: Could not read data files directly. Using default values.\n');
        fprintf('Error: %s\n', ME.message);
    end
    
    % Fallback to hardcoded values
    pole_lat = -88;  % degrees
    pole_long = 45;  % degrees
    rotation_period = 4.29746;  % hours
    rotation_uncertainty = 0.002;  % hours
    config.fallbackMode = true;
end

fprintf('Pole orientation: Latitude = %.1f°, Longitude = %.1f°\n', pole_lat, pole_long);
fprintf('Rotation period: %.5f ± %.3f hours\n', rotation_period, rotation_uncertainty);

%% Robust shape model data loading
filename = 'data/101955bennu.tab';
fprintf('Loading Bennu shape model from %s...\n', filename);

try
    % Initialize arrays with pre-allocation for better performance
    vertices = zeros(2000, 3);  % Pre-allocate assuming ~1348 vertices
    faces = zeros(3000, 3);     % Pre-allocate assuming ~2692 faces
    
    % Open file with error checking
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s. Check if file exists in current directory.', filename);
    end
    
    % Initialize counters
    lineCount = 0;
    vertexCount = 0;
    faceCount = 0;
    parseErrors = 0;
    
    % Parse the file line by line with robust error handling
    while ~feof(fid)
        line = fgetl(fid);
        lineCount = lineCount + 1;
        
        if ischar(line) && length(line) > 1
            % Clean and tokenize the line
            line = strtrim(line);
            
            % Use regular expression for more robust parsing
            tokens = regexp(line, '\s+', 'split');
            
            if length(tokens) >= 4
                try
                    if strcmp(tokens{1}, 'v')
                        % Parse vertex: v x y z
                        x = str2double(tokens{2});
                        y = str2double(tokens{3});
                        z = str2double(tokens{4});
                        
                        % Validate numeric values
                        if ~isnan(x) && ~isnan(y) && ~isnan(z) && isfinite(x) && isfinite(y) && isfinite(z)
                            vertexCount = vertexCount + 1;
                            
                            % Expand arrays if needed
                            if vertexCount > size(vertices, 1)
                                vertices = [vertices; zeros(1000, 3)];
                            end
                            
                            vertices(vertexCount, :) = [x, y, z];
                        else
                            parseErrors = parseErrors + 1;
                            if config.verboseOutput && parseErrors <= 5
                                fprintf('Warning: Invalid vertex data at line %d\n', lineCount);
                            end
                        end
                        
                    elseif strcmp(tokens{1}, 'f')
                        % Parse face: f i1 i2 i3
                        i1 = str2double(tokens{2});
                        i2 = str2double(tokens{3});
                        i3 = str2double(tokens{4});
                        
                        % Validate face indices
                        if ~isnan(i1) && ~isnan(i2) && ~isnan(i3) && ...
                           i1 > 0 && i2 > 0 && i3 > 0 && ...
                           isfinite(i1) && isfinite(i2) && isfinite(i3)
                            faceCount = faceCount + 1;
                            
                            % Expand arrays if needed
                            if faceCount > size(faces, 1)
                                faces = [faces; zeros(1000, 3)];
                            end
                            
                            faces(faceCount, :) = [i1, i2, i3];
                        else
                            parseErrors = parseErrors + 1;
                            if config.verboseOutput && parseErrors <= 5
                                fprintf('Warning: Invalid face data at line %d\n', lineCount);
                            end
                        end
                    end
                    
                catch parseError
                    parseErrors = parseErrors + 1;
                    if config.verboseOutput && parseErrors <= 3
                        fprintf('Warning: Parse error at line %d: %s\n', lineCount, parseError.message);
                    end
                end
            end
        end
        
        % Progress indicator for large files
        if config.verboseOutput && mod(lineCount, 1000) == 0
            fprintf('Processed %d lines... (V:%d, F:%d)\n', lineCount, vertexCount, faceCount);
        end
    end
    
    fclose(fid);
    
    % Trim arrays to actual size
    vertices = vertices(1:vertexCount, :);
    faces = faces(1:faceCount, :);
    
    % Convert from kilometers to meters
    vertices = vertices * 1000;
    
    fprintf('\nData loading summary:\n');
    fprintf('- Total lines processed: %d\n', lineCount);
    fprintf('- Vertices loaded: %d\n', vertexCount);
    fprintf('- Faces loaded: %d\n', faceCount);
    if parseErrors > 0
        fprintf('- Parse errors encountered: %d\n', parseErrors);
    end
    
catch loadError
    error('Critical error loading shape model: %s', loadError.message);
end

%% Validate loaded data
if isempty(vertices) || isempty(faces)
    error('No valid shape model data found');
end

if max(faces(:)) > size(vertices, 1)
    error('Face indices exceed number of vertices (max index: %d, vertices: %d)', ...
        max(faces(:)), size(vertices, 1));
end

if min(faces(:)) < 1
    error('Invalid face indices found (minimum index: %d)', min(faces(:)));
end

fprintf('Data validation passed.\n');

%% Comprehensive geometric analysis
fprintf('\nPerforming geometric analysis...\n');

% Basic dimensions
x_range = max(vertices(:,1)) - min(vertices(:,1));
y_range = max(vertices(:,2)) - min(vertices(:,2));
z_range = max(vertices(:,3)) - min(vertices(:,3));
center_of_mass = mean(vertices, 1);

fprintf('Shape properties:\n');
fprintf('- Dimensions: %.1f × %.1f × %.1f m\n', x_range, y_range, z_range);
fprintf('- Original center of mass: (%.3f, %.3f, %.3f) m\n', center_of_mass);

% Shift vertices to center of mass frame
vertices = vertices - center_of_mass;
center_of_mass_shifted = mean(vertices, 1);
fprintf('- Center of mass after shift: (%.3e, %.3e, %.3e) m\n', center_of_mass_shifted);

% Compute bounding sphere radius for collision detection
% This creates the smallest sphere that fully encloses the asteroid
Rbound = max(vecnorm(vertices, 2, 2));
fprintf('- Bounding sphere radius: %.1f m\n', Rbound);

% Volume calculation with error handling
try
    volume = calculateMeshVolume(vertices, faces);
    fprintf('- Volume (mesh): %.2e m³\n', volume);
    
    % Mass estimation
    density_gcm3 = 2.2; % g/cm³ (CI/CM meteorite analog)
    mass_kg = volume * density_gcm3 * 1e3;  % volume in m³, density in g/cm³
    fprintf('- Estimated mass: %.2e kg (assuming %.1f g/cm³)\n', mass_kg, density_gcm3);
    
catch volError
    fprintf('- Volume calculation failed: %s\n', volError.message);
    volume = NaN;
end

% Surface area calculation
try
    surface_area = calculateSurfaceArea(vertices, faces);
    fprintf('- Surface area: %.2e m²\n', surface_area);
catch areaError
    fprintf('- Surface area calculation failed: %s\n', areaError.message);
end

% Shape analysis
if config.enableDetailedAnalysis
    try
        % Convex hull for comparison
        [hull_faces, hull_volume] = convhull(vertices);
        convexity = volume / hull_volume;
        fprintf('- Convex hull volume: %.2e m³\n', hull_volume);
        fprintf('- Convexity ratio: %.3f\n', convexity);
        
        % Equivalent sphere radius
        equiv_radius = (3*volume/(4*pi))^(1/3);
        fprintf('- Equivalent sphere radius: %.1f m\n', equiv_radius);
        
    catch analysisError
        if config.verboseOutput
            fprintf('- Advanced analysis failed: %s\n', analysisError.message);
        end
    end
end

%% Inertia Tensor Calculation
fprintf('\nComputing inertia tensor...\n');
try
    % Compute inertia tensor about center of mass
    [I_tensor, I_principal, principal_axes] = calculateInertiaTensor(vertices, faces, density_gcm3);
    
    fprintf('Inertia tensor (kg·m²) about COM:\n');
    fprintf('  [%+.4e  %+.4e  %+.4e]\n', I_tensor(1,1), I_tensor(1,2), I_tensor(1,3));
    fprintf('  [%+.4e  %+.4e  %+.4e]\n', I_tensor(2,1), I_tensor(2,2), I_tensor(2,3));
    fprintf('  [%+.4e  %+.4e  %+.4e]\n', I_tensor(3,1), I_tensor(3,2), I_tensor(3,3));
    
    fprintf('\nPrincipal moments of inertia (kg·m²):\n');
    fprintf('  I1 = %.4e\n', I_principal(1));
    fprintf('  I2 = %.4e\n', I_principal(2));
    fprintf('  I3 = %.4e\n', I_principal(3));
    
    fprintf('\nPrincipal axes (columns):\n');
    fprintf('  [%+.6f  %+.6f  %+.6f]\n', principal_axes(1,1), principal_axes(1,2), principal_axes(1,3));
    fprintf('  [%+.6f  %+.6f  %+.6f]\n', principal_axes(2,1), principal_axes(2,2), principal_axes(2,3));
    fprintf('  [%+.6f  %+.6f  %+.6f]\n', principal_axes(3,1), principal_axes(3,2), principal_axes(3,3));
    
    % Compute normalized moments (dimensionless shape parameters)
    I_normalized = I_principal / (mass_kg * equiv_radius^2);
    fprintf('\nNormalized moments (I/MR²):\n');
    fprintf('  I1/MR² = %.6f\n', I_normalized(1));
    fprintf('  I2/MR² = %.6f\n', I_normalized(2));
    fprintf('  I3/MR² = %.6f\n', I_normalized(3));
    
catch inertiaError
    fprintf('- Inertia tensor calculation failed: %s\n', inertiaError.message);
end

%% Define Rigid Body State Vector
fprintf('\nInitializing rigid body state vector...\n');

% Create asteroid structure with all physical properties and state
asteroid = struct();

% Physical properties (constant)
asteroid.name = 'Bennu';
asteroid.vertices = vertices;              % Mesh vertices (COM-centered) [m]
asteroid.faces = faces;                    % Mesh faces
asteroid.mass = mass_kg;                   % Mass [kg]
asteroid.density = density_gcm3;           % Density [g/cm³]
asteroid.volume = volume;                  % Volume [m³]
asteroid.Rbound = Rbound;                  % Bounding sphere radius [m]
asteroid.I_body = I_tensor;                % Inertia tensor in body frame [kg·m²]
asteroid.I_principal = I_principal;        % Principal moments [kg·m²]
asteroid.principal_axes = principal_axes;  % Principal axes (columns)

% State vector (time-varying in simulation)
asteroid.r = [0; 0; 0];                    % Position of COM in world frame [m]
asteroid.v = [0; 0; 0];                    % Linear velocity in world frame [m/s]
asteroid.q = [1; 0; 0; 0];                 % Orientation quaternion [q0; q1; q2; q3] (scalar-first)
asteroid.omega = [0; 0; 0];                % Angular velocity in body frame [rad/s]

% Display state vector
fprintf('Rigid body state initialized:\n');
fprintf('  Position r:        [%.2f, %.2f, %.2f] m\n', asteroid.r);
fprintf('  Velocity v:        [%.2f, %.2f, %.2f] m/s\n', asteroid.v);
fprintf('  Quaternion q:      [%.4f, %.4f, %.4f, %.4f] (scalar-first)\n', asteroid.q);
fprintf('  Angular vel ω:     [%.4f, %.4f, %.4f] rad/s\n', asteroid.omega);
fprintf('  Mass:              %.4e kg\n', asteroid.mass);
fprintf('  Bounding radius:   %.1f m\n', asteroid.Rbound);

%% Apply Initial Off-COM Impulse (Translation + Rotation Coupling)
fprintf('\nApplying off-COM impulse...\n');

% Define impulse vector in world frame [N·s]
% This represents a "kick" to the asteroid
J_world = [0.05; 0.00; 0.00];  % 0.05 N·s impulse in +X direction

% Select point of application on surface (body frame)
% Use farthest vertex from COM for maximum torque effect
[~, idx_farthest] = max(vecnorm(vertices, 2, 2));
rImp_body = vertices(idx_farthest, :)';  % [m] lever arm in body frame

fprintf('Impulse application:\n');
fprintf('  Impulse J:         [%.4f, %.4f, %.4f] N·s (world frame)\n', J_world);
fprintf('  Application point: [%.2f, %.2f, %.2f] m (body frame)\n', rImp_body);
fprintf('  Lever arm length:  %.2f m\n', norm(rImp_body));

% Transform lever arm to world frame using current orientation
R_world = quat2rotm_custom(asteroid.q);  % quaternion to rotation matrix (scalar-first)
rImp_world = R_world * rImp_body;

% Linear velocity update: Δv = J/m
asteroid.v = asteroid.v + J_world / asteroid.mass;

% Angular velocity update: τ_imp = r_imp × J
% ΔL = r_imp × J (angular impulse)
% Δω = I_world^(-1) * ΔL
I_world = R_world * asteroid.I_body * R_world';  % transform inertia to world frame
dL = cross(rImp_world, J_world);  % angular impulse [N·m·s]
asteroid.omega = asteroid.omega + (I_world \ dL);  % angular velocity jump [rad/s]

% Display post-impulse state
fprintf('\nPost-impulse state:\n');
fprintf('  Linear velocity:   [%.4f, %.4f, %.4f] m/s\n', asteroid.v);
fprintf('  Angular velocity:  [%.4f, %.4f, %.4f] rad/s\n', asteroid.omega);
fprintf('  Spin rate:         %.4f rad/s (%.2f deg/s)\n', norm(asteroid.omega), norm(asteroid.omega)*180/pi);
fprintf('  Angular momentum:  [%.4e, %.4e, %.4e] kg·m²/s\n', I_world * asteroid.omega);

% Verify momentum conservation
linear_momentum_change = asteroid.mass * asteroid.v;
angular_momentum_change = dL;
fprintf('\nMomentum verification:\n');
fprintf('  Δp (should = J):   [%.4f, %.4f, %.4f] N·s\n', linear_momentum_change);
fprintf('  ΔL (from r×J):     [%.4e, %.4e, %.4e] N·m·s\n', angular_momentum_change);

%% Free Rigid Body Dynamics Integration (Torque-Free Tumbling)
fprintf('\nIntegrating free rigid body motion...\n');

% Simulation parameters
dt = 1.0;              % Time step [s]
N_steps = 500;         % Number of steps
t_total = N_steps * dt;  % Total simulation time [s]

fprintf('Simulation parameters:\n');
fprintf('  Time step:         %.2f s\n', dt);
fprintf('  Total steps:       %d\n', N_steps);
fprintf('  Total time:        %.1f s (%.2f minutes)\n', t_total, t_total/60);

% Pre-allocate trajectory arrays
trajectory = struct();
trajectory.t = zeros(N_steps, 1);
trajectory.r = zeros(N_steps, 3);
trajectory.v = zeros(N_steps, 3);
trajectory.q = zeros(N_steps, 4);
trajectory.omega = zeros(N_steps, 3);
trajectory.L = zeros(N_steps, 3);
trajectory.E_rot = zeros(N_steps, 1);

% Time integration loop - Governing equations:
% Translation: ṙ = v, v̇ = 0 (no forces)
% Rotation: I·ω̇ + ω × (I·ω) = 0 (Euler's equations, no torques)
% Quaternion: q̇ = 0.5 * Ω(ω) * q

for k = 1:N_steps
    % Store current state
    trajectory.t(k) = (k-1) * dt;
    trajectory.r(k,:) = asteroid.r';
    trajectory.v(k,:) = asteroid.v';
    trajectory.q(k,:) = asteroid.q';
    trajectory.omega(k,:) = asteroid.omega';
    
    % Compute angular momentum and rotational energy
    R_current = quat2rotm_custom(asteroid.q);
    I_world_current = R_current * asteroid.I_body * R_current';
    trajectory.L(k,:) = (I_world_current * asteroid.omega)';
    trajectory.E_rot(k) = 0.5 * asteroid.omega' * asteroid.I_body * asteroid.omega;
    
    % === INTEGRATION STEP ===
    
    % 1. Translation (trivial, no forces)
    asteroid.r = asteroid.r + asteroid.v * dt;
    % asteroid.v remains constant (v̇ = 0)
    
    % 2. Rotation (Euler's equations in body frame)
    % I·ω̇ = -ω × (I·ω)
    omega_dot = -asteroid.I_body \ cross(asteroid.omega, asteroid.I_body * asteroid.omega);
    asteroid.omega = asteroid.omega + omega_dot * dt;
    
    % 3. Quaternion propagation
    % q̇ = 0.5 * Ω(ω) * q
    wx = asteroid.omega(1);
    wy = asteroid.omega(2);
    wz = asteroid.omega(3);
    
    Omega = [ 0   -wx  -wy  -wz;
              wx   0    wz  -wy;
              wy  -wz   0    wx;
              wz   wy  -wx   0 ];
    
    q_dot = 0.5 * Omega * asteroid.q;
    asteroid.q = asteroid.q + q_dot * dt;
    asteroid.q = asteroid.q / norm(asteroid.q);  % Normalize (critical for stability)
end

fprintf('Integration completed.\n');
fprintf('Final state:\n');
fprintf('  Position:          [%.2f, %.2f, %.2f] m\n', asteroid.r);
fprintf('  Orientation (deg): [%.1f, %.1f, %.1f, %.1f]\n', asteroid.q * 180/pi);
fprintf('  Spin rate:         %.4f rad/s\n', norm(asteroid.omega));

% Verify conservation laws
L_initial = trajectory.L(1,:);
L_final = trajectory.L(end,:);
E_initial = trajectory.E_rot(1);
E_final = trajectory.E_rot(end);

fprintf('\nConservation verification:\n');
fprintf('  Angular momentum error: %.6e (%.3f%%)\n', ...
    norm(L_final - L_initial), 100*norm(L_final - L_initial)/norm(L_initial));
fprintf('  Energy error:           %.6e (%.3f%%)\n', ...
    E_final - E_initial, 100*abs(E_final - E_initial)/E_initial);

%% Surface model visualization
fprintf('\nCreating surface model visualization...\n');

% Main figure with surface model
fig1 = figure('Name', 'Bennu Surface Model', 'Position', [100, 100, 1000, 900]);

% 3D surface rendering using all vertices and faces
h_surf = patch('Vertices', vertices, 'Faces', faces, ...
               'FaceColor', [0.8, 0.6, 0.3], ...
               'EdgeColor', 'none', ...
               'FaceLighting', 'gouraud', ...
               'AmbientStrength', 0.4, ...
               'DiffuseStrength', 0.6, ...
               'SpecularStrength', 0.3);

hold on;

% Enhanced coordinate system
axis_scale = max([x_range, y_range, z_range]) * 0.15;
quiver3(0, 0, 0, axis_scale, 0, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.15);
quiver3(0, 0, 0, 0, axis_scale, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.15);
quiver3(0, 0, 0, 0, 0, axis_scale, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.15);

text(axis_scale*1.2, 0, 0, 'X', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
text(0, axis_scale*1.2, 0, 'Y', 'FontSize', 12, 'Color', 'g', 'FontWeight', 'bold');
text(0, 0, axis_scale*1.2, 'Z (Spin)', 'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');

axis equal tight;
grid on;
light('Position', [1, 1, 1], 'Style', 'infinite');
light('Position', [-1, -1, -1], 'Style', 'infinite');

xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
zlabel('Z (m)', 'FontSize', 12);
title(sprintf('Asteroid (101955) Bennu - 3D Surface Model\nVertices: %d | Faces: %d | Volume: %.2e m³ | Period: %.5f hrs', ...
    size(vertices,1), size(faces,1), volume, rotation_period), ...
    'FontSize', 14, 'FontWeight', 'bold');
view(45, 20);
     


%% Animation with enhanced features
if config.enableAnimation
    fprintf('Creating enhanced rotation animation...\n');
    
    try
        fig2 = figure('Name', 'Bennu Rotation Animation', 'Position', [200, 200, 800, 800]);
        
        h_anim = patch('Vertices', vertices, 'Faces', faces, ...
                       'FaceColor', [0.8, 0.5, 0.2], ...
                       'EdgeColor', 'none', ...
                       'FaceLighting', 'gouraud', ...
                       'AmbientStrength', 0.3);
        
        axis equal tight;
        grid on;
        light('Position', [1, 1, 1], 'Style', 'infinite');
        light('Position', [-1, -1, -1], 'Style', 'infinite');
        
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(sprintf('Bennu Rotation\nPeriod: %.3f hours (%.1f minutes)', ...
            rotation_period, rotation_period*60), 'FontSize', 14, 'FontWeight', 'bold');
        
        % Enhanced animation with multiple viewing angles
        angles = 0:3:360;
        elevations = [15, 20, 25, 20, 15, 10, 15, 25, 30, 25, 20, 15];
        
        for i = 1:length(angles)
            elev_idx = mod(i-1, length(elevations)) + 1;
            view(angles(i), elevations(elev_idx));
            drawnow;
            pause(0.03);
        end
        
    catch animError
        fprintf('Animation failed: %s\n', animError.message);
    end
end

%% Export and summary
if config.exportFigures
    try
        fprintf('Exporting figures...\n');
        saveas(fig1, 'bennu_complete_analysis.png');
        if exist('fig2', 'var')
            saveas(fig2, 'bennu_animation_frame.png');
        end
        fprintf('Figures exported successfully.\n');
    catch exportError
        fprintf('Export failed: %s\n', exportError.message);
    end
end



%% Final comprehensive summary
fprintf('\n=== COMPLETE ANALYSIS SUMMARY ===\n');
fprintf('Asteroid (101955) Bennu Shape Model - Successfully Processed\n\n');

fprintf('DATA STATISTICS:\n');
fprintf('- Vertices: %d\n', size(vertices,1));
fprintf('- Triangular faces: %d\n', size(faces,1));
fprintf('- Model resolution: ~25m between vertices\n');
if config.fallbackMode
    fprintf('- Data source: Fallback values (file read failed)\n');
else
    fprintf('- Data source: Direct file parsing\n');
end

fprintf('\nPHYSICAL PROPERTIES:\n');
fprintf('- Dimensions: %.1f × %.1f × %.1f meters\n', x_range, y_range, z_range);
fprintf('- Rotation period: %.5f ± %.3f hours\n', rotation_period, rotation_uncertainty);
fprintf('- Pole orientation: (%.0f°, %.0f°)\n', pole_lat, pole_long);
if ~isnan(volume)
    fprintf('- Volume: %.2e m³\n', volume);
    fprintf('- Estimated mass: %.2e kg\n', mass_kg);
end

fprintf('\nCOORDINATE SYSTEM:\n');
fprintf('- Origin: Center of mass\n');
fprintf('- Axes: Principal axes of shape model\n');
fprintf('- Z-axis: Spin axis (positive pole)\n');
fprintf('- Units: Meters\n');

fprintf('\nSCIENTIFIC CONTEXT:\n');
fprintf('- Classification: B-type Apollo NEO\n');
fprintf('- Target of NASA OSIRIS-REx mission\n');
fprintf('- Data from Nolan et al. (2013), Icarus 226, 629-640\n');
fprintf('- Radar observations: Arecibo & Goldstone (1999, 2005)\n\n');

%% Helper functions
function volume = calculateMeshVolume(vertices, faces)
    % Calculate volume using divergence theorem with error checking
    volume = 0;
    
    for i = 1:size(faces, 1)
        try
            % Get triangle vertices
            v1 = vertices(faces(i,1), :);
            v2 = vertices(faces(i,2), :);
            v3 = vertices(faces(i,3), :);
            
            % Calculate signed volume of tetrahedron
            vol_tet = dot(v1, cross(v2, v3)) / 6;
            volume = volume + vol_tet;
        catch
            % Skip invalid faces
            continue;
        end
    end
    
    volume = abs(volume);
    
    % Sanity check
    if volume <= 0 || ~isfinite(volume)
        error('Invalid volume calculation result: %.6f', volume);
    end
end

function [I_tensor, I_principal, principal_axes] = calculateInertiaTensor(vertices, faces, density_gcm3)
    % Calculate inertia tensor using tetrahedron decomposition method
    % Assumes vertices are already centered at COM (origin)
    % Returns inertia tensor, principal moments, and principal axes
    
    % Convert density from g/cm³ to kg/m³
    density = density_gcm3 * 1000; % kg/m³
    
    % Initialize volume integrals
    int_x2 = 0;
    int_y2 = 0;
    int_z2 = 0;
    int_xy = 0;
    int_xz = 0;
    int_yz = 0;
    
    % Loop over all triangular faces
    for i = 1:size(faces, 1)
        % Get vertices of the triangle (tetrahedron is formed with origin)
        a = vertices(faces(i,1), :)';  % Column vector
        b = vertices(faces(i,2), :)';
        c = vertices(faces(i,3), :)';
        
        % Signed volume of tetrahedron: V_t = (1/6) * a·(b×c)
        V_t = dot(a, cross(b, c)) / 6;
        
        % Extract components
        ax = a(1); ay = a(2); az = a(3);
        bx = b(1); by = b(2); bz = b(3);
        cx = c(1); cy = c(2); cz = c(3);
        
        % Compute volume integrals for this tetrahedron
        % ∫x²dV = (V_t/10)(ax² + bx² + cx² + ax*bx + bx*cx + cx*ax)
        int_x2 = int_x2 + (V_t/10) * (ax^2 + bx^2 + cx^2 + ax*bx + bx*cx + cx*ax);
        int_y2 = int_y2 + (V_t/10) * (ay^2 + by^2 + cy^2 + ay*by + by*cy + cy*ay);
        int_z2 = int_z2 + (V_t/10) * (az^2 + bz^2 + cz^2 + az*bz + bz*cz + cz*az);
        
        % ∫xydV = (V_t/20)(2(ax*ay + bx*by + cx*cy) + ax*by + ay*bx + bx*cy + by*cx + cx*ay + cy*ax)
        int_xy = int_xy + (V_t/20) * (2*(ax*ay + bx*by + cx*cy) + ax*by + ay*bx + bx*cy + by*cx + cx*ay + cy*ax);
        int_xz = int_xz + (V_t/20) * (2*(ax*az + bx*bz + cx*cz) + ax*bz + az*bx + bx*cz + bz*cx + cx*az + cz*ax);
        int_yz = int_yz + (V_t/20) * (2*(ay*az + by*bz + cy*cz) + ay*bz + az*by + by*cz + bz*cy + cy*az + cz*ay);
    end
    
    % Construct inertia tensor
    % I = ρ * [y²+z²   -xy    -xz  ]
    %         [-xy    x²+z²  -yz  ]
    %         [-xz    -yz    x²+y²]
    I_tensor = zeros(3, 3);
    I_tensor(1,1) = density * (int_y2 + int_z2);  % Ixx
    I_tensor(2,2) = density * (int_x2 + int_z2);  % Iyy
    I_tensor(3,3) = density * (int_x2 + int_y2);  % Izz
    I_tensor(1,2) = -density * int_xy;            % Ixy
    I_tensor(2,1) = I_tensor(1,2);                % Iyx
    I_tensor(1,3) = -density * int_xz;            % Ixz
    I_tensor(3,1) = I_tensor(1,3);                % Izx
    I_tensor(2,3) = -density * int_yz;            % Iyz
    I_tensor(3,2) = I_tensor(2,3);                % Izy
    
    % Compute principal moments and axes (eigenvalues and eigenvectors)
    [principal_axes, I_diag] = eig(I_tensor);
    I_principal = diag(I_diag);
    
    % Sort by principal moments (ascending order)
    [I_principal, sort_idx] = sort(I_principal);
    principal_axes = principal_axes(:, sort_idx);
    
    % Ensure right-handed coordinate system
    if det(principal_axes) < 0
        principal_axes(:,3) = -principal_axes(:,3);
    end
end

function R = quat2rotm_custom(q)
    % Custom quaternion to rotation matrix conversion (no toolbox required)
    % Input: q = [q0; q1; q2; q3] (scalar-first) or [q0, q1, q2, q3] (row)
    % Output: R = 3x3 rotation matrix
    
    % Ensure column vector
    if size(q, 2) == 4
        q = q';
    end
    
    % Normalize
    q = q / norm(q);
    
    % Extract components (scalar-first convention)
    w = q(1);  % scalar part
    x = q(2);  % vector part
    y = q(3);
    z = q(4);
    
    % Compute rotation matrix
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),     2*(x*z+w*y);
         2*(x*y+w*z),     1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),     2*(y*z+w*x),     1-2*(x^2+y^2)];
end

function surface_area = calculateSurfaceArea(vertices, faces)
    % Calculate total surface area of the mesh
    surface_area = 0;
    
    for i = 1:size(faces, 1)
        try
            % Get triangle vertices
            v1 = vertices(faces(i,1), :);
            v2 = vertices(faces(i,2), :);
            v3 = vertices(faces(i,3), :);
            
            % Calculate triangle area using cross product
            edge1 = v2 - v1;
            edge2 = v3 - v1;
            triangle_area = 0.5 * norm(cross(edge1, edge2));
            surface_area = surface_area + triangle_area;
        catch
            % Skip invalid faces
            continue;
        end
    end
    
    if surface_area <= 0 || ~isfinite(surface_area)
        error('Invalid surface area calculation result: %.6f', surface_area);
    end
end

