%% Asteroid (25143) Itokawa Shape Model - Complete & Robust Visualization
% Adapted from Bennu visualization code
% Data source: Gaskell Itokawa Shape Model (Hayabusa AMICA mission)
% Based on 775 Hayabusa AMICA images (Sept 11 - Nov 12, 2005)
% Model prepared: August 29, 2007 by Robert Gaskell
% File format: ICQ (Implicitly Connected Quadrilateral) and vertex-facet
% Reference: Gaskell et al. 2006, AIAA 2006-6660

clear all; close all; clc;

fprintf('=== Asteroid (25143) Itokawa Complete Shape Model Viewer ===\n\n');

%% Configuration and Settings
config = struct();
config.enableAnimation = false;
config.enableDetailedAnalysis = true;
config.exportFigures = false;
config.verboseOutput = true;
config.fallbackMode = false;
config.resolution = '512q';  % Options: '64q', '128q', '256q', '512q'

%% Load physical parameters
fprintf('Loading physical parameters...\n');

% Itokawa physical parameters (from Gaskell shape model)
% Pole definition from shape model:
% BODY2025143_POLE_RA  = 90.02564 degrees
% BODY2025143_POLE_DEC = -67.02704 degrees
% BODY2025143_PM = (129.73000, 712.14376110, 0)
% PM rate of 712.14376110 deg/day gives rotation period

pole_ra = 90.02564;  % degrees (Right Ascension)
pole_dec = -67.02704;  % degrees (Declination)
pm_rate = 712.14376110;  % degrees per day
rotation_period = 360.0 / pm_rate * 24;  % hours (360 deg / PM_rate * 24 hr/day)
rotation_uncertainty = 0.001;  % hours

% Calculate pole latitude and longitude for display
pole_lat = pole_dec;  % Declination is latitude
pole_long = pole_ra;  % Right Ascension is longitude

fprintf('Pole orientation (J2000): RA = %.5f°, Dec = %.5f°\n', pole_ra, pole_dec);
fprintf('Rotation period: %.5f hours (%.3f minutes)\n', rotation_period, rotation_period*60);

%% Robust shape model data loading
vertex_filename = sprintf('vertex/ver%s.tab', config.resolution);
fprintf('Loading Itokawa shape model from %s...\n', vertex_filename);

try
    % Open file with error checking
    fid = fopen(vertex_filename, 'r');
    if fid == -1
        error('Cannot open file: %s. Check if file exists in current directory.', vertex_filename);
    end
    
    % Read the header line to get counts
    header_line = fgetl(fid);
    header_values = sscanf(header_line, '%d');
    
    if length(header_values) >= 2
        num_vertices = header_values(1);
        num_faces = header_values(2);
        fprintf('Expected vertices: %d, faces: %d\n', num_vertices, num_faces);
    else
        error('Invalid header format in vertex file');
    end
    
    % Pre-allocate arrays
    vertices = zeros(num_vertices, 3);
    faces = zeros(num_faces, 3);
    
    % Initialize counters
    lineCount = 1;  % Already read header
    vertexCount = 0;
    faceCount = 0;
    parseErrors = 0;
    
    % Read vertices section
    fprintf('Reading vertices...\n');
    while vertexCount < num_vertices && ~feof(fid)
        line = fgetl(fid);
        lineCount = lineCount + 1;
        
        if ischar(line) && ~isempty(line)
            % Parse vertex line: index x y z
            values = sscanf(line, '%d %f %f %f');
            
            if length(values) == 4
                vertexCount = vertexCount + 1;
                vertices(vertexCount, :) = values(2:4)';
            else
                parseErrors = parseErrors + 1;
                if config.verboseOutput && parseErrors <= 5
                    fprintf('Warning: Could not parse vertex line %d: %s\n', lineCount, line);
                end
            end
        end
        
        % Progress indicator
        if config.verboseOutput && mod(vertexCount, 5000) == 0
            fprintf('  Loaded %d/%d vertices...\n', vertexCount, num_vertices);
        end
    end
    
    % Read faces section
    fprintf('Reading faces...\n');
    while faceCount < num_faces && ~feof(fid)
        line = fgetl(fid);
        lineCount = lineCount + 1;
        
        if ischar(line) && ~isempty(line)
            % Parse face line: index v1 v2 v3
            values = sscanf(line, '%d %d %d %d');
            
            if length(values) == 4
                faceCount = faceCount + 1;
                faces(faceCount, :) = values(2:4)';
            else
                parseErrors = parseErrors + 1;
                if config.verboseOutput && parseErrors <= 5
                    fprintf('Warning: Could not parse face line %d: %s\n', lineCount, line);
                end
            end
        end
        
        % Progress indicator
        if config.verboseOutput && mod(faceCount, 10000) == 0
            fprintf('  Loaded %d/%d faces...\n', faceCount, num_faces);
        end
    end
    
    fclose(fid);
    
    % Trim arrays to actual size (in case we got fewer than expected)
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
    
    % Mass estimation (Itokawa is S-type asteroid)
    density_gcm3 = 1.9; % g/cm³ (measured by Hayabusa mission)
    mass_kg = volume * density_gcm3 * 1e3;  % volume in m³, density in g/cm³
    fprintf('- Estimated mass: %.2e kg (using measured %.1f g/cm³)\n', mass_kg, density_gcm3);
    
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
asteroid.name = 'Itokawa';
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
J_world = [0.03; 0.00; 0.00];  % 0.03 N·s impulse in +X direction

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
fig1 = figure('Name', 'Itokawa Surface Model', 'Position', [100, 100, 1000, 900]);

% 3D surface rendering using all vertices and faces
h_surf = patch('Vertices', vertices, 'Faces', faces, ...
               'FaceColor', [0.7, 0.7, 0.6], ...
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
title(sprintf('Asteroid (25143) Itokawa - 3D Surface Model\nVertices: %d | Faces: %d | Volume: %.2e m³ | Period: %.3f hrs', ...
    size(vertices,1), size(faces,1), volume, rotation_period), ...
    'FontSize', 14, 'FontWeight', 'bold');
view(45, 20);

%% Animation with enhanced features
if config.enableAnimation
    fprintf('Creating enhanced rotation animation...\n');
    
    try
        fig2 = figure('Name', 'Itokawa Rotation Animation', 'Position', [200, 200, 800, 800]);
        
        h_anim = patch('Vertices', vertices, 'Faces', faces, ...
                       'FaceColor', [0.7, 0.7, 0.6], ...
                       'EdgeColor', 'none', ...
                       'FaceLighting', 'gouraud');
        
        axis equal tight;
        grid on;
        light('Position', [1, 1, 1], 'Style', 'infinite');
        light('Position', [-1, -1, -1], 'Style', 'infinite');
        
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(sprintf('Itokawa Rotation\nPeriod: %.3f hours (%.1f minutes)', ...
            rotation_period, rotation_period*60));
        
        % Enhanced animation with multiple viewing angles
        angles = 0:3:360;
        elevations = [15, 20, 25, 20, 15, 10, 15, 25, 30, 25, 20, 15];
        
        for i = 1:length(angles)
            view(angles(i), elevations(mod(i-1, length(elevations))+1));
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
        saveas(fig1, 'itokawa_complete_analysis.png');
        if exist('fig2', 'var')
            saveas(fig2, 'itokawa_animation_frame.png');
        end
        fprintf('Figures exported successfully.\n');
    catch exportError
        fprintf('Export failed: %s\n', exportError.message);
    end
end


%% Final comprehensive summary
fprintf('\n=== COMPLETE ANALYSIS SUMMARY ===\n');
fprintf('Asteroid (25143) Itokawa Shape Model - Successfully Processed\n\n');

fprintf('DATA STATISTICS:\n');
fprintf('- Resolution: %s\n', config.resolution);
fprintf('- Vertices: %d\n', size(vertices,1));
fprintf('- Triangular faces: %d\n', size(faces,1));

fprintf('\nPHYSICAL PROPERTIES:\n');
fprintf('- Dimensions: %.1f × %.1f × %.1f meters\n', x_range, y_range, z_range);
fprintf('- Rotation period: %.3f hours\n', rotation_period);
fprintf('- Pole orientation: (%.0f°, %.0f°)\n', pole_lat, pole_long);
if ~isnan(volume)
    fprintf('- Volume: %.2e m³\n', volume);
    fprintf('- Mass: %.2e kg\n', mass_kg);
    fprintf('- Density: %.1f g/cm³ (measured)\n', density_gcm3);
end

fprintf('\nCOORDINATE SYSTEM:\n');
fprintf('- Frame: Right-handed Cartesian body-fixed\n');
fprintf('- Z-axis: Rotation pole (RA=%.5f°, Dec=%.5f°)\n', pole_ra, pole_dec);
fprintf('- X-axis: Zero longitude through "black rock" feature\n');
fprintf('- Reference: Black rock at lat 3.357°S, lon 0°\n');
fprintf('- Units: Meters\n');

fprintf('\nSCIENTIFIC CONTEXT:\n');
fprintf('- Classification: S-type Apollo NEO\n');
fprintf('- Target of JAXA Hayabusa mission (2005-2010)\n');
fprintf('- First asteroid sample return mission\n');
fprintf('- Notable "peanut" shape with two lobes\n');
fprintf('- Shape model: 871 control points, 90,720 measurements\n');
fprintf('- Model accuracy: 35 cm RMS (20 cm per DOF)\n');
fprintf('- Data: 775 AMICA images (Sept-Nov 2005)\n\n');

%% Helper functions
function volume = calculateMeshVolume(vertices, faces)
    % Calculate volume using divergence theorem with error checking
    volume = 0;
    
    for i = 1:size(faces, 1)
        try
            % Get vertices of the triangle
            v1 = vertices(faces(i,1), :);
            v2 = vertices(faces(i,2), :);
            v3 = vertices(faces(i,3), :);
            
            % Calculate signed volume of tetrahedron formed by triangle and origin
            volume = volume + v1(1)*(v2(2)*v3(3) - v2(3)*v3(2)) + ...
                              v2(1)*(v3(2)*v1(3) - v3(3)*v1(2)) + ...
                              v3(1)*(v1(2)*v2(3) - v1(3)*v2(2));
        catch
            % Skip invalid faces
        end
    end
    
    volume = abs(volume) / 6;
    
    % Sanity check
    if volume <= 0 || ~isfinite(volume)
        error('Invalid volume calculation result: %.6f', volume);
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
            % Get vertices of the triangle
            v1 = vertices(faces(i,1), :);
            v2 = vertices(faces(i,2), :);
            v3 = vertices(faces(i,3), :);
            
            % Calculate edge vectors
            edge1 = v2 - v1;
            edge2 = v3 - v1;
            
            % Area is half the magnitude of the cross product
            cross_prod = cross(edge1, edge2);
            area = 0.5 * norm(cross_prod);
            surface_area = surface_area + area;
        catch
            % Skip invalid faces
        end
    end
    
    if surface_area <= 0 || ~isfinite(surface_area)
        error('Invalid surface area calculation result: %.6f', surface_area);
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


