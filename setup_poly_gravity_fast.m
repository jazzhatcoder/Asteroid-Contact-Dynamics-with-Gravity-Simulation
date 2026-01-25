%% FAST Setup: Preprocess Asteroids with REASONABLE Mesh Resolutions
% This uses appropriate mesh sizes for fast polyhedral gravity computation
% Run this INSTEAD of setup_poly_gravity.m

clear all; close all; clc;

fprintf('=== FAST Asteroid Preprocessing (Optimized Mesh Sizes) ===\n\n');

%% Configuration - Use moderate mesh resolution
N_mass_points = 50;  % Mass points per asteroid (was 100-200)

fprintf('Optimization strategy:\n');
fprintf('  - Use existing bennu() and Itokawa() functions (reasonable mesh sizes)\n');
fprintf('  - Target: ~50 mass points per body\n');
fprintf('  - Expected speedup: 50-100x vs high-res meshes\n\n');

%% Load and preprocess Bennu (using moderate resolution mesh)
fprintf('========================================\n');
fprintf('BENNU PREPROCESSING\n');
fprintf('========================================\n');

% Load Bennu mesh directly
bennu_raw = loadBennuMesh();

fprintf('Loaded Bennu mesh: %d vertices, %d faces\n', ...
    size(bennu_raw.vertices, 1), size(bennu_raw.faces, 1));

% Precompute polyhedral gravity structure
fprintf('Computing polyhedral gravity structure...\n');
polyBennu = precomputePolyhedron(bennu_raw.vertices, bennu_raw.faces, bennu_raw.density * 1000);

% Generate mass distribution using rejection sampling (no tetrahedralization needed)
fprintf('Generating %d mass points via rejection sampling...\n', N_mass_points);
bennu_raw.C_body = generateMassPoints(bennu_raw.vertices, bennu_raw.faces, N_mass_points);
bennu_raw.m_tet = ones(N_mass_points, 1) * bennu_raw.mass / N_mass_points;
bennu_raw.V_tet = bennu_raw.m_tet / (bennu_raw.density * 1000);

fprintf('✓ Bennu preprocessing complete: %d mass points\n\n', N_mass_points);

%% Load and preprocess Itokawa (using LOW resolution mesh)
fprintf('========================================\n');
fprintf('ITOKAWA PREPROCESSING\n');
fprintf('========================================\n');

% Use Itokawa() function which has configurable resolution
fprintf('Loading Itokawa with LOW resolution (64q)...\n');

% Temporarily modify Itokawa config to use low-res mesh
original_dir = cd('.');
try
    % Load low-resolution Itokawa mesh directly
    itokawa_raw = loadItokawa_lowres();
    
    fprintf('Loaded Itokawa mesh: %d vertices, %d faces\n', ...
        size(itokawa_raw.vertices, 1), size(itokawa_raw.faces, 1));
    
    % Precompute polyhedral gravity structure
    fprintf('Computing polyhedral gravity structure...\n');
    polyItokawa = precomputePolyhedron(itokawa_raw.vertices, itokawa_raw.faces, itokawa_raw.density * 1000);
    
    % Generate mass distribution
    fprintf('Generating %d mass points via rejection sampling...\n', N_mass_points);
    itokawa_raw.C_body = generateMassPoints(itokawa_raw.vertices, itokawa_raw.faces, N_mass_points);
    itokawa_raw.m_tet = ones(N_mass_points, 1) * itokawa_raw.mass / N_mass_points;
    itokawa_raw.V_tet = itokawa_raw.m_tet / (itokawa_raw.density * 1000);
    
    fprintf('✓ Itokawa preprocessing complete: %d mass points\n\n', N_mass_points);
    
catch ME
    cd(original_dir);
    error('Failed to load Itokawa: %s', ME.message);
end
cd(original_dir);

%% Save preprocessed data
fprintf('========================================\n');
fprintf('SAVING PREPROCESSED DATA\n');
fprintf('========================================\n');

save('asteroid_poly_data_fast.mat', 'polyBennu', 'polyItokawa', 'bennu_raw', 'itokawa_raw');

fprintf('✓ Saved to asteroid_poly_data_fast.mat\n');
fprintf('  polyBennu:    %d vertices, %d faces, %d edges\n', ...
    size(polyBennu.V, 1), size(polyBennu.F, 1), size(polyBennu.edges, 1));
fprintf('  polyItokawa:  %d vertices, %d faces, %d edges\n', ...
    size(polyItokawa.V, 1), size(polyItokawa.F, 1), size(polyItokawa.edges, 1));
fprintf('  bennu_raw:    %d mass points\n', length(bennu_raw.m_tet));
fprintf('  itokawa_raw:  %d mass points\n', length(itokawa_raw.m_tet));

fprintf('\n========================================\n');
fprintf('FAST PREPROCESSING COMPLETE\n');
fprintf('========================================\n');
fprintf('Mesh sizes optimized for fast computation:\n');
fprintf('  Bennu:   ~%d vertices (moderate resolution)\n', size(polyBennu.V, 1));
fprintf('  Itokawa: ~%d vertices (LOW resolution - ver64q)\n', size(polyItokawa.V, 1));
fprintf('  Mass points: %d per body\n', N_mass_points);
fprintf('  Expected simulation time: 5-30 minutes (was 6+ hours)\n\n');
fprintf('Next: Run the simulation with this data file:\n');
fprintf('  >> two_body_poly_simulation_fast\n');

%% Helper functions

function data = loadBennuMesh()
    % Load Bennu from existing data file
    filename = 'data/101955bennu.tab';
    vertices = zeros(2000, 3);
    faces = zeros(3000, 3);
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open Bennu file: %s', filename);
    end
    
    vertexCount = 0;
    faceCount = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        
        tokens = strsplit(strtrim(line));
        if isempty(tokens), continue; end
        
        if strcmp(tokens{1}, 'v')
            vertexCount = vertexCount + 1;
            vertices(vertexCount, :) = [str2double(tokens{2}), ...
                                        str2double(tokens{3}), ...
                                        str2double(tokens{4})];
        elseif strcmp(tokens{1}, 'f')
            faceCount = faceCount + 1;
            faces(faceCount, :) = [str2double(tokens{2}), ...
                                   str2double(tokens{3}), ...
                                   str2double(tokens{4})];
        end
    end
    fclose(fid);
    
    vertices = vertices(1:vertexCount, :) * 1000;  % km to m
    faces = faces(1:faceCount, :);
    
    % Center at COM
    vertices = vertices - mean(vertices, 1);
    
    % Compute volume and mass
    volume = calculateMeshVolume(vertices, faces);
    density_gcm3 = 2.2;
    mass = volume * density_gcm3 * 1e3;
    
    data.vertices = vertices;
    data.faces = faces;
    data.volume = volume;
    data.density = density_gcm3;
    data.mass = mass;
    data.Rbound = max(vecnorm(vertices, 2, 2));
end

function data = loadItokawa_lowres()
    % Load LOW-resolution Itokawa mesh (ver64q.tab)
    filename = 'vertex/ver64q.tab';
    
    fprintf('  Reading file: %s\n', filename);
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open Itokawa file: %s. Ensure ver64q.tab exists.', filename);
    end
    
    % Read header
    header_line = fgetl(fid);
    header_values = sscanf(header_line, '%d');
    num_vertices = header_values(1);
    num_faces = header_values(2);
    
    fprintf('  Expected: %d vertices, %d faces\n', num_vertices, num_faces);
    
    % Pre-allocate
    vertices = zeros(num_vertices, 3);
    faces = zeros(num_faces, 3);
    
    % Read vertices
    for i = 1:num_vertices
        line = fgetl(fid);
        values = sscanf(line, '%d %f %f %f');
        vertices(i, :) = values(2:4)';
    end
    
    % Read faces
    for i = 1:num_faces
        line = fgetl(fid);
        values = sscanf(line, '%d %d %d %d');
        faces(i, :) = values(2:4)';
    end
    
    fclose(fid);
    
    % Convert to meters and center at COM
    vertices = vertices * 1000;  % km to m
    vertices = vertices - mean(vertices, 1);
    
    % Compute volume and mass
    volume = calculateMeshVolume(vertices, faces);
    density_gcm3 = 1.9;
    mass = volume * density_gcm3 * 1e3;
    
    data.vertices = vertices;
    data.faces = faces;
    data.volume = volume;
    data.density = density_gcm3;
    data.mass = mass;
    data.Rbound = max(vecnorm(vertices, 2, 2));
end

function volume = calculateMeshVolume(vertices, faces)
    volume = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i,1), :);
        v2 = vertices(faces(i,2), :);
        v3 = vertices(faces(i,3), :);
        volume = volume + det([v1; v2; v3]);
    end
    volume = abs(volume) / 6;
end

function points = generateMassPoints(V, F, N)
    % Generate N points inside mesh using rejection sampling
    bbox_min = min(V, [], 1);
    bbox_max = max(V, [], 1);
    
    points = zeros(N, 3);
    count = 0;
    max_attempts = N * 1000;
    attempts = 0;
    
    while count < N && attempts < max_attempts
        attempts = attempts + 1;
        p = bbox_min + rand(1, 3) .* (bbox_max - bbox_min);
        
        if isInsideMesh(p, V, F)
            count = count + 1;
            points(count, :) = p;
            
            if mod(count, 10) == 0
                fprintf('  Generated %d/%d points\n', count, N);
            end
        end
    end
    
    if count < N
        warning('Only generated %d/%d points after %d attempts', count, N, attempts);
        points = points(1:count, :);
    end
end

function inside = isInsideMesh(p, V, F)
    % Ray casting test
    ray_dir = [1, 0, 0];
    count = 0;
    
    for i = 1:size(F, 1)
        v0 = V(F(i,1), :);
        v1 = V(F(i,2), :);
        v2 = V(F(i,3), :);
        
        if rayTriangleIntersect(p, ray_dir, v0, v1, v2)
            count = count + 1;
        end
    end
    
    inside = mod(count, 2) == 1;
end

function hit = rayTriangleIntersect(orig, dir, v0, v1, v2)
    % Möller-Trumbore algorithm
    edge1 = v1 - v0;
    edge2 = v2 - v0;
    h = cross(dir, edge2);
    a = dot(edge1, h);
    
    if abs(a) < 1e-10
        hit = false;
        return;
    end
    
    f = 1 / a;
    s = orig - v0;
    u = f * dot(s, h);
    
    if u < 0 || u > 1
        hit = false;
        return;
    end
    
    q = cross(s, edge1);
    v = f * dot(dir, q);
    
    if v < 0 || u + v > 1
        hit = false;
        return;
    end
    
    t = f * dot(edge2, q);
    hit = t > 1e-10;
end
