%% Setup Script: Preprocess Asteroids for Polyhedral Gravity Simulation
% This script must be run ONCE before the two-body simulation
% It creates preprocessed data structures for shape-based gravity

clear all; close all; clc;

fprintf('=== Asteroid Preprocessing for Polyhedral Gravity ===\n\n');

%% Configuration
% Tetrahedralization mesh size (adjust for performance vs accuracy trade-off)
% Larger Hmax = fewer tetrahedra = faster but less accurate
% Smaller Hmax = more tetrahedra = slower but more accurate
Hmax_bennu = 30;    % meters (for Bennu, ~500m diameter)
Hmax_itokawa = 40;  % meters (for Itokawa, ~600m size)

fprintf('Mesh parameters:\n');
fprintf('  Bennu Hmax:   %.1f m\n', Hmax_bennu);
fprintf('  Itokawa Hmax: %.1f m\n', Hmax_itokawa);
fprintf('\n');

%% Load and preprocess Bennu
fprintf('========================================\n');
fprintf('BENNU PREPROCESSING\n');
fprintf('========================================\n');

% Load Bennu mesh
bennu_raw = loadBennuMesh();

% Precompute polyhedral gravity geometry
polyBennu = precomputePolyhedron(bennu_raw.vertices, bennu_raw.faces, bennu_raw.density * 1000);

% Tetrahedralize volume for mass distribution
try
    [bennu_raw.C_body, bennu_raw.V_tet, bennu_raw.m_tet] = ...
        tetrahedralizeMesh(bennu_raw.vertices, bennu_raw.faces, bennu_raw.density * 1000, Hmax_bennu);
    
    fprintf('\n✓ Bennu tetrahedralization successful!\n');
    
catch ME
    fprintf('ℹ Bennu tetrahedralization not available (PDE Toolbox required)\n');
    fprintf('Falling back to simple mass point distribution...\n');
    
    % Fallback: distribute mass points uniformly
    N_points = 200;
    bennu_raw.C_body = generateMassPoints(bennu_raw.vertices, bennu_raw.faces, N_points);
    bennu_raw.m_tet = ones(N_points, 1) * bennu_raw.mass / N_points;
    bennu_raw.V_tet = bennu_raw.m_tet / (bennu_raw.density * 1000);
end

fprintf('\n');

%% Load and preprocess Itokawa
fprintf('========================================\n');
fprintf('ITOKAWA PREPROCESSING\n');
fprintf('========================================\n');

% Load Itokawa mesh
itokawa_raw = loadItokawaMesh();

% Precompute polyhedral gravity geometry
polyItokawa = precomputePolyhedron(itokawa_raw.vertices, itokawa_raw.faces, itokawa_raw.density * 1000);

% Tetrahedralize volume for mass distribution
try
    [itokawa_raw.C_body, itokawa_raw.V_tet, itokawa_raw.m_tet] = ...
        tetrahedralizeMesh(itokawa_raw.vertices, itokawa_raw.faces, itokawa_raw.density * 1000, Hmax_itokawa);
    
    fprintf('\n✓ Itokawa tetrahedralization successful!\n');
    
catch ME
    fprintf('ℹ Itokawa tetrahedralization not available (PDE Toolbox required)\n');
    fprintf('Falling back to simple mass point distribution...\n');
    
    % Fallback: distribute mass points uniformly
    N_points = 200;
    itokawa_raw.C_body = generateMassPoints(itokawa_raw.vertices, itokawa_raw.faces, N_points);
    itokawa_raw.m_tet = ones(N_points, 1) * itokawa_raw.mass / N_points;
    itokawa_raw.V_tet = itokawa_raw.m_tet / (itokawa_raw.density * 1000);
end

fprintf('\n');

%% Save preprocessed data
fprintf('========================================\n');
fprintf('SAVING PREPROCESSED DATA\n');
fprintf('========================================\n');

save('asteroid_poly_data.mat', 'polyBennu', 'polyItokawa', 'bennu_raw', 'itokawa_raw');

fprintf('✓ Saved to asteroid_poly_data.mat\n');
fprintf('  polyBennu:    Polyhedral gravity structure\n');
fprintf('  polyItokawa:  Polyhedral gravity structure\n');
fprintf('  bennu_raw:    Mass distribution (C_body, m_tet)\n');
fprintf('  itokawa_raw:  Mass distribution (C_body, m_tet)\n');

fprintf('\n========================================\n');
fprintf('PREPROCESSING COMPLETE\n');
fprintf('========================================\n');
fprintf('You can now run two_body_poly_simulation.m\n\n');

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

function data = loadItokawaMesh()
    % Load Itokawa from existing data file
    filename = 'vertex/ver512q.tab';
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open Itokawa file: %s', filename);
    end
    
    header_line = fgetl(fid);
    header_values = sscanf(header_line, '%d');
    num_vertices = header_values(1);
    num_faces = header_values(2);
    
    vertices = zeros(num_vertices, 3);
    faces = zeros(num_faces, 3);
    
    for i = 1:num_vertices
        line = fgetl(fid);
        values = sscanf(line, '%d %f %f %f');
        vertices(i, :) = values(2:4)';
    end
    
    for i = 1:num_faces
        line = fgetl(fid);
        values = sscanf(line, '%d %d %d %d');
        faces(i, :) = values(2:4)';
    end
    
    fclose(fid);
    
    vertices = vertices * 1000;  % km to m
    vertices = vertices - mean(vertices, 1);
    
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
    % Fallback: generate N points inside the mesh using rejection sampling
    fprintf('  Generating %d mass points via rejection sampling...\n', N);
    
    % Bounding box
    bbox_min = min(V, [], 1);
    bbox_max = max(V, [], 1);
    
    points = zeros(N, 3);
    count = 0;
    
    while count < N
        % Random point in bounding box
        p = bbox_min + rand(1, 3) .* (bbox_max - bbox_min);
        
        % Check if inside mesh (simple ray casting)
        if isInsideMesh(p, V, F)
            count = count + 1;
            points(count, :) = p;
        end
    end
    
    fprintf('  ✓ Generated %d points\n', N);
end

function inside = isInsideMesh(p, V, F)
    % Simple inside test: count ray intersections
    % Ray from p in +X direction
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
    % Möller-Trumbore ray-triangle intersection
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
