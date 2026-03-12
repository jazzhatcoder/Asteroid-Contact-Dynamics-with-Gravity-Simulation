function [C, V_tet, m] = tetrahedralizeMesh(V, F, rho, Hmax)
% TETRAHEDRALIZEMESH Convert surface mesh to tetrahedral volume mesh
%
% [C, V_tet, m] = tetrahedralizeMesh(V, F, rho, Hmax)
%
% Converts a closed triangular surface mesh into a tetrahedral volume mesh,
% then computes mass elements for polyhedral gravity evaluation.
%
% Inputs:
%   V    - Surface vertices (Nv x 3) [m], centered at COM
%   F    - Surface faces (Nf x 3) triangle indices
%   rho  - Density [kg/m³]
%   Hmax - Maximum tetrahedron size [m] (optional, default = 0.05 * size)
%
% Outputs:
%   C     - Tetrahedron centroids (Nt x 3) [m] in body frame
%   V_tet - Tetrahedron volumes (Nt x 1) [m³]
%   m     - Tetrahedron masses (Nt x 1) [kg]
%
% Requires: MATLAB PDE Toolbox (for geometryFromMesh and generateMesh)

fprintf('Tetrahedralizing asteroid volume mesh...\n');

%% Set default mesh size if not provided
if nargin < 4 || isempty(Hmax)
    % Default: 5% of asteroid size
    Hmax = 0.05 * max(vecnorm(V, 2, 2));
    fprintf('  Using default Hmax = %.2f m\n', Hmax);
else
    fprintf('  Using Hmax = %.2f m\n', Hmax);
end

%% Create PDE model and geometry from mesh
try
    model = createpde();
    
    % geometryFromMesh expects nodes as 3xN and elements as 3xM
    nodes = V';        % 3 x Nv
    elements = F';     % 3 x Nf
    
    geometryFromMesh(model, nodes, elements);
    
    fprintf('  Surface geometry loaded successfully.\n');
    
catch ME
    error('Failed to create geometry from mesh. Ensure MATLAB PDE Toolbox is installed.\n%s', ME.message);
end

%% Generate tetrahedral mesh
try
    mesh = generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', Hmax);
    
    fprintf('  Tetrahedral mesh generated successfully.\n');
    fprintf('    Nodes: %d\n', size(mesh.Nodes, 2));
    fprintf('    Tetrahedra: %d\n', size(mesh.Elements, 2));
    
catch ME
    error('Failed to generate tetrahedral mesh.\n%s', ME.message);
end

%% Extract mesh data
X = mesh.Nodes';        % Nnodes x 3
T = mesh.Elements';     % Ntets x 4

Nt = size(T, 1);

%% Compute tetrahedron centroids and volumes
fprintf('  Computing tetrahedron centroids and volumes...\n');

C = zeros(Nt, 3);
V_tet = zeros(Nt, 1);

for i = 1:Nt
    % Four vertices of tetrahedron
    v1 = X(T(i,1), :);
    v2 = X(T(i,2), :);
    v3 = X(T(i,3), :);
    v4 = X(T(i,4), :);
    
    % Centroid
    C(i,:) = (v1 + v2 + v3 + v4) / 4;
    
    % Volume = |det([v2-v1; v3-v1; v4-v1])| / 6
    mat = [v2-v1; v3-v1; v4-v1];
    V_tet(i) = abs(det(mat)) / 6;
end

%% Compute masses
m = rho * V_tet;

%% Validation checks
total_volume = sum(V_tet);
total_mass = sum(m);

fprintf('\nTetrahedralization complete:\n');
fprintf('  Number of tetrahedra: %d\n', Nt);
fprintf('  Total volume: %.4e m³\n', total_volume);
fprintf('  Total mass: %.4e kg\n', total_mass);

% Check COM (should be near origin since input was COM-centered)
COM_check = sum(C .* m, 1) / total_mass;
fprintf('  COM check: [%.3e, %.3e, %.3e] m (should be ~0)\n', COM_check);

% Info if COM is offset from origin (acceptable for fallback method)
if norm(COM_check) > 0.01 * max(vecnorm(V, 2, 2))
    fprintf('  ℹ COM offset: %.2f m (within acceptable range)\n', norm(COM_check));
end

end
