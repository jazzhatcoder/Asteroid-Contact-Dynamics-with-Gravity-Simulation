function g = polyGravityAtPoints(poly, P)
% POLYGRAVITYATPOINTS Evaluate polyhedral gravity at field points
%
% g = polyGravityAtPoints(poly, P)
%
% Computes gravitational acceleration using Werner-Scheeres polyhedral
% gravity model for constant-density polyhedra.
%
% Inputs:
%   poly - Polyhedron structure from precomputePolyhedron()
%   P    - Field points (Np x 3) [m] in world coordinates
%
% Output:
%   g    - Gravitational acceleration (Np x 3) [m/s²]
%
% Reference:
%   Werner & Scheeres (1997), "Exterior gravitation of a polyhedron
%   derived and compared with harmonic and mascon gravitation
%   representations of asteroid 4769 Castalia"
%   Celestial Mechanics and Dynamical Astronomy, 65, 313-344

G = 6.67430e-11;  % Gravitational constant [m³/(kg·s²)]

Np = size(P, 1);
g = zeros(Np, 3);

% Loop over all field points
for idx = 1:Np
    r = P(idx, :)';  % Field point (3x1)
    
    % Initialize accumulation
    edge_sum = zeros(3, 1);
    face_sum = zeros(3, 1);
    
    %% Edge contributions
    Ne = size(poly.edges, 1);
    
    for e = 1:Ne
        % Edge vertices
        v1 = poly.V(poly.edges(e,1), :)';
        v2 = poly.V(poly.edges(e,2), :)';
        
        % Vectors from field point to edge vertices
        r1 = v1 - r;
        r2 = v2 - r;
        
        % Edge vector
        e_vec = poly.edgeVecs(e, :)';
        
        % Distances
        r1_mag = norm(r1);
        r2_mag = norm(r2);
        
        % Log factor (edge integral)
        % L_e = ln((r1 + r2 + Le) / (r1 + r2 - Le))
        Le = poly.edgeLengths(e);
        
        % Avoid singularities
        if r1_mag < 1e-10 || r2_mag < 1e-10
            continue;
        end
        
        numerator = r1_mag + r2_mag + Le;
        denominator = r1_mag + r2_mag - Le;
        
        if denominator > 1e-10
            L_e = log(numerator / denominator);
        else
            L_e = 0;  % Skip if degenerate
        end
        
        % Edge dyad contribution
        % g += E_e * L_e (dyad is 3x3, L_e is scalar)
        E_e = squeeze(poly.edgeDyads(e, :, :));
        edge_sum = edge_sum + E_e * e_vec * L_e;
    end
    
    %% Face contributions
    Nf = size(poly.F, 1);
    
    for f = 1:Nf
        % Face vertices
        v1 = poly.V(poly.F(f,1), :)';
        v2 = poly.V(poly.F(f,2), :)';
        v3 = poly.V(poly.F(f,3), :)';
        
        % Vectors from field point to face vertices
        r1 = v1 - r;
        r2 = v2 - r;
        r3 = v3 - r;
        
        % Distances
        r1_mag = norm(r1);
        r2_mag = norm(r2);
        r3_mag = norm(r3);
        
        % Avoid singularities
        if r1_mag < 1e-10 || r2_mag < 1e-10 || r3_mag < 1e-10
            continue;
        end
        
        % Solid angle omega_f (Oosterom & Strackee formula)
        % omega = 2 * atan2(r1·(r2×r3), r1·r2·r3 + r1_mag*r2_mag*r3_mag + ...)
        
        numerator = dot(r1, cross(r2, r3));
        denominator = r1_mag * r2_mag * r3_mag + ...
                      r1_mag * dot(r2, r3) + ...
                      r2_mag * dot(r1, r3) + ...
                      r3_mag * dot(r1, r2);
        
        omega_f = 2 * atan2(numerator, denominator);
        
        % Face dyad contribution
        % g -= F_f * omega_f (dyad is 3x3, omega_f is scalar)
        F_f = squeeze(poly.faceDyads(f, :, :));
        n_f = poly.faceNormals(f, :)';
        face_sum = face_sum + F_f * n_f * omega_f;
    end
    
    %% Combine contributions
    % g(r) = G*rho * (edge_sum - face_sum)
    g(idx, :) = (G * poly.rho * (edge_sum - face_sum))';
end

end
