function [F, tau] = mutualForceTorque(polyA, rA, qA, rB, qB, CB_body, mB)
% MUTUALFORCETORQUE Compute gravitational force and torque on body B due to polyhedron A
%
% [F, tau] = mutualForceTorque(polyA, rA, qA, rB, qB, CB_body, mB)
%
% Evaluates polyhedral gravity of A at mass points of B, then sums to get
% net force and torque on B.
%
% Inputs:
%   polyA   - Polyhedron structure for asteroid A (from precomputePolyhedron)
%   rA      - Position of A's COM in world frame (3x1) [m]
%   qA      - Orientation quaternion of A (4x1, scalar-first)
%   rB      - Position of B's COM in world frame (3x1) [m]
%   qB      - Orientation quaternion of B (4x1, scalar-first)
%   CB_body - Mass element centroids of B in body frame (N x 3) [m]
%   mB      - Mass element masses (N x 1) [kg]
%
% Outputs:
%   F   - Net gravitational force on B (3x1) [N]
%   tau - Net gravitational torque on B about its COM (3x1) [N·m]
%
% Note: Gravity field is computed in A's body frame, so we must transform
% all positions appropriately.

%% Transform B's mass points to A's body frame
% This is the key coordinate transformation for shape-shape gravity

% B's mass points in world frame
RB = quat2rotm_custom(qB);
CB_world = (RB * CB_body')' + rB';  % N x 3 (each row is a point)

% A's orientation
RA = quat2rotm_custom(qA);

% Transform to A's body frame (where polyA is defined)
% P_A_body = R_A' * (P_world - r_A)
CB_relative = CB_world - rA';  % N x 3
CB_in_A_frame = (RA' * CB_relative')';  % N x 3

%% Evaluate A's polyhedral gravity at B's mass points
% This is where the shape-shape magic happens!
gA = polyGravityAtPoints(polyA, CB_in_A_frame);  % N x 3 acceleration [m/s²]

% Transform acceleration back to world frame
gA_world = (RA * gA')';  % N x 3

%% Compute forces on each mass element
F_points = mB .* gA_world;  % N x 3 force vectors [N]

%% Sum net force on B
F = sum(F_points, 1)';  % 3x1 net force

%% Compute torque about B's COM
% Lever arm from B's COM to each mass point (in world frame)
lever = CB_world - rB';  % N x 3

% Torque = sum of r × F
tau = sum(cross(lever, F_points, 2), 1)';  % 3x1 net torque

end


function R = quat2rotm_custom(q)
% Custom quaternion to rotation matrix (scalar-first convention)
    if size(q, 2) == 4
        q = q';
    end
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),     2*(x*z+w*y);
         2*(x*y+w*z),     1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),     2*(y*z+w*x),     1-2*(x^2+y^2)];
end
