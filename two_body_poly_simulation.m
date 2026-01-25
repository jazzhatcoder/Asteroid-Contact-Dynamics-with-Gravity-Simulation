%% Two-Body Asteroid Simulation with Shape-Based Polyhedral Gravity
% Simulates Bennu and Itokawa with Werner-Scheeres polyhedral gravity model
% Shape-based gravity provides realistic torques and near-field accuracy

clear all; close all; clc;

fprintf('=== Two-Body Simulation with Polyhedral Gravity ===\n\n');

%% Load preprocessed polyhedral gravity data
fprintf('Loading preprocessed polyhedral gravity data...\n');

% Try fast version first (optimized mesh sizes)
if exist('asteroid_poly_data_fast.mat', 'file')
    fprintf('Found FAST preprocessed data (optimized)\n');
    load('asteroid_poly_data_fast.mat', 'polyBennu', 'polyItokawa', 'bennu_raw', 'itokawa_raw');
elseif exist('asteroid_poly_data.mat', 'file')
    fprintf('Found standard preprocessed data\n');
    load('asteroid_poly_data.mat', 'polyBennu', 'polyItokawa', 'bennu_raw', 'itokawa_raw');
else
    error(['No preprocessed data found!\n' ...
           'Run setup_poly_gravity_fast.m first to create optimized meshes.']);
end

fprintf('✓ Loaded polyhedral gravity structures\n');
fprintf('  Bennu:   %d vertices, %d faces, %d edges, %d mass points\n', ...
    size(bennu_raw.vertices, 1), size(bennu_raw.faces, 1), ...
    size(polyBennu.edges, 1), length(bennu_raw.m_tet));
fprintf('  Itokawa: %d vertices, %d faces, %d edges, %d mass points\n', ...
    size(itokawa_raw.vertices, 1), size(itokawa_raw.faces, 1), ...
    size(polyItokawa.edges, 1), length(itokawa_raw.m_tet));

%% Create asteroid structures
bennu = struct();
bennu.name = 'Bennu';
bennu.mass = bennu_raw.mass;
bennu.density = bennu_raw.density;
bennu.volume = bennu_raw.volume;
bennu.Rbound = bennu_raw.Rbound;
bennu.vertices = bennu_raw.vertices;
bennu.faces = bennu_raw.faces;
bennu.I_body = calculateInertiaTensor(bennu_raw.vertices, bennu_raw.faces, bennu_raw.density);
bennu.poly = polyBennu;
bennu.C_body = bennu_raw.C_body;  % Mass point locations (body frame)
bennu.m_tet = bennu_raw.m_tet;    % Mass point masses

itokawa = struct();
itokawa.name = 'Itokawa';
itokawa.mass = itokawa_raw.mass;
itokawa.density = itokawa_raw.density;
itokawa.volume = itokawa_raw.volume;
itokawa.Rbound = itokawa_raw.Rbound;
itokawa.vertices = itokawa_raw.vertices;
itokawa.faces = itokawa_raw.faces;
itokawa.I_body = calculateInertiaTensor(itokawa_raw.vertices, itokawa_raw.faces, itokawa_raw.density);
itokawa.poly = polyItokawa;
itokawa.C_body = itokawa_raw.C_body;  % Mass point locations (body frame)
itokawa.m_tet = itokawa_raw.m_tet;    % Mass point masses

fprintf('\nBennu properties:\n');
fprintf('  Mass:   %.4e kg\n', bennu.mass);
fprintf('  Volume: %.4e m³\n', bennu.volume);
fprintf('  Rbound: %.1f m\n', bennu.Rbound);

fprintf('\nItokawa properties:\n');
fprintf('  Mass:   %.4e kg\n', itokawa.mass);
fprintf('  Volume: %.4e m³\n', itokawa.volume);
fprintf('  Rbound: %.1f m\n', itokawa.Rbound);

%% Set initial conditions for two-body problem
fprintf('\nSetting initial conditions...\n');

% Bennu initial state
bennu.r = [-500; 0; 0];        % Position [m] (left side)
bennu.v = [0.5; 0.2; 0];       % Velocity [m/s] (moving right and up)
bennu.q = [1; 0; 0; 0];        % Orientation (identity)
bennu.omega = [0.01; 0.02; 0.03]; % Angular velocity [rad/s] (tumbling)

% Itokawa initial state
itokawa.r = [500; 0; 0];       % Position [m] (right side)
itokawa.v = [-0.3; -0.1; 0];   % Velocity [m/s] (moving left and down)
itokawa.q = [1; 0; 0; 0];      % Orientation (identity)
itokawa.omega = [0.02; -0.01; 0.015]; % Angular velocity [rad/s] (tumbling)

fprintf('Initial separation: %.1f m\n', norm(itokawa.r - bennu.r));
fprintf('Initial relative velocity: %.3f m/s\n', norm(itokawa.v - bennu.v));

%% Physical constants
G = 6.67430e-11;  % Gravitational constant [m³/(kg·s²)]
e_restitution = 0.5;  % Coefficient of restitution (0 = inelastic, 1 = elastic)
mu_friction = 0.6;    % Coefficient of friction (typical for rocky surfaces)
v_resting_threshold = 0.01;  % Resting contact threshold [m/s]

%% Simulation parameters
dt = 0.5;              % Time step [s]
N_steps = 2000;        % Number of steps (FULL simulation)
t_total = N_steps * dt;  % Total simulation time [s]

fprintf('\nSimulation parameters:\n');
fprintf('  Gravitational constant: %.5e m³/(kg·s²)\n', G);
fprintf('  Coefficient of restitution: %.2f\n', e_restitution);
fprintf('  Coefficient of friction: %.2f\n', mu_friction);
fprintf('  Time step:              %.2f s\n', dt);
fprintf('  Total steps:            %d\n', N_steps);
fprintf('  Total time:             %.1f s (%.2f minutes)\n', t_total, t_total/60);

% Check if mass points are already reasonable
fprintf('\nMass point configuration:\n');
fprintf('  Bennu mass points:   %d\n', length(bennu.m_tet));
fprintf('  Itokawa mass points: %d\n', length(itokawa.m_tet));

% Only subsample if needed (if more than 100 points per body)
if length(bennu.m_tet) > 100 || length(itokawa.m_tet) > 100
    fprintf('  ⚠ Mass points > 100, subsampling for performance...\n');
    
    subsample_bennu = max(1, floor(length(bennu.m_tet) / 50));
    subsample_itokawa = max(1, floor(length(itokawa.m_tet) / 50));

    bennu.C_body = bennu.C_body(1:subsample_bennu:end, :);
    bennu.m_tet = bennu.m_tet(1:subsample_bennu:end);
    bennu.m_tet = bennu.m_tet * (bennu.mass / sum(bennu.m_tet));

    itokawa.C_body = itokawa.C_body(1:subsample_itokawa:end, :);
    itokawa.m_tet = itokawa.m_tet(1:subsample_itokawa:end);
    itokawa.m_tet = itokawa.m_tet * (itokawa.mass / sum(itokawa.m_tet));

    fprintf('  After subsampling: Bennu %d, Itokawa %d\n', ...
        length(bennu.m_tet), length(itokawa.m_tet));
else
    fprintf('  ✓ Mass points already optimized, no subsampling needed\n');
end

%% Pre-allocate trajectory arrays
trajectory = struct();
trajectory.t = zeros(N_steps, 1);

% Bennu trajectory
trajectory.bennu.r = zeros(N_steps, 3);
trajectory.bennu.v = zeros(N_steps, 3);
trajectory.bennu.q = zeros(N_steps, 4);
trajectory.bennu.omega = zeros(N_steps, 3);
trajectory.bennu.F_grav = zeros(N_steps, 3);
trajectory.bennu.tau_grav = zeros(N_steps, 3);

% Itokawa trajectory
trajectory.itokawa.r = zeros(N_steps, 3);
trajectory.itokawa.v = zeros(N_steps, 3);
trajectory.itokawa.q = zeros(N_steps, 4);
trajectory.itokawa.omega = zeros(N_steps, 3);
trajectory.itokawa.F_grav = zeros(N_steps, 3);
trajectory.itokawa.tau_grav = zeros(N_steps, 3);

% System properties
trajectory.separation = zeros(N_steps, 1);
trajectory.contact = false(N_steps, 1);
trajectory.contact_normal = zeros(N_steps, 3);
trajectory.contact_point = zeros(N_steps, 3);
trajectory.bennu_lever = zeros(N_steps, 3);
trajectory.itokawa_lever = zeros(N_steps, 3);
trajectory.vrel = zeros(N_steps, 3);
trajectory.vn = zeros(N_steps, 1);
trajectory.vt = zeros(N_steps, 3);
trajectory.closing = false(N_steps, 1);
trajectory.Jn = zeros(N_steps, 1);
trajectory.Jt = zeros(N_steps, 1);
trajectory.sticking = false(N_steps, 1);
trajectory.E_total = zeros(N_steps, 1);
trajectory.L_total = zeros(N_steps, 3);

fprintf('\nStarting integration...\n');
fprintf('Using shape-based polyhedral gravity (Werner-Scheeres model)\n');
fprintf('Progress: [');
for i=1:50; fprintf(' '); end
fprintf('] 0%%\n');

tic_start = tic;

%% Time integration loop with polyhedral gravity
% Governing equations:
% Translation: ṙ = v, v̇ = F/m (gravitational force from distributed mass)
% Rotation: I·ω̇ = τ_grav - ω × (I·ω) (Euler's equations WITH gravitational torque!)
% Quaternion: q̇ = 0.5 * Ω(ω) * q

for k = 1:N_steps
    % Store current state
    trajectory.t(k) = (k-1) * dt;
    
    trajectory.bennu.r(k,:) = bennu.r';
    trajectory.bennu.v(k,:) = bennu.v';
    trajectory.bennu.q(k,:) = bennu.q';
    trajectory.bennu.omega(k,:) = bennu.omega';
    
    trajectory.itokawa.r(k,:) = itokawa.r';
    trajectory.itokawa.v(k,:) = itokawa.v';
    trajectory.itokawa.q(k,:) = itokawa.q';
    trajectory.itokawa.omega(k,:) = itokawa.omega';
    
    % === STEP 4.1: COMPUTE MUTUAL GRAVITATIONAL FORCE AND TORQUE (SHAPE-BASED) ===
    
    % Separation vector (from bennu to itokawa)
    r12 = itokawa.r - bennu.r;
    d = norm(r12);
    
    % Store separation distance
    trajectory.separation(k) = d;
    
    % Compute gravitational force and torque using polyhedral gravity model
    % Force on Bennu due to Itokawa's shape
    [F_on_bennu, tau_on_bennu] = mutualForceTorque(...
        itokawa.poly, ...          % Source polyhedron (Itokawa)
        itokawa.r, itokawa.q, ...  % Source pose
        bennu.r, bennu.q, ...      % Target pose (Bennu)
        bennu.C_body, bennu.m_tet);  % Target mass distribution
    
    % Force on Itokawa due to Bennu's shape
    [F_on_itokawa, tau_on_itokawa] = mutualForceTorque(...
        bennu.poly, ...              % Source polyhedron (Bennu)
        bennu.r, bennu.q, ...        % Source pose
        itokawa.r, itokawa.q, ...    % Target pose (Itokawa)
        itokawa.C_body, itokawa.m_tet);  % Target mass distribution
    
    % Store forces and torques
    trajectory.bennu.F_grav(k,:) = F_on_bennu';
    trajectory.bennu.tau_grav(k,:) = tau_on_bennu';
    trajectory.itokawa.F_grav(k,:) = F_on_itokawa';
    trajectory.itokawa.tau_grav(k,:) = tau_on_itokawa';
    
    % === STEP 5.1: DETECT CONTACT (BROAD PHASE) ===
    % Check if bounding spheres overlap
    if d <= (bennu.Rbound + itokawa.Rbound)
        trajectory.contact(k) = true;
        
        % === STEP 5.2: DEFINE CONTACT GEOMETRY ===
        % Contact normal (from bennu → itokawa)
        if d > 0
            contact_normal = r12 / d;
        else
            contact_normal = [1; 0; 0];  % Default if exactly overlapping
        end
        
        % Contact point - use actual mesh geometry
        R_bennu_current = quat2rotm_custom(bennu.q);
        vertices_world_bennu = (R_bennu_current * bennu.vertices')' + bennu.r';
        projections_bennu = vertices_world_bennu * contact_normal;
        [~, idx_bennu] = max(projections_bennu);
        contact_point_bennu = vertices_world_bennu(idx_bennu, :)';
        
        R_itokawa_current = quat2rotm_custom(itokawa.q);
        vertices_world_itokawa = (R_itokawa_current * itokawa.vertices')' + itokawa.r';
        projections_itokawa = vertices_world_itokawa * (-contact_normal);
        [~, idx_itokawa] = max(projections_itokawa);
        contact_point_itokawa = vertices_world_itokawa(idx_itokawa, :)';
        
        contact_point = 0.5 * (contact_point_bennu + contact_point_itokawa);
        
        % Lever arms
        r1c = contact_point - bennu.r;
        r2c = contact_point - itokawa.r;
        
        % Store contact geometry
        trajectory.contact_normal(k,:) = contact_normal';
        trajectory.contact_point(k,:) = contact_point';
        trajectory.bennu_lever(k,:) = r1c';
        trajectory.itokawa_lever(k,:) = r2c';
        
        % === STEP 5.3: RELATIVE VELOCITY AT CONTACT POINT ===
        vp1 = bennu.v + cross(bennu.omega, r1c);
        vp2 = itokawa.v + cross(itokawa.omega, r2c);
        vrel = vp2 - vp1;
        
        % === STEP 6.1: DECOMPOSE RELATIVE VELOCITY ===
        vn = dot(vrel, contact_normal);
        vt = vrel - vn * contact_normal;
        
        % === STEP 6.2: COMPUTE TANGENTIAL DIRECTION ===
        if norm(vt) > 1e-8
            t = vt / norm(vt);
        else
            t = [0; 0; 0];
        end
        
        trajectory.vrel(k,:) = vrel';
        trajectory.vn(k) = vn;
        trajectory.vt(k,:) = vt';
        
        % === STEP 6.3: CANDIDATE TANGENTIAL IMPULSE ===
        Jt_star = 0;
        
        if norm(vt) > 1e-8
            R_bennu = quat2rotm_custom(bennu.q);
            R_itokawa = quat2rotm_custom(itokawa.q);
            I_world_bennu = R_bennu * bennu.I_body * R_bennu';
            I_world_itokawa = R_itokawa * itokawa.I_body * R_itokawa';
            
            r1_cross_t = cross(r1c, t);
            r2_cross_t = cross(r2c, t);
            
            term1_t = cross(I_world_bennu \ r1_cross_t, r1c);
            term2_t = cross(I_world_itokawa \ r2_cross_t, r2c);
            
            denominator_t = (1/bennu.mass) + (1/itokawa.mass) + dot(t, term1_t + term2_t);
            Jt_star = -norm(vt) / denominator_t;
        end
        
        % === STEP 6.6: RESTING CONTACT DETECTION ===
        is_resting = abs(vn) < v_resting_threshold;
        
        % Check if bodies are closing or in resting contact
        if vn < 0 || (is_resting && vn < v_resting_threshold)
            trajectory.closing(k) = true;
            
            % === STEP 5.4: COMPUTE NORMAL IMPULSE ===
            R_bennu = quat2rotm_custom(bennu.q);
            R_itokawa = quat2rotm_custom(itokawa.q);
            I_world_bennu = R_bennu * bennu.I_body * R_bennu';
            I_world_itokawa = R_itokawa * itokawa.I_body * R_itokawa';
            
            r1_cross_n = cross(r1c, contact_normal);
            r2_cross_n = cross(r2c, contact_normal);
            
            term1 = cross(I_world_bennu \ r1_cross_n, r1c);
            term2 = cross(I_world_itokawa \ r2_cross_n, r2c);
            
            denominator = (1/bennu.mass) + (1/itokawa.mass) + dot(contact_normal, term1 + term2);
            
            if is_resting
                e_contact = 0.0;
            else
                e_contact = e_restitution;
            end
            
            Jn = -(1 + e_contact) * vn / denominator;
            trajectory.Jn(k) = Jn;
            
            % === STEP 6.4: APPLY COULOMB FRICTION LIMIT ===
            Jt = 0;
            
            if norm(vt) > 1e-8
                Jt_max = mu_friction * Jn;
                
                if abs(Jt_star) <= Jt_max
                    Jt = Jt_star;
                    trajectory.sticking(k) = true;
                else
                    Jt = sign(Jt_star) * Jt_max;
                    trajectory.sticking(k) = false;
                end
                
                trajectory.Jt(k) = Jt;
            end
            
            % === STEP 6.5: APPLY COMBINED IMPULSE ===
            if norm(vt) > 1e-8
                Jt_vec = Jt * t;
            else
                Jt_vec = [0; 0; 0];
            end
            
            Jn_vec = Jn * contact_normal;
            J_total = Jn_vec + Jt_vec;
            
            % Linear velocity updates
            bennu.v = bennu.v - J_total / bennu.mass;
            itokawa.v = itokawa.v + J_total / itokawa.mass;
            
            % Angular velocity updates
            dL_world_bennu = cross(r1c, J_total);
            dL_world_itokawa = cross(r2c, J_total);
            
            dL_body_bennu = R_bennu' * dL_world_bennu;
            dL_body_itokawa = R_itokawa' * dL_world_itokawa;
            
            bennu.omega = bennu.omega - (bennu.I_body \ dL_body_bennu);
            itokawa.omega = itokawa.omega + (itokawa.I_body \ dL_body_itokawa);
        end
    end
    
    % === COMPUTE TOTAL ENERGY AND ANGULAR MOMENTUM ===
    
    % Kinetic energy
    KE_bennu = 0.5 * bennu.mass * (bennu.v' * bennu.v);
    KE_itokawa = 0.5 * itokawa.mass * (itokawa.v' * itokawa.v);
    
    % Rotational energy
    RE_bennu = 0.5 * bennu.omega' * bennu.I_body * bennu.omega;
    RE_itokawa = 0.5 * itokawa.omega' * itokawa.I_body * itokawa.omega;
    
    % Gravitational potential energy (approximate as point-mass at this separation)
    PE = -G * bennu.mass * itokawa.mass / d;
    
    % Total energy
    trajectory.E_total(k) = KE_bennu + KE_itokawa + RE_bennu + RE_itokawa + PE;
    
    % Total angular momentum
    r_com = (bennu.mass * bennu.r + itokawa.mass * itokawa.r) / (bennu.mass + itokawa.mass);
    L_orbital_bennu = cross(bennu.r - r_com, bennu.mass * bennu.v);
    L_orbital_itokawa = cross(itokawa.r - r_com, itokawa.mass * itokawa.v);
    
    R_bennu = quat2rotm_custom(bennu.q);
    R_itokawa = quat2rotm_custom(itokawa.q);
    L_spin_bennu = R_bennu * (bennu.I_body * bennu.omega);
    L_spin_itokawa = R_itokawa * (itokawa.I_body * itokawa.omega);
    
    trajectory.L_total(k,:) = (L_orbital_bennu + L_orbital_itokawa + L_spin_bennu + L_spin_itokawa)';
    
    % === INTEGRATION STEP ===
    
    % --- BENNU ---
    % Translational motion
    bennu.v = bennu.v + (F_on_bennu / bennu.mass) * dt;
    bennu.r = bennu.r + bennu.v * dt;
    
    % Rotational motion (WITH GRAVITATIONAL TORQUE!)
    % Transform gravitational torque to body frame
    tau_body_bennu = R_bennu' * tau_on_bennu;
    
    % Euler's equations: I·ω̇ = τ_grav - ω × (I·ω)
    omega_dot_bennu = bennu.I_body \ (tau_body_bennu - cross(bennu.omega, bennu.I_body * bennu.omega));
    bennu.omega = bennu.omega + omega_dot_bennu * dt;
    
    % Quaternion propagation
    wx = bennu.omega(1); wy = bennu.omega(2); wz = bennu.omega(3);
    Omega_bennu = [ 0   -wx  -wy  -wz;
                    wx   0    wz  -wy;
                    wy  -wz   0    wx;
                    wz   wy  -wx   0 ];
    q_dot_bennu = 0.5 * Omega_bennu * bennu.q;
    bennu.q = bennu.q + q_dot_bennu * dt;
    bennu.q = bennu.q / norm(bennu.q);
    
    % --- ITOKAWA ---
    % Translational motion
    itokawa.v = itokawa.v + (F_on_itokawa / itokawa.mass) * dt;
    itokawa.r = itokawa.r + itokawa.v * dt;
    
    % Rotational motion (WITH GRAVITATIONAL TORQUE!)
    % Transform gravitational torque to body frame
    tau_body_itokawa = R_itokawa' * tau_on_itokawa;
    
    % Euler's equations: I·ω̇ = τ_grav - ω × (I·ω)
    omega_dot_itokawa = itokawa.I_body \ (tau_body_itokawa - cross(itokawa.omega, itokawa.I_body * itokawa.omega));
    itokawa.omega = itokawa.omega + omega_dot_itokawa * dt;
    
    % Quaternion propagation
    wx = itokawa.omega(1); wy = itokawa.omega(2); wz = itokawa.omega(3);
    Omega_itokawa = [ 0   -wx  -wy  -wz;
                      wx   0    wz  -wy;
                      wy  -wz   0    wx;
                      wz   wy  -wx   0 ];
    q_dot_itokawa = 0.5 * Omega_itokawa * itokawa.q;
    itokawa.q = itokawa.q + q_dot_itokawa * dt;
    itokawa.q = itokawa.q / norm(itokawa.q);
    
    % Enhanced progress indicator with progress bar
    if mod(k, 50) == 0 || k == N_steps
        elapsed = toc(tic_start);
        pct = k / N_steps;
        eta = elapsed / pct - elapsed;
        
        % Progress bar
        fprintf('\rProgress: [');
        num_bars = floor(pct * 50);
        for i=1:num_bars; fprintf('='); end
        if num_bars < 50; fprintf('>'); end
        for i=(num_bars+2):50; fprintf(' '); end
        fprintf('] %.1f%% | Step %d/%d | Sep: %.1fm | Elapsed: %.1fs | ETA: %.1fs', ...
            pct*100, k, N_steps, d, elapsed, eta);
    end
end

fprintf('\n');

fprintf('Integration completed.\n');

%% Final state
fprintf('\nFinal state:\n');
fprintf('Bennu:\n');
fprintf('  Position:   [%.2f, %.2f, %.2f] m\n', bennu.r);
fprintf('  Velocity:   [%.3f, %.3f, %.3f] m/s\n', bennu.v);
fprintf('  Spin rate:  %.4f rad/s (%.2f deg/s)\n', norm(bennu.omega), norm(bennu.omega)*180/pi);

fprintf('Itokawa:\n');
fprintf('  Position:   [%.2f, %.2f, %.2f] m\n', itokawa.r);
fprintf('  Velocity:   [%.3f, %.3f, %.3f] m/s\n', itokawa.v);
fprintf('  Spin rate:  %.4f rad/s (%.2f deg/s)\n', norm(itokawa.omega), norm(itokawa.omega)*180/pi);

fprintf('\nFinal separation: %.1f m\n', trajectory.separation(end));

% Contact summary
num_contacts = sum(trajectory.contact);
num_closing = sum(trajectory.closing);
fprintf('\nCollision Summary:\n');
fprintf('  Contacts detected:  %d / %d timesteps (%.1f%%)\n', ...
    num_contacts, N_steps, 100*num_contacts/N_steps);
if num_closing > 0
    fprintf('  Collisions resolved: %d\n', num_closing);
    fprintf('  Sticking contacts:   %d (%.1f%%)\n', ...
        sum(trajectory.sticking), 100*sum(trajectory.sticking)/num_closing);
    fprintf('  Sliding contacts:    %d (%.1f%%)\n', ...
        num_closing - sum(trajectory.sticking), ...
        100*(num_closing - sum(trajectory.sticking))/num_closing);
end

%% Conservation verification
E_initial = trajectory.E_total(1);
E_final = trajectory.E_total(end);
L_initial = trajectory.L_total(1,:);
L_final = trajectory.L_total(end,:);

fprintf('\nConservation verification:\n');
fprintf('  Energy change:          %.6e J (%.3f%%)\n', ...
    E_final - E_initial, 100*abs(E_final - E_initial)/abs(E_initial));
fprintf('  Angular momentum error: %.6e kg·m²/s (%.3f%%)\n', ...
    norm(L_final - L_initial), 100*norm(L_final - L_initial)/norm(L_initial));

%% Visualization
fprintf('\nCreating visualizations...\n');

% Figure 1: Trajectories
figure('Name', 'Two-Body Polyhedral Gravity Simulation', 'Position', [50, 50, 1400, 900]);

% 3D trajectory with COM paths
subplot(2,2,1);
plot3(trajectory.bennu.r(:,1), trajectory.bennu.r(:,2), trajectory.bennu.r(:,3), ...
      'r-', 'LineWidth', 1.5, 'DisplayName', 'Bennu COM');
hold on;
plot3(trajectory.itokawa.r(:,1), trajectory.itokawa.r(:,2), trajectory.itokawa.r(:,3), ...
      'b-', 'LineWidth', 1.5, 'DisplayName', 'Itokawa COM');
plot3(trajectory.bennu.r(1,1), trajectory.bennu.r(1,2), trajectory.bennu.r(1,3), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot3(trajectory.itokawa.r(1,1), trajectory.itokawa.r(1,2), trajectory.itokawa.r(1,3), ...
      'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('COM Trajectories');
legend('Location', 'best');
axis equal;
view(45, 30);

% Separation vs time
subplot(2,2,2);
plot(trajectory.t, trajectory.separation, 'k-', 'LineWidth', 2);
hold on;
yline(bennu.Rbound + itokawa.Rbound, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Contact Threshold');
if any(trajectory.contact)
    contact_times = trajectory.t(trajectory.contact);
    contact_seps = trajectory.separation(trajectory.contact);
    scatter(contact_times, contact_seps, 50, 'r', 'filled', 'DisplayName', 'Contact');
end
grid on;
xlabel('Time (s)');
ylabel('Separation (m)');
title('Distance Between Asteroids');
legend('Location', 'best');

% Spin rates
subplot(2,2,3);
spin_bennu = vecnorm(trajectory.bennu.omega, 2, 2);
spin_itokawa = vecnorm(trajectory.itokawa.omega, 2, 2);

plot(trajectory.t, spin_bennu * 180/pi, 'r-', 'LineWidth', 2, 'DisplayName', 'Bennu');
hold on;
plot(trajectory.t, spin_itokawa * 180/pi, 'b-', 'LineWidth', 2, 'DisplayName', 'Itokawa');
grid on;
xlabel('Time (s)');
ylabel('Spin Rate (deg/s)');
title('Angular Velocity (Gravity + Friction Effects)');
legend('Location', 'best');

% Energy
subplot(2,2,4);
plot(trajectory.t, (trajectory.E_total - E_initial)/abs(E_initial) * 100, 'g-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Energy Change (%)');
title('Total Energy');

sgtitle('Polyhedral Gravity Simulation: Shape-Based Dynamics', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 2: Gravitational torques
figure('Name', 'Gravitational Torques', 'Position', [100, 100, 1200, 600]);

subplot(2,2,1);
plot(trajectory.t, vecnorm(trajectory.bennu.F_grav, 2, 2), 'r-', 'LineWidth', 2);
hold on;
plot(trajectory.t, vecnorm(trajectory.itokawa.F_grav, 2, 2), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Force (N)');
title('Gravitational Force Magnitude');
legend('On Bennu', 'On Itokawa');

subplot(2,2,2);
plot(trajectory.t, vecnorm(trajectory.bennu.tau_grav, 2, 2), 'r-', 'LineWidth', 2);
hold on;
plot(trajectory.t, vecnorm(trajectory.itokawa.tau_grav, 2, 2), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Torque (N·m)');
title('Gravitational Torque Magnitude (NEW!)');
legend('On Bennu', 'On Itokawa');

subplot(2,2,3);
semilogy(trajectory.t, vecnorm(trajectory.bennu.F_grav, 2, 2), 'r-', 'LineWidth', 2);
hold on;
semilogy(trajectory.t, vecnorm(trajectory.itokawa.F_grav, 2, 2), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Force (N)');
title('Gravitational Force (log scale)');
legend('On Bennu', 'On Itokawa');

subplot(2,2,4);
semilogy(trajectory.t, max(vecnorm(trajectory.bennu.tau_grav, 2, 2), 1e-20), 'r-', 'LineWidth', 2);
hold on;
semilogy(trajectory.t, max(vecnorm(trajectory.itokawa.tau_grav, 2, 2), 1e-20), 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Torque (N·m)');
title('Gravitational Torque (log scale)');
legend('On Bennu', 'On Itokawa');

%% Figure 3: Full 3D Mesh Animation
fprintf('\nCreating 3D mesh animation...\n');

% Animation parameters
frame_skip = max(1, floor(N_steps / 200));  % Show ~200 frames max
frame_indices = 1:frame_skip:N_steps;
if frame_indices(end) ~= N_steps
    frame_indices = [frame_indices, N_steps];  % Include final frame
end

% Create animation figure
fig_anim = figure('Name', '3D Asteroid Collision Animation', 'Position', [150, 50, 1200, 900]);

% Setup video writer
video_filename = 'asteroid_collision_animation.mp4';
video_writer = VideoWriter(video_filename, 'MPEG-4');
video_writer.FrameRate = 30;  % 30 fps
video_writer.Quality = 95;     % High quality
open(video_writer);
fprintf('  Saving animation to: %s\n', video_filename);

fprintf('  Animating %d frames...\\n', length(frame_indices));
fprintf('  Progress: [');
for i=1:50; fprintf(' '); end
fprintf('] 0%%\n');

for frame_idx = 1:length(frame_indices)
    k = frame_indices(frame_idx);
    
    % Progress bar for animation
    if mod(frame_idx, max(1, floor(length(frame_indices)/50))) == 0 || frame_idx == length(frame_indices)
        pct = frame_idx / length(frame_indices);
        fprintf('\r  Progress: [');
        num_bars = floor(pct * 50);
        for i=1:num_bars; fprintf('='); end
        if num_bars < 50; fprintf('>'); end
        for i=(num_bars+2):50; fprintf(' '); end
        fprintf('] %.1f%% | Frame %d/%d', pct*100, frame_idx, length(frame_indices));
    end
    
    clf;
    
    % Get current state
    r_bennu = trajectory.bennu.r(k, :)';
    q_bennu = trajectory.bennu.q(k, :)';
    r_itokawa = trajectory.itokawa.r(k, :)';
    q_itokawa = trajectory.itokawa.q(k, :)';
    
    % Compute rotation matrices from quaternions
    R_bennu = quat2rotm_custom(q_bennu);
    R_itokawa = quat2rotm_custom(q_itokawa);
    
    % Transform Bennu mesh vertices to world frame
    V_bennu_world = (R_bennu * bennu.vertices' + r_bennu)';
    
    % Transform Itokawa mesh vertices to world frame
    V_itokawa_world = (R_itokawa * itokawa.vertices' + r_itokawa)';
    
    % Plot Bennu mesh
    patch('Vertices', V_bennu_world, 'Faces', bennu.faces, ...
          'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', ...
          'AmbientStrength', 0.3, 'DiffuseStrength', 0.6, ...
          'SpecularStrength', 0.3, 'SpecularExponent', 10);
    hold on;
    
    % Plot Itokawa mesh
    patch('Vertices', V_itokawa_world, 'Faces', itokawa.faces, ...
          'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', ...
          'AmbientStrength', 0.3, 'DiffuseStrength', 0.6, ...
          'SpecularStrength', 0.3, 'SpecularExponent', 10);
    
    % Plot COM points
    plot3(r_bennu(1), r_bennu(2), r_bennu(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot3(r_itokawa(1), r_itokawa(2), r_itokawa(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    % Plot trajectories up to current time
    plot3(trajectory.bennu.r(1:k,1), trajectory.bennu.r(1:k,2), trajectory.bennu.r(1:k,3), ...
          'r-', 'LineWidth', 1.5, 'Color', [1 0 0 0.3]);
    plot3(trajectory.itokawa.r(1:k,1), trajectory.itokawa.r(1:k,2), trajectory.itokawa.r(1:k,3), ...
          'b-', 'LineWidth', 1.5, 'Color', [0 0 1 0.3]);
    
    % If in contact, show contact point and normal
    if trajectory.contact(k)
        cp = trajectory.contact_point(k, :);
        cn = trajectory.contact_normal(k, :);
        plot3(cp(1), cp(2), cp(3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
        quiver3(cp(1), cp(2), cp(3), cn(1)*50, cn(2)*50, cn(3)*50, ...
                'g', 'LineWidth', 3, 'MaxHeadSize', 2);
        text(cp(1), cp(2), cp(3)+100, 'CONTACT', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Lighting and appearance
    light('Position', [1 1 1], 'Style', 'infinite');
    light('Position', [-1 -1 -1], 'Style', 'infinite');
    
    % Axis settings
    axis equal;
    grid on;
    xlabel('X (m)', 'FontSize', 12);
    ylabel('Y (m)', 'FontSize', 12);
    zlabel('Z (m)', 'FontSize', 12);
    
    % Title with time and separation info
    title_str = sprintf('Time: %.1f s | Separation: %.1f m | Energy Change: %.2f%%', ...
                       trajectory.t(k), trajectory.separation(k), ...
                       100*(trajectory.E_total(k) - E_initial)/abs(E_initial));
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    
    % Set view angle
    view(45, 20);
    
    % Adjust axis limits to show both asteroids with some margin
    all_verts = [V_bennu_world; V_itokawa_world];
    x_range = [min(all_verts(:,1)), max(all_verts(:,1))];
    y_range = [min(all_verts(:,2)), max(all_verts(:,2))];
    z_range = [min(all_verts(:,3)), max(all_verts(:,3))];
    
    margin = 200;  % 200 m margin
    xlim([x_range(1)-margin, x_range(2)+margin]);
    ylim([y_range(1)-margin, y_range(2)+margin]);
    zlim([z_range(1)-margin, z_range(2)+margin]);
    
    % Force rendering and capture frame
    drawnow;
    frame = getframe(fig_anim);
    writeVideo(video_writer, frame);
end

% Close video writer
close(video_writer);

fprintf('\n');
fprintf('Animation complete!\n');
fprintf('  ✓ Video saved as: %s\n', video_filename);
fprintf('  ✓ Figure shows final state - you can rotate the view with the mouse.\n');

fprintf('\nVisualization complete.\n');

%% Helper functions

function I_tensor = calculateInertiaTensor(vertices, faces, density_gcm3)
    % Calculate inertia tensor from mesh
    density = density_gcm3 * 1000;  % Convert to kg/m³
    int_x2 = 0; int_y2 = 0; int_z2 = 0;
    int_xy = 0; int_xz = 0; int_yz = 0;
    
    for i = 1:size(faces, 1)
        a = vertices(faces(i,1), :)';
        b = vertices(faces(i,2), :)';
        c = vertices(faces(i,3), :)';
        V_t = dot(a, cross(b, c)) / 6;
        
        ax = a(1); ay = a(2); az = a(3);
        bx = b(1); by = b(2); bz = b(3);
        cx = c(1); cy = c(2); cz = c(3);
        
        int_x2 = int_x2 + (V_t/10) * (ax^2 + bx^2 + cx^2 + ax*bx + bx*cx + cx*ax);
        int_y2 = int_y2 + (V_t/10) * (ay^2 + by^2 + cy^2 + ay*by + by*cy + cy*ay);
        int_z2 = int_z2 + (V_t/10) * (az^2 + bz^2 + cz^2 + az*bz + bz*cz + cz*az);
        int_xy = int_xy + (V_t/20) * (2*(ax*ay + bx*by + cx*cy) + ax*by + ay*bx + bx*cy + by*cx + cx*ay + cy*ax);
        int_xz = int_xz + (V_t/20) * (2*(ax*az + bx*bz + cx*cz) + ax*bz + az*bx + bx*cz + bz*cx + cx*az + cz*ax);
        int_yz = int_yz + (V_t/20) * (2*(ay*az + by*bz + cy*cz) + ay*bz + az*by + by*cz + bz*cy + cy*az + cz*ay);
    end
    
    I_tensor = zeros(3, 3);
    I_tensor(1,1) = density * (int_y2 + int_z2);
    I_tensor(2,2) = density * (int_x2 + int_z2);
    I_tensor(3,3) = density * (int_x2 + int_y2);
    I_tensor(1,2) = -density * int_xy; I_tensor(2,1) = I_tensor(1,2);
    I_tensor(1,3) = -density * int_xz; I_tensor(3,1) = I_tensor(1,3);
    I_tensor(2,3) = -density * int_yz; I_tensor(3,2) = I_tensor(2,3);
end

function R = quat2rotm_custom(q)
    % Convert quaternion to rotation matrix
    if size(q, 2) == 4
        q = q';
    end
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),     2*(x*z+w*y);
         2*(x*y+w*z),     1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),     2*(y*z+w*x),     1-2*(x^2+y^2)];
end
