%% Two-Body Asteroid Simulation with Mutual Gravity
% Simulates Bennu and Itokawa with gravitational interaction
% Point-mass gravity acts at each asteroid's center of mass

clear all; close all; clc;

fprintf('=== Two-Body Asteroid Simulation with Mutual Gravity ===\n\n');

%% Load and initialize Bennu
fprintf('Loading Bennu asteroid...\n');
bennu = loadBennu();
fprintf('Bennu loaded: Mass = %.4e kg, Rbound = %.1f m\n', bennu.mass, bennu.Rbound);

%% Load and initialize Itokawa
fprintf('\nLoading Itokawa asteroid...\n');
itokawa = loadItokawa();
fprintf('Itokawa loaded: Mass = %.4e kg, Rbound = %.1f m\n', itokawa.mass, itokawa.Rbound);

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
v_resting_threshold = 0.01;  % Resting contact threshold [m/s] - below this, use contact dynamics

%% Simulation parameters
dt = 0.5;              % Time step [s]
N_steps = 2000;        % Number of steps
t_total = N_steps * dt;  % Total simulation time [s]

fprintf('\nSimulation parameters:\n');
fprintf('  Gravitational constant: %.5e m³/(kg·s²)\n', G);
fprintf('  Coefficient of restitution: %.2f\n', e_restitution);
fprintf('  Time step:              %.2f s\n', dt);
fprintf('  Total steps:            %d\n', N_steps);
fprintf('  Total time:             %.1f s (%.2f minutes)\n', t_total, t_total/60);

%% Pre-allocate trajectory arrays
trajectory = struct();
trajectory.t = zeros(N_steps, 1);

% Bennu trajectory
trajectory.bennu.r = zeros(N_steps, 3);
trajectory.bennu.v = zeros(N_steps, 3);
trajectory.bennu.q = zeros(N_steps, 4);
trajectory.bennu.omega = zeros(N_steps, 3);
trajectory.bennu.F = zeros(N_steps, 3);

% Itokawa trajectory
trajectory.itokawa.r = zeros(N_steps, 3);
trajectory.itokawa.v = zeros(N_steps, 3);
trajectory.itokawa.q = zeros(N_steps, 4);
trajectory.itokawa.omega = zeros(N_steps, 3);
trajectory.itokawa.F = zeros(N_steps, 3);

% System properties
trajectory.separation = zeros(N_steps, 1);
trajectory.contact = false(N_steps, 1);  % Broad-phase collision detection flag
trajectory.contact_normal = zeros(N_steps, 3);  % Contact normal (from bennu to itokawa)
trajectory.contact_point = zeros(N_steps, 3);   % Contact point in world frame
trajectory.bennu_lever = zeros(N_steps, 3);     % Lever arm from bennu COM to contact point
trajectory.itokawa_lever = zeros(N_steps, 3);   % Lever arm from itokawa COM to contact point
trajectory.vrel = zeros(N_steps, 3);            % Relative velocity at contact point
trajectory.vn = zeros(N_steps, 1);              % Normal component of relative velocity
trajectory.vt = zeros(N_steps, 3);              % Tangential component of relative velocity
trajectory.closing = false(N_steps, 1);         % True if bodies are closing (vn < 0)
trajectory.Jn = zeros(N_steps, 1);              % Normal impulse magnitude [N·s]
trajectory.Jt = zeros(N_steps, 1);              % Tangential impulse magnitude [N·s]
trajectory.sticking = false(N_steps, 1);        % True if contact sticks (friction sufficient)
trajectory.E_total = zeros(N_steps, 1);
trajectory.L_total = zeros(N_steps, 3);

fprintf('\nStarting integration...\n');

%% Time integration loop with mutual gravity
% Governing equations:
% Translation: ṙ = v, v̇ = F/m (gravitational force)
% Rotation: I·ω̇ + ω × (I·ω) = 0 (Euler's equations, no gravitational torque at COM)
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
    
    % === STEP 4.1: COMPUTE MUTUAL GRAVITATIONAL FORCE ===
    
    % Separation vector (from bennu to itokawa)
    r12 = itokawa.r - bennu.r;
    d = norm(r12);
    
    % Unit vector
    if d > 0
        n12 = r12 / d;
        
        % Gravitational force magnitude: F = G * m1 * m2 / d²
        F_mag = G * bennu.mass * itokawa.mass / (d^2);
        
        % Force on bennu (toward itokawa)
        F_bennu = F_mag * n12;
        
        % Force on itokawa (toward bennu) - Newton's 3rd law
        F_itokawa = -F_bennu;
    else
        % Avoid singularity if asteroids overlap
        F_bennu = [0; 0; 0];
        F_itokawa = [0; 0; 0];
    end
    
    % Store forces
    trajectory.bennu.F(k,:) = F_bennu';
    trajectory.itokawa.F(k,:) = F_itokawa';
    
    % Store separation distance
    trajectory.separation(k) = d;
    
    % === STEP 5.1: DETECT CONTACT (BROAD PHASE) ===
    % Check if bounding spheres overlap
    if d <= (bennu.Rbound + itokawa.Rbound)
        trajectory.contact(k) = true;  % Possible contact detected
        
        % === STEP 5.2: DEFINE CONTACT GEOMETRY ===
        % Contact normal (from bennu → itokawa)
        contact_normal = n12;  % Already computed: r12 / d
        
        % Contact point - use actual mesh geometry for realistic torque
        % Find contact vertex on Bennu's surface (farthest in contact direction)
        R_bennu_current = quat2rotm_custom(bennu.q);
        vertices_world_bennu = (R_bennu_current * bennu.vertices')' + bennu.r';
        projections_bennu = vertices_world_bennu * contact_normal;
        [~, idx_bennu] = max(projections_bennu);
        contact_point_bennu = vertices_world_bennu(idx_bennu, :)';
        
        % Find contact vertex on Itokawa's surface (farthest opposite direction)
        R_itokawa_current = quat2rotm_custom(itokawa.q);
        vertices_world_itokawa = (R_itokawa_current * itokawa.vertices')' + itokawa.r';
        projections_itokawa = vertices_world_itokawa * (-contact_normal);
        [~, idx_itokawa] = max(projections_itokawa);
        contact_point_itokawa = vertices_world_itokawa(idx_itokawa, :)';
        
        % Contact point is midpoint between the two contact vertices
        contact_point = 0.5 * (contact_point_bennu + contact_point_itokawa);
        
        % Lever arms relative to COM
        r1c = contact_point - bennu.r;      % bennu COM to contact point
        r2c = contact_point - itokawa.r;    % itokawa COM to contact point
        
        % Store contact geometry
        trajectory.contact_normal(k,:) = contact_normal';
        trajectory.contact_point(k,:) = contact_point';
        trajectory.bennu_lever(k,:) = r1c';
        trajectory.itokawa_lever(k,:) = r2c';
        
        % === STEP 5.3: RELATIVE VELOCITY AT CONTACT POINT ===
        % Rigid-body velocity at contact point: v_p = v + ω × r
        vp1 = bennu.v + cross(bennu.omega, r1c);
        vp2 = itokawa.v + cross(itokawa.omega, r2c);
        
        % Relative velocity at contact point
        vrel = vp2 - vp1;
        
        % === STEP 6.1: DECOMPOSE RELATIVE VELOCITY ===
        % Normal component (positive = separating, negative = closing)
        vn = dot(vrel, contact_normal);
        
        % Tangential component (sliding velocity)
        vt = vrel - vn * contact_normal;
        
        % === STEP 6.2: COMPUTE TANGENTIAL DIRECTION ===
        % Unit tangent direction (friction direction)
        if norm(vt) > 1e-8
            t = vt / norm(vt);   % unit tangent direction
        else
            t = [0; 0; 0];       % no sliding (pure normal collision)
        end
        
        % Store relative velocity data
        trajectory.vrel(k,:) = vrel';
        trajectory.vn(k) = vn;
        trajectory.vt(k,:) = vt';
        
        % === STEP 6.3: CANDIDATE TANGENTIAL IMPULSE (STICKING) ===
        % Compute tangential impulse that would fully cancel sliding
        % This is the "sticking" impulse (before friction limit)
        Jt_star = 0;  % Default: no tangential impulse
        
        if norm(vt) > 1e-8
            % Transform inertia tensors to world frame (if not already done)
            R_bennu = quat2rotm_custom(bennu.q);
            R_itokawa = quat2rotm_custom(itokawa.q);
            I_world_bennu = R_bennu * bennu.I_body * R_bennu';
            I_world_itokawa = R_itokawa * itokawa.I_body * R_itokawa';
            
            % Compute denominator for tangential impulse
            % denominator = 1/m₁ + 1/m₂ + t·[(I₁⁻¹(r₁×t))×r₁ + (I₂⁻¹(r₂×t))×r₂]
            r1_cross_t = cross(r1c, t);
            r2_cross_t = cross(r2c, t);
            
            term1_t = cross(I_world_bennu \ r1_cross_t, r1c);
            term2_t = cross(I_world_itokawa \ r2_cross_t, r2c);
            
            denominator_t = (1/bennu.mass) + (1/itokawa.mass) + dot(t, term1_t + term2_t);
            
            % Tangential impulse magnitude to cancel all sliding
            % Jt_star = -|vt| / denominator (negative because we want to stop the motion)
            Jt_star = -norm(vt) / denominator_t;
        end
        
        % === STEP 6.6: RESTING CONTACT DETECTION ===
        % Check if this is a resting contact (low relative velocity)
        % For resting contacts, use contact dynamics instead of collision response
        is_resting = abs(vn) < v_resting_threshold;
        
        % Check if bodies are closing OR in resting contact with penetration
        if vn < 0 || (is_resting && vn < v_resting_threshold)
            trajectory.closing(k) = true;  % Bodies are approaching or in contact
            
            % === STEP 5.4: COMPUTE NORMAL IMPULSE ===
            % Transform inertia tensors to world frame
            R_bennu = quat2rotm_custom(bennu.q);
            R_itokawa = quat2rotm_custom(itokawa.q);
            I_world_bennu = R_bennu * bennu.I_body * R_bennu';
            I_world_itokawa = R_itokawa * itokawa.I_body * R_itokawa';
            
            % Compute denominator terms
            % n · [(I1^-1(r1 × n)) × r1 + (I2^-1(r2 × n)) × r2]
            r1_cross_n = cross(r1c, contact_normal);
            r2_cross_n = cross(r2c, contact_normal);
            
            term1 = cross(I_world_bennu \ r1_cross_n, r1c);
            term2 = cross(I_world_itokawa \ r2_cross_n, r2c);
            
            denominator = (1/bennu.mass) + (1/itokawa.mass) + dot(contact_normal, term1 + term2);
            
            % Select coefficient of restitution based on contact type
            if is_resting
                % Resting contact: use e=0 (perfectly inelastic) to prevent bouncing
                % Apply impulse to maintain non-penetration constraint: vn ≥ 0
                e_contact = 0.0;
            else
                % Impact: use normal coefficient of restitution
                e_contact = e_restitution;
            end
            
            % Normal impulse magnitude: J_n = -(1 + e) * v_n / denominator
            Jn = -(1 + e_contact) * vn / denominator;
            
            % Store impulse magnitude
            trajectory.Jn(k) = Jn;
            
            % === STEP 6.4: APPLY COULOMB FRICTION LIMIT ===
            % Real friction is bounded by: |Jt| ≤ μ * Jn
            % This creates two regimes: sticking vs sliding
            
            % Initialize tangential impulse
            Jt = 0;  % Default: no tangential impulse
            
            if norm(vt) > 1e-8  % Only apply friction if there's tangential motion
                % Friction limit (maximum tangential impulse magnitude)
                Jt_max = mu_friction * Jn;
                
                % Check sticking vs sliding
                if abs(Jt_star) <= Jt_max
                    % Case 1: Sticking contact
                    % Friction is strong enough to cancel tangential motion
                    Jt = Jt_star;  % Apply full candidate impulse
                    trajectory.sticking(k) = true;
                else
                    % Case 2: Sliding contact
                    % Friction saturates at Coulomb limit
                    % Apply maximum friction in direction opposing motion
                    Jt = sign(Jt_star) * Jt_max;  % Clamp to friction limit (preserve sign!)
                    trajectory.sticking(k) = false;
                end
                
                % Store tangential impulse magnitude
                trajectory.Jt(k) = Jt;
            end
            
            % ===== VALIDATION DIAGNOSTICS =====
            % Capture pre-impulse state for validation
            v1_pre = bennu.v;
            v2_pre = itokawa.v;
            omega1_pre = bennu.omega;
            omega2_pre = itokawa.omega;
            
            % Pre-impulse kinetic energy (for energy dissipation check)
            KE_pre = 0.5 * bennu.mass * (v1_pre' * v1_pre) + ...
                     0.5 * itokawa.mass * (v2_pre' * v2_pre) + ...
                     0.5 * omega1_pre' * bennu.I_body * omega1_pre + ...
                     0.5 * omega2_pre' * itokawa.I_body * omega2_pre;
            
            % Pre-impulse linear momentum (should be conserved)
            p_pre = bennu.mass * v1_pre + itokawa.mass * v2_pre;
            
            % === STEP 6.5: APPLY COMBINED IMPULSE (NORMAL + TANGENTIAL) ===
            % Total impulse vector at contact point:
            % J = Jn*n + Jt*t (normal component + tangential component)
            
            % Compute tangential impulse vector
            if norm(vt) > 1e-8
                Jt_vec = Jt * t;  % Tangential impulse in direction t
            else
                Jt_vec = [0; 0; 0];  % No tangential component
            end
            
            % Normal impulse vector
            Jn_vec = Jn * contact_normal;
            
            % Combined impulse vector (world frame)
            J_total = Jn_vec + Jt_vec;
            
            % Linear velocity updates
            % Newton's 3rd law: impulse on bennu is -J, on itokawa is +J
            bennu.v = bennu.v - J_total / bennu.mass;
            itokawa.v = itokawa.v + J_total / itokawa.mass;
            
            % Angular velocity updates: ΔL = r × J
            % Compute angular impulse in world frame
            dL_world_bennu = cross(r1c, J_total);
            dL_world_itokawa = cross(r2c, J_total);
            
            % Transform angular impulse to body frame (omega is in body frame!)
            dL_body_bennu = R_bennu' * dL_world_bennu;
            dL_body_itokawa = R_itokawa' * dL_world_itokawa;
            
            % Update angular velocities: Δω = I_body^(-1) * ΔL_body
            bennu.omega = bennu.omega - (bennu.I_body \ dL_body_bennu);
            itokawa.omega = itokawa.omega + (itokawa.I_body \ dL_body_itokawa);
            
            % ===== VALIDATION CHECKS =====
            % Post-impulse state
            v1_post = bennu.v;
            v2_post = itokawa.v;
            omega1_post = bennu.omega;
            omega2_post = itokawa.omega;
            
            % Post-impulse kinetic energy
            KE_post = 0.5 * bennu.mass * (v1_post' * v1_post) + ...
                      0.5 * itokawa.mass * (v2_post' * v2_post) + ...
                      0.5 * omega1_post' * bennu.I_body * omega1_post + ...
                      0.5 * omega2_post' * itokawa.I_body * omega2_post;
            
            % Post-impulse linear momentum
            p_post = bennu.mass * v1_post + itokawa.mass * v2_post;
            
            % Validation outputs (first collision only)
            if trajectory.t(k) == trajectory.t(find(trajectory.closing, 1, 'first'))
                fprintf('\n===== COLLISION VALIDATION (t = %.2f s) =====\n', trajectory.t(k));
                fprintf('✔ VELOCITY CHANGES:\n');
                fprintf('  Bennu Δv:      [%.4f, %.4f, %.4f] m/s\n', v1_post - v1_pre);
                fprintf('  Itokawa Δv:    [%.4f, %.4f, %.4f] m/s\n', v2_post - v2_pre);
                fprintf('  Bennu Δω:      [%.4f, %.4f, %.4f] rad/s\n', omega1_post - omega1_pre);
                fprintf('  Itokawa Δω:    [%.4f, %.4f, %.4f] rad/s\n', omega2_post - omega2_pre);
                fprintf('  |Δω_bennu|:    %.4f rad/s (%.2f deg/s)\n', ...
                    norm(omega1_post - omega1_pre), norm(omega1_post - omega1_pre)*180/pi);
                fprintf('  |Δω_itokawa|:  %.4f rad/s (%.2f deg/s)\n', ...
                    norm(omega2_post - omega2_pre), norm(omega2_post - omega2_pre)*180/pi);
                
                fprintf('\n✔ MOMENTUM CONSERVATION:\n');
                fprintf('  p_before:      [%.6e, %.6e, %.6e] kg·m/s\n', p_pre);
                fprintf('  p_after:       [%.6e, %.6e, %.6e] kg·m/s\n', p_post);
                fprintf('  Δp:            [%.6e, %.6e, %.6e] kg·m/s\n', p_post - p_pre);
                fprintf('  |Δp|/|p|:      %.6e (should be ~0)\n', norm(p_post - p_pre)/norm(p_pre));
                
                fprintf('\n✔ ENERGY DISSIPATION:\n');
                fprintf('  KE_before:     %.6e J\n', KE_pre);
                fprintf('  KE_after:      %.6e J\n', KE_post);
                fprintf('  ΔKE:           %.6e J (should be negative)\n', KE_post - KE_pre);
                fprintf('  Energy ratio:  %.4f (should be < 1 for e=%.2f)\n', ...
                    KE_post/KE_pre, e_restitution);
                
                fprintf('\n✔ CONTACT GEOMETRY:\n');
                fprintf('  Contact point: [%.2f, %.2f, %.2f] m\n', contact_point);
                fprintf('  Lever arm 1:   [%.2f, %.2f, %.2f] m (length: %.2f m)\n', ...
                    r1c, norm(r1c));
                fprintf('  Lever arm 2:   [%.2f, %.2f, %.2f] m (length: %.2f m)\n', ...
                    r2c, norm(r2c));
                fprintf('  Impulse |J|:   %.4f N·s\n', Jn);
                fprintf('  Closing speed: %.4f m/s\n', -vn);
                
                if KE_post > KE_pre
                    fprintf('\n⚠ WARNING: Energy increased! Check impulse sign.\n');
                end
                if norm(p_post - p_pre)/norm(p_pre) > 1e-6
                    fprintf('\n⚠ WARNING: Momentum not conserved! Check impulse application.\n');
                end
                if norm(omega1_post - omega1_pre) < 1e-6 && norm(r1c) > 1e-3
                    fprintf('\n⚠ WARNING: No spin change despite off-center contact!\n');
                end
            end
        end
    end
    
    % === STEP 5.6: RESUME NORMAL INTEGRATION ===
    % After impulse application, continue with Step 4 (gravity + integration)
    % The velocities have been updated instantaneously
    % Now proceed with normal gravitational forces and time integration
    % Collisions may happen again in subsequent timesteps
    
    % === COMPUTE TOTAL ENERGY AND ANGULAR MOMENTUM ===
    
    % Kinetic energy
    KE_bennu = 0.5 * bennu.mass * (bennu.v' * bennu.v);
    KE_itokawa = 0.5 * itokawa.mass * (itokawa.v' * itokawa.v);
    
    % Rotational energy
    RE_bennu = 0.5 * bennu.omega' * bennu.I_body * bennu.omega;
    RE_itokawa = 0.5 * itokawa.omega' * itokawa.I_body * itokawa.omega;
    
    % Gravitational potential energy
    PE = -G * bennu.mass * itokawa.mass / d;
    
    % Total energy
    trajectory.E_total(k) = KE_bennu + KE_itokawa + RE_bennu + RE_itokawa + PE;
    
    % Total angular momentum (orbital + spin)
    % Orbital angular momentum: L_orbital = r × (m*v)
    r_com = (bennu.mass * bennu.r + itokawa.mass * itokawa.r) / (bennu.mass + itokawa.mass);
    L_orbital_bennu = cross(bennu.r - r_com, bennu.mass * bennu.v);
    L_orbital_itokawa = cross(itokawa.r - r_com, itokawa.mass * itokawa.v);
    
    % Spin angular momentum (transform from body to world frame)
    % L_world = R * (I_body * omega_body)
    R_bennu = quat2rotm_custom(bennu.q);
    R_itokawa = quat2rotm_custom(itokawa.q);
    L_spin_bennu = R_bennu * (bennu.I_body * bennu.omega);
    L_spin_itokawa = R_itokawa * (itokawa.I_body * itokawa.omega);
    
    trajectory.L_total(k,:) = (L_orbital_bennu + L_orbital_itokawa + L_spin_bennu + L_spin_itokawa)';
    
    % === INTEGRATION STEP ===
    
    % --- BENNU ---
    % STEP 4.2: Update translational motion (semi-implicit Euler)
    % v̇ = F/m (Newton's 2nd law)
    bennu.v = bennu.v + (F_bennu / bennu.mass) * dt;  % Update velocity first
    bennu.r = bennu.r + bennu.v * dt;                 % Then position with new velocity
    
    % STEP 4.3: Update rotation (unchanged from Step 3 - no gravitational torque at COM)
    % Euler's equations in body frame: I·ω̇ = -ω × (I·ω)
    omega_dot_bennu = -bennu.I_body \ cross(bennu.omega, bennu.I_body * bennu.omega);
    bennu.omega = bennu.omega + omega_dot_bennu * dt;
    
    % Quaternion propagation: q̇ = 0.5·Ω(ω)·q
    wx = bennu.omega(1); wy = bennu.omega(2); wz = bennu.omega(3);
    Omega_bennu = [ 0   -wx  -wy  -wz;
                    wx   0    wz  -wy;
                    wy  -wz   0    wx;
                    wz   wy  -wx   0 ];
    q_dot_bennu = 0.5 * Omega_bennu * bennu.q;
    bennu.q = bennu.q + q_dot_bennu * dt;
    bennu.q = bennu.q / norm(bennu.q);  % Normalize
    
    % --- ITOKAWA ---
    % STEP 4.2: Update translational motion (semi-implicit Euler)
    % v̇ = F/m (Newton's 2nd law)
    itokawa.v = itokawa.v + (F_itokawa / itokawa.mass) * dt;  % Update velocity first
    itokawa.r = itokawa.r + itokawa.v * dt;                   % Then position with new velocity
    
    % STEP 4.3: Update rotation (unchanged from Step 3 - no gravitational torque at COM)
    % Euler's equations in body frame: I·ω̇ = -ω × (I·ω)
    omega_dot_itokawa = -itokawa.I_body \ cross(itokawa.omega, itokawa.I_body * itokawa.omega);
    itokawa.omega = itokawa.omega + omega_dot_itokawa * dt;
    
    % Quaternion propagation: q̇ = 0.5·Ω(ω)·q
    wx = itokawa.omega(1); wy = itokawa.omega(2); wz = itokawa.omega(3);
    Omega_itokawa = [ 0   -wx  -wy  -wz;
                      wx   0    wz  -wy;
                      wy  -wz   0    wx;
                      wz   wy  -wx   0 ];
    q_dot_itokawa = 0.5 * Omega_itokawa * itokawa.q;
    itokawa.q = itokawa.q + q_dot_itokawa * dt;
    itokawa.q = itokawa.q / norm(itokawa.q);  % Normalize
    
    % Progress indicator
    if mod(k, 200) == 0
        fprintf('  Step %d/%d: separation = %.1f m, E_total = %.4e J\n', ...
            k, N_steps, d, trajectory.E_total(k));
    end
end

fprintf('Integration completed.\n');

%% Final state
fprintf('\nFinal state:\n');
fprintf('Bennu:\n');
fprintf('  Position:   [%.2f, %.2f, %.2f] m\n', bennu.r);
fprintf('  Velocity:   [%.3f, %.3f, %.3f] m/s\n', bennu.v);
fprintf('  Spin rate:  %.4f rad/s\n', norm(bennu.omega));

fprintf('Itokawa:\n');
fprintf('  Position:   [%.2f, %.2f, %.2f] m\n', itokawa.r);
fprintf('  Velocity:   [%.3f, %.3f, %.3f] m/s\n', itokawa.v);
fprintf('  Spin rate:  %.4f rad/s\n', norm(itokawa.omega));

fprintf('\nFinal separation: %.1f m\n', trajectory.separation(end));
fprintf('Minimum separation: %.1f m (at t = %.1f s)\n', ...
    min(trajectory.separation), trajectory.t(find(trajectory.separation == min(trajectory.separation), 1)));

% Contact detection summary
num_contacts = sum(trajectory.contact);
num_closing = sum(trajectory.closing);
fprintf('\nCollision Detection (Broad Phase):\n');
fprintf('  Bennu bounding radius:    %.1f m\n', bennu.Rbound);
fprintf('  Itokawa bounding radius:  %.1f m\n', itokawa.Rbound);
fprintf('  Contact threshold:        %.1f m\n', bennu.Rbound + itokawa.Rbound);
fprintf('  Contacts detected:        %d / %d timesteps (%.1f%%)\n', ...
    num_contacts, N_steps, 100*num_contacts/N_steps);
if num_contacts > 0
    first_contact_idx = find(trajectory.contact, 1, 'first');
    fprintf('  First contact at:         t = %.1f s (step %d)\n', ...
        trajectory.t(first_contact_idx), first_contact_idx);
    fprintf('  Closing contacts:         %d / %d (%.1f%%)\n', ...
        num_closing, num_contacts, 100*num_closing/num_contacts);
    if num_closing > 0
        closing_indices = find(trajectory.closing);
        max_closing_speed = -min(trajectory.vn(closing_indices));
        max_impulse = max(trajectory.Jn(closing_indices));
        total_impulse = sum(trajectory.Jn(closing_indices));
        fprintf('  Max closing speed:        %.4f m/s\n', max_closing_speed);
        fprintf('  Max impulse magnitude:    %.4f N·s\n', max_impulse);
        fprintf('  Total impulse delivered:  %.4f N·s\n', total_impulse);
        fprintf('  Collisions resolved:      %d\n', num_closing);
        
        % Friction statistics
        num_with_friction = sum(abs(trajectory.Jt(closing_indices)) > 0);
        num_sticking = sum(trajectory.sticking(closing_indices));
        num_sliding = num_with_friction - num_sticking;
        
        fprintf('\nFriction Statistics:\n');
        fprintf('  Contacts with friction:   %d / %d (%.1f%%)\n', ...
            num_with_friction, num_closing, 100*num_with_friction/num_closing);
        if num_with_friction > 0
            fprintf('  Sticking contacts:        %d (%.1f%%)\n', ...
                num_sticking, 100*num_sticking/num_with_friction);
            fprintf('  Sliding contacts:         %d (%.1f%%)\n', ...
                num_sliding, 100*num_sliding/num_with_friction);
            fprintf('  Max tangential impulse:   %.4f N·s\n', ...
                max(abs(trajectory.Jt(closing_indices))));
            fprintf('  Total tangential impulse: %.4f N·s\n', ...
                sum(abs(trajectory.Jt(closing_indices))));
        end
        
        % Resting contact detection
        vn_at_closing = trajectory.vn(closing_indices);
        num_resting = sum(abs(vn_at_closing) < v_resting_threshold);
        if num_resting > 0
            fprintf('\nResting Contact Statistics:\n');
            fprintf('  Resting contacts:         %d / %d (%.1f%%)\n', ...
                num_resting, num_closing, 100*num_resting/num_closing);
        end
    end
else
    fprintf('  No contacts detected during simulation.\n');
end

%% Conservation verification
E_initial = trajectory.E_total(1);
E_final = trajectory.E_total(end);
L_initial = trajectory.L_total(1,:);
L_final = trajectory.L_total(end,:);

fprintf('\nConservation verification:\n');
fprintf('  Energy error:           %.6e J (%.3f%%)\n', ...
    E_final - E_initial, 100*abs(E_final - E_initial)/abs(E_initial));
fprintf('  Angular momentum error: %.6e kg·m²/s (%.3f%%)\n', ...
    norm(L_final - L_initial), 100*norm(L_final - L_initial)/norm(L_initial));

%% Visualization
fprintf('\nCreating visualizations...\n');

% Figure 1: Trajectories
figure('Name', 'Two-Body Trajectories', 'Position', [50, 50, 1400, 900]);

% 3D trajectory plot
subplot(2,2,1);
plot3(trajectory.bennu.r(:,1), trajectory.bennu.r(:,2), trajectory.bennu.r(:,3), ...
      'r-', 'LineWidth', 2, 'DisplayName', 'Bennu');
hold on;
plot3(trajectory.itokawa.r(:,1), trajectory.itokawa.r(:,2), trajectory.itokawa.r(:,3), ...
      'b-', 'LineWidth', 2, 'DisplayName', 'Itokawa');
plot3(trajectory.bennu.r(1,1), trajectory.bennu.r(1,2), trajectory.bennu.r(1,3), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot3(trajectory.itokawa.r(1,1), trajectory.itokawa.r(1,2), trajectory.itokawa.r(1,3), ...
      'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('3D Trajectories');
legend('Location', 'best');
axis equal;
view(45, 30);

% Separation vs time
subplot(2,2,2);
plot(trajectory.t, trajectory.separation, 'k-', 'LineWidth', 2);
hold on;
yline(bennu.Rbound + itokawa.Rbound, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Contact Threshold');
% Highlight contact regions
if any(trajectory.contact)
    contact_times = trajectory.t(trajectory.contact);
    contact_seps = trajectory.separation(trajectory.contact);
    scatter(contact_times, contact_seps, 50, 'r', 'filled', 'DisplayName', 'Contact Detected');
end
grid on;
xlabel('Time (s)');
ylabel('Separation (m)');
title('Distance Between Asteroids (Broad Phase)');
legend('Location', 'best');

% Energy conservation
subplot(2,2,3);
% Compute spin rates (magnitude of angular velocity)
spin_bennu = vecnorm(trajectory.bennu.omega, 2, 2);
spin_itokawa = vecnorm(trajectory.itokawa.omega, 2, 2);

plot(trajectory.t, spin_bennu * 180/pi, 'r-', 'LineWidth', 2, 'DisplayName', 'Bennu');
hold on;
plot(trajectory.t, spin_itokawa * 180/pi, 'b-', 'LineWidth', 2, 'DisplayName', 'Itokawa');

% Mark collision events
if any(trajectory.closing)
    collision_times = trajectory.t(trajectory.closing);
    for i = 1:length(collision_times)
        xline(collision_times(i), 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

grid on;
xlabel('Time (s)');
ylabel('Spin Rate (deg/s)');
title('Angular Velocity (Friction Effects)');
legend('Location', 'best');

% Angular momentum conservation
subplot(2,2,4);
L_error = vecnorm(trajectory.L_total - L_initial, 2, 2);
plot(trajectory.t, L_error / norm(L_initial) * 100, 'm-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Angular Momentum Error (%)');
title('Total Angular Momentum Conservation');

sgtitle('Two-Body Simulation: Bennu and Itokawa with Mutual Gravity + Friction', 'FontSize', 14, 'FontWeight', 'bold');

% Figure 2: Force analysis
figure('Name', 'Gravitational Forces', 'Position', [100, 100, 1200, 600]);

subplot(1,2,1);
plot(trajectory.t, vecnorm(trajectory.bennu.F, 2, 2), 'r-', 'LineWidth', 2);
hold on;
plot(trajectory.t, vecnorm(trajectory.itokawa.F, 2, 2), 'b--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Force Magnitude (N)');
title('Gravitational Force Magnitude');
legend('Force on Bennu', 'Force on Itokawa');

subplot(1,2,2);
semilogy(trajectory.t, vecnorm(trajectory.bennu.F, 2, 2), 'r-', 'LineWidth', 2);
hold on;
semilogy(trajectory.t, vecnorm(trajectory.itokawa.F, 2, 2), 'b--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Force Magnitude (N)');
title('Gravitational Force Magnitude (log scale)');
legend('Force on Bennu', 'Force on Itokawa');

% Figure 3: Energy dissipation and friction effects
figure('Name', 'Energy & Friction Analysis', 'Position', [150, 150, 1200, 600]);

subplot(1,2,1);
% Energy conservation
plot(trajectory.t, (trajectory.E_total - E_initial)/abs(E_initial) * 100, 'g-', 'LineWidth', 2);
hold on;

% Mark collision events with energy loss
if any(trajectory.closing)
    collision_times = trajectory.t(trajectory.closing);
    for i = 1:length(collision_times)
        xline(collision_times(i), 'r:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

grid on;
xlabel('Time (s)');
ylabel('Energy Error (%)');
title('Total Energy (Dissipation from Collisions + Friction)');
legend('Energy', 'Location', 'best');

subplot(1,2,2);
% Impulse magnitudes over time
if any(trajectory.closing)
    closing_indices = find(trajectory.closing);
    Jn_closing = trajectory.Jn(closing_indices);
    Jt_closing = abs(trajectory.Jt(closing_indices));
    t_closing = trajectory.t(closing_indices);
    
    stem(t_closing, Jn_closing, 'r', 'LineWidth', 1.5, 'DisplayName', 'Normal Impulse');
    hold on;
    stem(t_closing, Jt_closing, 'b', 'LineWidth', 1.5, 'DisplayName', 'Tangential Impulse (Friction)');
    
    grid on;
    xlabel('Time (s)');
    ylabel('Impulse Magnitude (N·s)');
    title('Collision Impulses (Normal vs Tangential)');
    legend('Location', 'best');
end

fprintf('Visualization complete.\n');

%% Helper functions

function asteroid = loadBennu()
    % Load Bennu shape model and compute properties
    filename = 'data/101955bennu.tab';
    
    % Initialize arrays
    vertices = zeros(2000, 3);
    faces = zeros(3000, 3);
    
    % Open and parse file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open Bennu file: %s', filename);
    end
    
    vertexCount = 0;
    faceCount = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && length(line) > 1
            tokens = regexp(strtrim(line), '\s+', 'split');
            if length(tokens) >= 4
                if strcmp(tokens{1}, 'v')
                    vertexCount = vertexCount + 1;
                    if vertexCount > size(vertices, 1)
                        vertices = [vertices; zeros(1000, 3)];
                    end
                    vertices(vertexCount, :) = [str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4})];
                elseif strcmp(tokens{1}, 'f')
                    faceCount = faceCount + 1;
                    if faceCount > size(faces, 1)
                        faces = [faces; zeros(1000, 3)];
                    end
                    faces(faceCount, :) = [str2double(tokens{2}), str2double(tokens{3}), str2double(tokens{4})];
                end
            end
        end
    end
    fclose(fid);
    
    % Trim arrays
    vertices = vertices(1:vertexCount, :);
    faces = faces(1:faceCount, :);
    
    % Convert to meters
    vertices = vertices * 1000;
    
    % Shift to COM
    vertices = vertices - mean(vertices, 1);
    
    % Compute properties
    asteroid = struct();
    asteroid.name = 'Bennu';
    asteroid.vertices = vertices;
    asteroid.faces = faces;
    asteroid.Rbound = max(vecnorm(vertices, 2, 2));
    
    % Volume and mass
    asteroid.volume = calculateMeshVolume(vertices, faces);
    asteroid.density = 2.2;  % g/cm³
    asteroid.mass = asteroid.volume * asteroid.density * 1e3;  % kg
    
    % Inertia tensor
    [asteroid.I_body, asteroid.I_principal, asteroid.principal_axes] = ...
        calculateInertiaTensor(vertices, faces, asteroid.density);
end

function asteroid = loadItokawa()
    % Load Itokawa shape model and compute properties
    filename = 'vertex/ver512q.tab';
    
    % Open and parse file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open Itokawa file: %s', filename);
    end
    
    % Read header
    header_line = fgetl(fid);
    header_values = sscanf(header_line, '%d');
    num_vertices = header_values(1);
    num_faces = header_values(2);
    
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
    
    % Convert to meters
    vertices = vertices * 1000;
    
    % Shift to COM
    vertices = vertices - mean(vertices, 1);
    
    % Compute properties
    asteroid = struct();
    asteroid.name = 'Itokawa';
    asteroid.vertices = vertices;
    asteroid.faces = faces;
    asteroid.Rbound = max(vecnorm(vertices, 2, 2));
    
    % Volume and mass
    asteroid.volume = calculateMeshVolume(vertices, faces);
    asteroid.density = 1.9;  % g/cm³
    asteroid.mass = asteroid.volume * asteroid.density * 1e3;  % kg
    
    % Inertia tensor
    [asteroid.I_body, asteroid.I_principal, asteroid.principal_axes] = ...
        calculateInertiaTensor(vertices, faces, asteroid.density);
end

function volume = calculateMeshVolume(vertices, faces)
    volume = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i,1), :);
        v2 = vertices(faces(i,2), :);
        v3 = vertices(faces(i,3), :);
        volume = volume + v1(1)*(v2(2)*v3(3) - v2(3)*v3(2)) + ...
                          v2(1)*(v3(2)*v1(3) - v3(3)*v1(2)) + ...
                          v3(1)*(v1(2)*v2(3) - v1(3)*v2(2));
    end
    volume = abs(volume) / 6;
end

function [I_tensor, I_principal, principal_axes] = calculateInertiaTensor(vertices, faces, density_gcm3)
    density = density_gcm3 * 1000;
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
    
    [principal_axes, I_diag] = eig(I_tensor);
    I_principal = diag(I_diag);
    [I_principal, sort_idx] = sort(I_principal);
    principal_axes = principal_axes(:, sort_idx);
    
    if det(principal_axes) < 0
        principal_axes(:,3) = -principal_axes(:,3);
    end
end

function R = quat2rotm_custom(q)
    if size(q, 2) == 4
        q = q';
    end
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1-2*(y^2+z^2),   2*(x*y-w*z),     2*(x*z+w*y);
         2*(x*y+w*z),     1-2*(x^2+z^2),   2*(y*z-w*x);
         2*(x*z-w*y),     2*(y*z+w*x),     1-2*(x^2+y^2)];
end
