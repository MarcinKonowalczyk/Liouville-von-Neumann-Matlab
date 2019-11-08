% This script performs a representative Liouville von Neumann simulation by
% propagation of the density matrix. The quantum mechanical system used for
% simulation consists of three spins: electrons (A,B) and a nucleus (C).
% Only one of the electrons is coupled to the nucleus (A-C) with a
% hyperfine coupling specified by `hfc`. The system is also subject to an
% external magnetic fields specified by `B0`. The calculation runs for the
% time points specified by `T`.
%
% The code is intended to be used to learn about the basics of spin
% chemistry, not as a tool for simulation. It is heavily commented, 
% and to use it you should go though it line-by-line to understand what it
% does.
%
% Below we explain different modes for the observable function and the
% initial spin density which we've deemed interesting and have therefore
% programmed it.
% 
% OBSERVABLE
% Option which specifies the function
%  singlet - Singlet probability
%  triplet - Triplet probability
%  espin   - Total spin (squared) of the electrons
%  total   - Total spin (squared) of all the nuclei
%  trace   - The trace of the density matrix
%
% RH0
% Option which specifies the initial spin density matrix
%  singlet   - 1/2 in SingletAlpha and 12 in SingletBeta
%  triplet   - 1/6th in each of the triplet sublevels
%  inactive  - 1/2 in T+Alpha (3x Alpha) and 1/2 in T-Beta (3x Beta)
%  triplet+  - 1/2 in Triplet+Alpha and 1/2 in Triplet+Beta
%  triplet0  - 1/2 in Triplet0Alpha and 1/2 in Triplet0Beta
%  singleta  - 100% in SingletAlpha
%  coherence - 1/2 in the coherence between SingletAlpha and Triplet0Beta,
%              and 1/2 in the coherence between SingletBeta and
%              Triplet0Alpha. This is a 90 degree clockwise rotation of
%              the singlet case, and it evolves in an equivalent fasion
%              through the coherences.
% 
% Written by Marcin Konowalczyk and Gabriel Moise
% Timmel Group @ Oxford University; 2017
 
close all; clear; clc;
 
%% Things you can change
% ==========================================
% Simulation options
hfc = 1.0; % Hyperfine coupling in mT
B0 = 0.05; % Magnetic field in mT
T = linspace(0,3e-6,10000); % Time vector
opt.observable = 's0'; % 'singlet', 'triplet', 'espin', 'total', 'trace'
opt.rh0 = 'singlet'; % 'singlet', 'triplet', 'inactive', 'triplet+', 'triplet0', 'singleta', 'coherence'
            
% Calculation options
opt.fixtrace = true; % Fix trace of the density matrix to 1 after each timestep
opt.check = false; % Do sanity checks (for advanced users only)
 
% Visualisation options
opt.verbosity = false; % Write to the console
opt.rhoPlot = false; % Ctrl + C to stop the execution
opt.rhoPlotSkip = 10; % Do phoPlot every n frames
% ==========================================
 
%% Basic spin matrices (in alpha/beta basis)
Ix = 0.5 * [0  1; 1   0];
Iy = 0.5 * [0 -1i;1i  0];
Iz = 0.5 * [1  0; 0  -1];
Ie = eye(2);
 
tkron = @(x,y,z) kron(kron(x,y),z); % Successive kronecker product of three elements
tol = 10*eps; % Tolerance for asserts
 
%% Transform to Singlet Triplet basis
s2 = 1/sqrt(2);
%      T+ S  T0   T-
Pe = [ 1  0   0   0
    0 s2  s2   0
    0 -s2 s2   0
    0  0   0   1];
P = kron(Pe,eye(2)); % Nucleus stays in alpha/beta basis
ST = @(A) P*A*P'; % Transfomration to the ST basis
 
% Check that P is unitary
if opt.check
    condition = ST(eye(8)) - eye(8) < tol;
    assert(all(condition(:)),'P is not unitary');
end
 
%% Product spin matrices
% Note: In Matlab you cannot (easily) work on the vector operators.
% Therefore, we work with the components thereof.
 
% Components of SA vector operator
SAx = ST(tkron(Ix,Ie,Ie));
SAy = ST(tkron(Iy,Ie,Ie));
SAz = ST(tkron(Iz,Ie,Ie));
SAp = SAx + 1i*SAy; % Raising and lowering operators
SAm = SAx - 1i*SAy;
 
% Components of SB vector operator
SBx = ST(tkron(Ie,Ix,Ie));
SBy = ST(tkron(Ie,Iy,Ie));
SBz = ST(tkron(Ie,Iz,Ie));
SBp = SBx + 1i*SBy;
SBm = SBx - 1i*SBy;
 
% Components of SC vector operator
SCx = ST(tkron(Ie,Ie,Ix));
SCy = ST(tkron(Ie,Ie,Iy));
SCz = ST(tkron(Ie,Ie,Iz));
 
% Total spin angular momentum (squared) of each species
SA2 = SAx^2 + SAy^2 + SAz^2;
SB2 = SBx^2 + SBy^2 + SBz^2;
SC2 = SCx^2 + SCy^2 + SCz^2;
 
% Product operators
SASB = SAx*SBx + SAy*SBy + SAz*SBz;
SASC = SAx*SCx + SAy*SCy + SAz*SCz;
SBSC = SBx*SCx + SBy*SCy + SBz*SCz;
 
% Total electron spin operator (squared)
SAB2 = SA2 + SB2 + 2*SASB; % Method #1
%SAB2 = SA2 + SB2 + 2*SAz*SBz + SAp*SBm + SBp*SAm; % Method #2
 
% Components of the SAB vector operator
SABx = SAx + SBx;
SABy = SAy + SBy;
SABz = SAz + SBz;
 
% Product of the vector operator SAB SC
SABSC = SABx*SCx + SABy*SCy + SABz*SCz;
 
% Total spin operator F2 = (SA + SB + SC)^2 <- SA, SB, SC are vector operators
%F2 = SA2 + SB2 + SC2 + 2*SASB + 2*SASC + 2*SBSC; % Method #1
F2 = SAB2 + SC2 + 2*SABSC; % Method 2
 
gn = 267.513e3; % Hydrogen nucleus gyromagnetic ratio in rad/s/mT
ge = 1.760859644e8; % Electron gyromagnetic ratio in rad/s/mT
 
%% Make the Hamiltonian
a = ge * hfc; % Hyperfine coupling constant
omega0 = ge*B0; % Electron Larmor frequency
omega0n = -gn*B0; % Nuclear Larmor frequency
 
H_mag = omega0.*(SAz + SBz) + omega0n.*SCz; % Electron and nuclear Zeeman
H_hyperfine = a*(SAx*SCx + SAy*SCy + SAz*SCz);
 
H = H_mag + H_hyperfine;
 
[V,E] = eig(H); % Get the eigenvalues and eigenvectors
E = diag(E);
 
%% Define the projection operators, propagators and the observable
% Singlet projection operator in the aba basis
PS = 0.25*eye(8) - tkron(Ix,Ix,Ie) - tkron(Iy,Iy,Ie) - tkron(Iz,Iz,Ie);
PT = eye(8) - PS;
PS = ST(PS); % Project into ST basis
PT = ST(PT);
 
%Create the initial density matrix
switch lower(opt.rh0)
    case {'singlet', 's', 's0'}
        % Initial density matrix in all singlet sublevels (Sa and Sb)
        %rh0 = PS/2; % Method 1
        rh0 = zeros(8); rh0(5,5) = 0.5; rh0(6,6) = 0.5; % Method 2
    case {'triplet', 't'}
        % Initial density matrix in all triplet sublevels
        rh0 = PT/6;
    case 'inactive'
        % Half of the population in T+a and half in T-b
        % (this should not evolve in time)
        rh0 = zeros(8); rh0(1,1) = 0.5; rh0(8,8) = 0.5;
    case {'triplet+', 't+'}
        % Initial density matrix in T+
        rh0 = zeros(8); rh0(1,1) = 0.5; rh0(2,2) = 0.5;
    case {'triplet0', 't0'}
        % Initial density matrix in T0
        rh0 = zeros(8); rh0(3,3) = 0.5; rh0(4,4) = 0.5;
    case {'singleta', 'sa'}
        % Initial density matrix in Sa
        rh0 = zeros(8); rh0(5,5) = 1;
    case {'coherence', 'cs', 'cs0'}
        % Initial density matrix in coherence-equivalent rho of S0
        rh0 = zeros(8); rh0(5,4) = 0.5; rh0(6,3) = 0.5;
    otherwise
        error('Invalid starting point specified'); 
end
 
% Sanity checks
if opt.check
    assert(abs(trace(rh0)-1) < tol,'Invalid starting point trace')
    condition = diag(rh0) >= -tol;
    assert(all(condition(:)),'Invalid sign of the starting point')
end
    
% Calcualte the propagator matrices
dt = mean(diff(T));
Um = expm(-1i*H*dt); Up = expm(1i*H*dt);
rhp = rh0; % Spin density at previous timestep
 
% Define the function to calculate the observable on the density matrix
switch lower(opt.observable)
    case {'singlet', 's', 's0'}
        % Singlet yield
        observable = @(rh) real(trace(PS*rh));
        %observable = @(rh) real(rh(5,5) + rh(6,6));
    case {'triplet', 't'}
        % Three different approaches to triplet yield
        observable = @(rh) 1-real(trace(PS*rh));
        %observable = @(rh) real(trace(PT*rh));
        %observable = @(rh) real(trace(rh(1,1)+rh(2,2)+rh(3,3)+rh(4,4)+rh(7,7)+rh(8,8)));
    case 'espin'
        % Total SA ans SB spin
        observable = @(rh) real(trace(SAB2*rh));
    case {'total', 'f'}
        % Total spin
        observable = @(rh) real(trace(F2*rh));
    case 'trace'
        % Trace of the density matrix
        observable = @(rh) (1-real(trace(rh)))/eps;
    otherwise
        error('Invalid observable specified');
end
 
%% Iterate through all the time points
sig = NaN(1,length(T)); % Prealocate the matrix for the signal
sig(1) = observable(rh0); % Calculate the observable of the initial density matrix
Tsig = (0:length(T)-1)*dt; % Timesteps at which the 'sig' was simulated
for i = 2:length(T);
    t = T(i);
    if opt.verbosity && ~mod(i,length(T)/10); fprintf('%3.f%%\n',i./length(T).*100); end;
    
    % Propagate the density matrix by dt
    rht = Um*rhp*Up;
    if opt.fixtrace && opt.check; rht = rht/trace(rht); end; % Fix the roundof erros on trace
    rhp = rht;
    
    % Plot intermediate density matrix and the observable
    if opt.rhoPlot && ~mod(i,opt.rhoPlotSkip);
        figure(1);
        % Plot density matrix
        subplot(2,1,1);
        imagesc(abs(rht)); zlim([0 1]);
        view(0,90); axis square; axis tight;
        ticklabels = {'T+a','T+b','T0a','T0b','Sa','Sb','T-a','T-b'};
        set(gca,'XTickLabel',ticklabels,'YTickLabel',ticklabels);
        title('\rho(t)');
        % Plot the observable
        subplot(2,1,2);
        plot(Tsig,sig); xlim([0 max(Tsig)]);
        grid on; title('Observable as a funciton of time');
        drawnow; % Update the plot
    end
    
    % Sanity checks
    if opt.check
        % The trace of the density matrix must be equal to one (within machine precision)
        % Error accummulates therefore multiply `tol` by `i`
        assert(abs(trace(rht) - 1) < i*tol,'Invalid trace %d', abs(trace(rht) - 1)/tol)
 
        % Diagonal elements of the density matrix must be positive (within machine precision)
        condition = diag(rht) >= -tol;
        assert(all(condition(:)),'Invalid sign')
    end
    
    % Calculate the observable at time t
    sig(i) = observable(rht);
end
 
% Interpolate from the simulation points to the desired time points
sig = interp1(Tsig,sig,T,'pchip'); % Cubic interpolate the desired T points
 
if opt.verbosity; fprintf('Done\n'); end;
 
%% Plot
% Close the rhoPlot figure
if opt.rhoPlot;
    f = figure(1);
    close(f);
end
 
% Make the final plot of the observable as a function of time
figure;
plot(1e6*T,sig);
grid on;
xlabel('time / us')
switch lower(opt.observable)
    case {'singlet', 's', 's0'}
        ylabel('Singlet fraction')
        ylim([0 1]);
    case {'triplet', 't'}
        ylabel('Triplet fraction')
        ylim([0 1]);
    case 'espin'
        ylabel('Total electron spin')
        ylim([0 2]);
    case 'total'
        ylabel('Total spin')
        ylim([0 2.75]);
    case 'trace'
        ylabel('1 - Trace of the density matrix / eps')
end
title('Liouville von Neumann simulation by density matrix propagation')