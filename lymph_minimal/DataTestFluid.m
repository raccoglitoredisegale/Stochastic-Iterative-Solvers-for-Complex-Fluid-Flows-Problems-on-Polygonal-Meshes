%% Unsteady Stokes problem in [0 1]^2 with mixed conditions

Data.name = 'DataTestFluid';

Data.TagElFluid   = 1;         % Element tag
Data.TagBcFluid   = [2 3 4 5]; % Boundary tag
% The paper applies traction (Library 'D') on Bottom/Left [2/5] 
% and velocity (Library 'N') on Right/Top [3/4]
Data.LabBcFluid   = 'DNND';    % (D)irichlet (N)eumann 

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 75;       % number of elements for a new mesh
Data.MeshFromFile = false;     % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'UnitSquareMixed_20_el'; %filename for mesh 

%% Discretization properties                            
%% Time integration
Data.t0 = 0;
Data.T  = 0.25;      % Paper specifies T = 0.25 for this test
Data.dt = 1e-3;     % Paper specifies dt = 1e-3 for this test

Data.timeint   = 'CN';

%% Space discretization
Data.degree  = 1;        % Polynomial degree 
Data.penalty_coeff = 10; % Penalty coefficient

%% Visualization settings
Data.PlotExact   = true;
Data.PlotGridSol = true;
Data.VisualizationStep = 0.1;
Data.PlotIniCond  = false;
Data.ComputeVelAndPres = false;

%% Save Solution settings
Data.timesave  = round(Data.T / Data.dt);

%% properties fluid material
% Material parameter
Data.mu_f  = {@(x,y) 1 + 0.*x.*y};

% =========================================================================
% FORCING TERMS FOR THE CONTINUOUS EQUATION
% =========================================================================

% Spatial Forcing term: -\nabla(\nabla \cdot \sigma)
Data.source_sigma   = {@(x,y)  (pi^2) * sin(pi*x).*sin(pi*y); 
                       @(x,y) -(pi^2) * cos(pi*x).*cos(pi*y);
                       @(x,y)  (pi^2) * cos(pi*x).*cos(pi*y);
                       @(x,y) -(pi^2) * sin(pi*x).*sin(pi*y)};
Data.source_sigma_t = {@(t) sin(2*t)};

% Temporal Forcing term: (1/\mu) * \partial_t dev(\sigma)
Data.source_sigma_d   = {@(x,y)  sin(pi*x).*sin(pi*y); 
                         @(x,y)  0.*x.*y;
                         @(x,y)  0.*x.*y;
                         @(x,y) -sin(pi*x).*sin(pi*y)};
Data.source_sigma_d_t = {@(t) 2*cos(2*t)};

% =========================================================================
% BOUNDARY CONDITIONS
% =========================================================================

% Dirichlet BC (Library 'D': imposes sigma.n for the paper's traction bounds)
Data.DirBCsigma      = {@(x,y)  sin(pi*x).*sin(pi*y);
                        @(x,y)  0.*x.*y;
                        @(x,y)  0.*x.*y;
                        @(x,y) -sin(pi*x).*sin(pi*y)};                     

% Neumann BC (Library 'N': imposes div(sigma) for the paper's velocity bounds)
% Partial derivative with respect to x (\partial_x \sigma)
Data.NeuBCsigma_x  = {@(x,y)  pi*cos(pi*x).*sin(pi*y); 
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y;
                      @(x,y) -pi*cos(pi*x).*sin(pi*y)};
 
% Partial derivative with respect to y (\partial_y \sigma)
Data.NeuBCsigma_y  = {@(x,y)  pi*sin(pi*x).*cos(pi*y); 
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y;
                      @(x,y) -pi*sin(pi*x).*cos(pi*y)};

% =========================================================================
% EXACT SOLUTION (For initial conditions & error tracking)
% =========================================================================

Data.sigma_ex    = {@(x,y)  sin(pi*x).*sin(pi*y);
                    @(x,y)  0.*x.*y;
                    @(x,y)  0.*x.*y;
                    @(x,y) -sin(pi*x).*sin(pi*y)};
Data.sigma_t_ex  = {@(t) sin(2*t)}; 
Data.sigma_dt_ex = {@(t) 2*cos(2*t)};

% forcing term for velocity and pressure recovery (not used in Sec 5.1)
Data.source_vel   = {@(x,y) 0.*x.*y; 
                     @(x,y) 0.*x.*y};
Data.source_vel_t = {@(t) 0*t};

Data.source_vel_d   = {@(x,y) 0.*x.*y; 
                       @(x,y) 0.*x.*y};
Data.source_vel_d_t = {@(t) 0*t};

Data.vel0 = {@(x,y)     0.*x.*y; 
             @(x,y)     0.*x.*y};