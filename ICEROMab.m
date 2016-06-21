% obj = ICEROMab(ab, mind, thick, filepath, [options], [NNMoptions])
%
% Uses the ABINT class to generate NLROMs of beam and plate structures. 


% Based on implicitcondensation_v5 by Rob Kuether
%
% Dependencies: npermutek.m     % Generates permutations
%               icloading.m     % Generate load cases
%               expansion.m     % Expand mebrane deformations
%               noregress_constrained.m % Fit nonlinear coefficients (constrained)
%               noregress.m     % Fit nonlinear coefficients (unconstrained)
%

classdef ICEROMab
    properties(GetAccess = 'public', SetAccess = 'public')
        abint; % ABINT object associated with this ROM
        mind; % Modeshape indices for the ROM
        thick; % Geometric thickness of the model
        options; % Option settings for the fit routine
        linDispLC;
        filepath; 
        springs;
        
        linModel; % The linear model [structure]
        phiReduced; % Reduced modal basis of the NLROM
        freqsReduced; % The frequencies of the reduced modal basis
        util; % Things we need/want to know but don't want to save. [structure]
        NLROM; % Structure of NLROM properties for later use
        x_rat; % No idea what this is
        phi_m; % membrane basis vectors (NDof x Nmembrane modes)
        nums_m; % permutation matrix describing q_membrane terms
        NNMoptions; % Structure for the NNM continuation options
    end
    
    methods
        function obj = ICEROMab(ab, thick, filepath, ...
                options, NNMoptions, springs)
            if nargin < 4
                options = struct;
            end
            if nargin < 5
                NNMoptions = struct;
            end
            if nargin < 6
                springs = [];
            end
            obj.NNMoptions = ICEROMab.parseNNMOptions(NNMoptions);
            obj.options = ICEROMab.parseOptions(options);
            obj.springs = springs;
            obj.abint = ab;
            obj.thick = thick;
            obj.filepath = filepath;
            
        end
        
        
        function obj = linearModel(obj, boundary, nmodes)
            %    obj = linearModel(obj, boundary, nmodes)
            %
            %    Get modes, mass, and stiffness matrix of the Abaqus model
            %    using specified boundary conditions.
            
            if (nargin < 3)
                nmodes = 20;
            end
            
            lm.bcDofAbs = [];
            for i = 1:size(boundary, 1)
                node = boundary(i, 1);
                bcRow = find(boundary(i, 2:7));
                bcRowAbs = obj.abint.getAbsDofNumber(node + bcRow(:)/10);
                lm.bcDofAbs = [lm.bcDofAbs; bcRowAbs];
            end
            
            
            lm.internalDofAbs = (1:length(obj.abint.dof))';
            for i = 1:length(lm.bcDofAbs)
                lm.internalDofAbs(lm.internalDofAbs == lm.bcDofAbs(i)) = [];
            end
            
            fprintf('Running Modal Analysis on %s\n', obj.abint.modelname);
            [lm.phi, lm.fn, lm.DOF] = obj.abint.modal(boundary, nmodes, obj.springs); 
            
            % Get mass and stiffness matrix
            fprintf('Generating Mass/Stiffness Matrices for %s\n', obj.abint.modelname);
            [lm.M, lm.K] = obj.abint.matrix(boundary, obj.springs);
            
            obj.linModel = lm;
        end
        
        
        function obj = genLoadCases(obj, mind, deflect)
            obj.mind = mind;
            
            % Run preliminary load cases; linear
            nM = length(obj.mind);
            Ml = obj.linModel.M(obj.linModel.internalDofAbs, obj.linModel.internalDofAbs);
            phi = obj.linModel.phi(obj.linModel.internalDofAbs, obj.mind);
            omega = obj.linModel.fn(obj.mind)*2*pi;
            
            
            if nargin < 3 % If no deflection vector supplied
                iceLoad = bsxfun(@times, obj.thick*Ml*phi, omega'.^2./max(abs(phi)));

                dispIce = obj.abint.static(obj.linModel.bcDofAbs, iceLoad, obj.linModel.internalDofAbs, 'nonlinear', obj.springs);
                
                maxIce = max(abs(dispIce))';
                phiNl = diag(pinv(phi)*dispIce);
                fprintf('Displacement Ratios; Initial:\n');
                for i = 1:nM
                    fprintf('\tMode %i: %.2f%% \n', i, (maxIce(i))/obj.thick*100);
                end
                % Now fit to a quadratic, cubic, or cubic + quadratic 1-mode
                % model
                frHat = diag(phi'*iceLoad);
                bEst = (diag(phi'*iceLoad) - (omega).^2.*phiNl)./phiNl.^3;

                % Use these numbers to estimate modal force level to obtain +/-
                % 15% of the linear response
                % f = frHat/omega^2; c = b/omega^2
                cubicEq = @(f, c) fzero(@(q) q + c*q^3 - f, 0);
                deflect = zeros(nM, 1);
                for i = 1:nM
                    c = bEst(i)/omega(i)^2;
                    startPoint = frHat(i)/omega(i)^2;
                    force = fzero(@(f) (f - cubicEq(f, c))^2/f^2 - 0.15^2, startPoint);
                    deflect(i) = max(abs(phi(:, i)))*force;
                end

            end
            
            % Otherwise just use the supplied deflection vector
            if (length(deflect) ~= length(obj.mind))
                error('Need 1 deflection value per mode');
            end
            deflect = deflect(:)*obj.thick;
            % Check results of this procedure
            iceLoad2 = bsxfun(@times, bsxfun(@times, deflect', Ml*phi), omega'.^2./max(abs(phi)));
            dispIce2 = obj.abint.static(obj.linModel.bcDofAbs, iceLoad2, obj.linModel.internalDofAbs, 'nonlinear', obj.springs);
            maxIce = max(abs(dispIce2))';

            fprintf('Displacement Ratios; Updated:\n');
            for i = 1:nM
                fprintf('\tMode %i: %.2f%% \n', i, (maxIce(i))/deflect(i)*100);
            end
            
            obj.linDispLC = deflect;
            
        end
        
        function obj = runLoadCases(obj, deflect)
            if (nargin == 2)
                if (length(deflect) ~= length(obj.mind)); error('Wrong size on deflections'); end
                obj.linDispLC = deflect;
            end
            
            omega = obj.linModel.fn(obj.mind)*2*pi;
            phi = obj.linModel.phi(obj.linModel.internalDofAbs, obj.mind);
            Ml = obj.linModel.M(obj.linModel.internalDofAbs, obj.linModel.internalDofAbs);
            
            % Obtain the loads
            [util.P, util.F] = icloading(omega/2/pi, ...
                phi, obj.linDispLC, Ml, [], 1, obj.options.quads, obj.options.triple, ...
                obj.options.rf, 0, obj.options.cf);
            

            xCheck = max(abs(phi*(bsxfun(@rdivide, phi'*util.F, omega.^2))));
            
            % (2) Run applied force static analysis for each load case
            util.disp = obj.abint.static(obj.linModel.bcDofAbs, util.F, obj.linModel.internalDofAbs, 'nonlinear', obj.springs);
            [obj.phi_m, obj.nums_m] = expansion(util.disp, phi); %#ok
            
            
            
            % (3) Perform NL coefficient fit
            [Knl, tlist, util.dispErr, util.forceErr] = noregress_constrained(omega/2/pi, phi, util.F, ...
                    util.disp, obj.options.quads, obj.options.triple, [], phi'*Ml);
                
%             [Knl, tlist, util.dispErr, util.forceErr] = noregress(omega/2/pi, phi, util.F, ...
%                     util.disp, obj.options.quads, obj.options.triple, [], phi'*Ml);
            [NLROM.N1, NLROM.N2, NLROM.Nlist]= Knl2Nash(Knl,tlist);
            NLROM.Knl = Knl;    
            NLROM.tlist = tlist;
            NLROM.Khat = diag(omega.^2);
            NLROM.Mhat = eye(size(NLROM.Khat));    
                
            obj.NLROM = NLROM;
            obj.util = util;
        end
        
        function diagnostics(obj, modeRange)
            if nargin < 2
                modeRange = [1, length(obj.linModel.fn)];
            end
            % Diagnostics
            % First load cases are positive and negative single mode forces:
            Lambda = diag((2*pi*obj.linModel.fn).^2);
            singleModes = 1:2*length(obj.mind);
            Q = Lambda\obj.linModel.phi'*obj.util.F(:, singleModes);
            
            Qnl = obj.linModel.phi'*obj.linModel.M*obj.util.disp(:, singleModes);
            
            ratio = Q./Qnl;
            [~, indices] = max(abs(ratio), [], 1);
            
            % Summary of single-mode displacement ratios
            fprintf('Summary of single-mode displacement ratios:\n')
            fprintf('\tMode\tLoad\tRatio [%%]\n');
            
            for i = 1:2*length(obj.mind)
                keyboard
                fprintf('\t%i\t\t%i\t\t%.2f\n', obj.mind(indices(i)), ...
                    obj.util.P(i, 1), 1/ratio(indices(i), i)*100);
            end
                    
            % Plot modes excited by each modal force.  If modes not in the ROM
            % are excited then the model may fail.
            mind_maxdisp = sign(obj.util.P(1:2*length(obj.mind),1)).*vec(obj.mind(abs(obj.util.P(1:2*length(obj.mind),1))));
            q = obj.linModel.phi.'*obj.linModel.M*obj.util.disp(:, 1:2*length(obj.mind));
            figure;
            plot(1:length(q), q(:, 1:length(obj.mind)), '.-', 'linewidth', 2); grid on;
            hold on; 
            plot(1:length(q), q(:, length(obj.mind) + 1:end), '+--', 'linewidth', 2);  
            hold off;
            legend(num2str(mind_maxdisp)); xlim(modeRange);
            xlabel('\bfMode Number'); ylabel('\bfModal Displacement');
            title('\bfModes Excited by Each Modal Force');
            
%             figure;
%             plot(1:size(obj.util.disp,2), obj.linModel.phiTM*obj.util.disp,'.');
%             legend(num2str(obj.linModel.fn.'));
%             xlabel('Load Case'); ylabel('Modal Displacement');
%             title('\bfModes Excited by Each Load Case');
        end
        
       
        function [fnl, N1_u, N2_u] = nlForce(obj, u)
            % Compute the nonlinear force resulting from a displacement u
            % according to the ROM's Nash-form coefficients.
            %
            % The Nash-form coefficients are given by N1, N2, and Nlist.
            % From Knl2Nash.m:
            %   In theory:
            %    Knl(u)= [1/2 * N1(u)+ 1/3 * N2(u,u)] * u
            %     and
            %    dKnl(u) = [N1(u) + N2(u,u)]*u;
            %      where N1(u) is N1 evaluated at u 
            %      where N2(u,u) is N2 evaluated at the Nlist quad combos of u
            %
            %    EXAMPLE evaluation of N1(u) and N2(u,u):
            %     u_nl = u(Nlist(:,1)).*u(Nlist(:,2));
            %     for jj=1:Neqn;
            %       N1_u(:,jj)= N1(:,:,jj) * u;
            %       N2_u(:,jj)= N2(:,:,jj) * u_nl;
            %     end;
            N1 = obj.NLROM.N1; N2 = obj.NLROM.N2; Nlist = obj.NLROM.Nlist;
            u_nl = u(Nlist(:, 1)).*u(Nlist(:, 2));
            n = length(u); N1_u = zeros(n); N2_u = zeros(n);
            for j = 1:n
                N1_u(:, j) = N1(:, :, j)*u;
                N2_u(:, j) = N2(:, :, j)*u_nl;
            end
            fnl = (1/2*N1_u + 1/3*N2_u)*u;
            %            keyboard
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% UTILITY METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function save(obj)
            rom = obj; %#ok
            mindstr = strrep(num2str(obj.mind), '  ', '-');
            filename = sprintf('%s_qu%itr%i_%s', obj.abint.modelname, ...
                          obj.options.quads, obj.options.triple, mindstr);
            save(fullfile(obj.filepath, filename), 'rom');
        end
        
        function sys = defsys(obj, mode)
            %    sys = defsys(obj)
            %
            % Create the "sys" object expected by NNMcont using the
            % included "NNMoptions" structure
            
            % Linear System:
            sys.Klin = obj.NLROM.Khat;
            sys.Mlin = obj.NLROM.Mhat;
            
            % Nonlinearities:
            sys.nl = NL_ROM(obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist);
            
            % Options
            opt = obj.NNMoptions;
            sys.norm = opt.norm;
            sys.NumMeth = opt.NumMeth;
            sys.TFS = opt.TFS;
            sys.PrecS = opt.PrecS;
            sys.itoptS = opt.itoptS;
            sys.itmaxS = opt.itoptS;
            sys.RelTolode = opt.RelTolode;
            sys.AbsTolode = opt.AbsTolode;
            sys.NRprec = opt.NRprec;
            sys.alpha_m = opt.alpha_m;
            sys.h = opt.h;
            sys.hMax = opt.hMax;
            sys.hMin = opt.hMin;
            sys.betamax = opt.betamax;
            sys.betamin = opt.betamin;
            sys.wstop = opt.wstop;
            sys.shootinPeriod = opt.shootingPeriod;
            
            % Mode and filename
            mindstr = strrep(num2str(obj.mind), '  ', '-');
            sys.filename = sprintf('%s_%s_NNM-%i', obj.abint.modelname, mindstr, obj.mind(mode));
            sys.mode = mode;
            
        
        end
    end
    
    methods(Static)
        
        function opt = parseOptions(options)
            % Options
            % Scaling Method for Applied Loads? 1- Constant force, [0-
            % Constant displacement]
            try opt.cf = options.cf; catch err; opt.cf = 0; end
            % Reduce scale factors of 2 and 3 combos? [1- Yes], 0- No
            try opt.rf = options.rf; catch err; opt.rf = 1; end
            % do constrained fit in no_regress? [1- Yes], 0- No
            try opt.constr = options.constr; catch err; opt.constr = 1; end
            % include quadratics in no_regress? [1- Yes], 0- No
            try opt.quads = options.quads; catch err; opt.quads = 1; end
            % include triple cubics in no_regress? [1- Yes], 0- No
            try opt.triple = options.triple; catch err; opt.triple = 1; end
            % include membrane expansion? [1- Yes], 0- No
            try opt.memb = options.memb; catch err; opt.memb = 1; end
            % Overwrite linear model? [1 - Yes], 0 - No
            try opt.ovrwrt = options.overwrite; catch err; opt.ovrwrt = 1; end
        end
        
        function opt = parseNNMOptions(opt)
            % Options for the NNMcont continuation code
            if ~isfield(opt, 'norm')
                opt.norm = 5e-8; % Make this smaller if having convergence issues
            end
            if ~isfield(opt,'NumMeth')
                opt.NumMeth='NEWMARK'; % Numerical method
            end
            if ~isfield(opt,'TFS')
                opt.TFS=500; % Time steps per period integration
            end
            if ~isfield(opt,'PrecS'),
                opt.PrecS = 1E-6; % 
            end
            if ~isfield(opt,'itoptS'),
                opt.itoptS = 3; % No idea what this does
            end
            if ~isfield(opt,'itmaxS'),
                opt.itmaxS = 1; % "h_max???"
            end
            if ~isfield(opt, 'RelTolode'),
                opt.RelTolode = 1.0e-7; % ODE relative tolerance
            end
            if ~isfield(opt,'AbsTolode')
                opt.AbsTolode = 1.0e-7; % ODE absolute tolerance
            end
            if ~isfield(opt,'NRprec'),
                opt.NRprec = 1.0e-7; % Don't know what this does 
            end
            if ~isfield(opt,'alpha_m'),
                opt.alpha_m = 0; % Don't know what this does
            end
            if ~isfield(opt,'alpha_f'); 
                opt.alpha_f = 0.01; % Don't know what this does
            end
            
            if ~isfield(opt,'h'); opt.h = 1E-5; end
            if ~isfield(opt, 'hMax'); opt.hMax = 1E-2; end
            if ~isfield(opt, 'hMin'); opt.hMin = 1E-9; end
            if ~isfield(opt, 'betamax'); opt.betamax = 85; end
            if ~isfield(opt, 'betamin'); opt.betamin = 0; end
            if ~isfield(opt, 'wstop'); opt.wstop = []; end
            if ~isfield(opt, 'shootingPeriod'); opt.shootingPeriod = 1; end
            if strcmpi(opt.shootingPeriod, 'HALF'); opt.shootingPeriod = 1; end
            if strcmpi(opt.shootingPeriod, 'FULL'); opt.shootingPeriod = 2; end
            if (opt.shootingPeriod ~= 1) && (opt.shootingPeriod ~= 2)
                error('Incorrect specification for shootingPeriod')
            end
        end
    end
end




