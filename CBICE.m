% CIBCE(abint, mind, thick, linDispLC, workingDir, filepath, options, NNMoptions)
%
%
% Based on implicitcondensation_v5 by Rob Kuether
%
% Dependencies: npermutek.m     % Generates permutations
%               icloading.m     % Generate load cases
%               expansion.m     % Expand mebrane deformations
%               noregress_constrained.m % Fit nonlinear coefficients (constrained)
%               noregress.m     % Fit nonlinear coefficients (unconstrained)
%

classdef CBICE
    properties(GetAccess = 'public', SetAccess = 'public')
        
        abint; % ABINT interface for this CBROM
        NLid; % Unique ID for the nonlinear ROM
        mind; % Modeshape indices for the ROM
        thick; % Geometric thickness of the model
        options; % Option settings for the fit routine
        matDir; % Base directory for file storage
        part; % Part designation (depends on makeup of the FE model)
        linDispLC;
        workingDir; % FEA working directory
        filepath;
        linModel; % The linear model [structure]
        phiReduced; % Reduced modal basis of the NLROM
        freqsReduced; % The frequencies of the reduced modal basis
        util; % Things we need/want to know but don't want to save. [structure]
        NLROM; % Structure of NLROM properties for later use
        phi_m; % membrane basis vectors (NDof x Nmembrane modes)
        nums_m; % permutation matrix describing q_membrane terms
        NNMoptions; % Structure for the NNM continuation options
        
        normFactor; % Normalization factor for each flexible mode
        
    end
    
    methods
        function obj = CBICE(abint, thick, workingDir, ...
                filepath, options, NNMoptions)
            if nargin < 7
                options = struct;
            end
            if nargin < 8
                NNMoptions = struct;
            end
            obj.NNMoptions = ICEROM.parseNNMOptions(NNMoptions);
            obj.abint = abint;
            obj.filepath = filepath;
            obj.thick = thick;
            obj.options = ICEROM.parseOptions(options);
            obj.workingDir = workingDir;
            
        end
        
        
        function [obj, omegaN, MkbHat, MbbHat, KbbHat] = makeCB(obj, nModes, bcs, boundaries, method)
            %    [obj, omegaN, MkbHat, MbbHat, KbbHat] = makeCB(obj, nModes, bcs, boundaries, method)
            %
            % Make a Craig-Bampton substructure of a model. Requires the
            % number of elastic modes for retention (will change to
            % frequency range later), the actual boundary conditions of the
            % part, and the boundary nodes used for attachment to generate
            % fixed-interface modes.
            %
            % "Method" is "abaqus" [default] or "matlab". "matlab" uses
            % MATLAB's built-in routines to perform modal analyses of the
            % parts; "abaqus" calls out to Abaqus FEA software.
            
            if ~(strcmpi(method, 'matlab') || strcmpi(method, 'abaqus'))
                error('Invalid method selected for CBICE.makeCB');
            end
            % If matrices don't exist yet, get 'em
            if isempty(obj.linModel)
                [obj.linModel.M, obj.linModel.K] = obj.abint.matrix(bcs); 
            end
            
            % Sort out which DOF are which (store in ABSOLUTE DOF
            % numbering)
            % Base boundary conditions
            bcDof = [];
            for i = 1:size(bcs, 1)
                node = bcs(i, 1);
                constraints = find(bcs(i, 2:7));
                bcDof = [bcDof; node + constraints(:)/10];
            end
            % Relate to absolute DOF numbers
            bcDofAbs = obj.abint.getAbsDofNumber(bcDof);
            
            % Connector boundary DOF
            connDof = [];
            for i = 1:size(boundaries, 1)
                node = boundaries(i, 1);
                constraints = find(boundaries(i, 2:7));
                connDof = [connDof; node + constraints(:)/10];
            end
            % Relate to absolute DOF numbers
            connDofAbs = obj.abint.getAbsDofNumber(connDof);
            
            % Get internal DOF numbering
            internalDofAbs = (1:length(obj.abint.dof))';
            remove = [bcDofAbs; connDofAbs];
            for i = 1:length(remove)
                internalDofAbs(internalDofAbs == remove(i)) = [];
            end
            
            % Store these in "linModel"
            obj.linModel.bcDofAbs = bcDofAbs;
            obj.linModel.connDofAbs = connDofAbs;
            obj.linModel.internalDofAbs = internalDofAbs;
            dofOrder = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs; obj.linModel.bcDofAbs];
            [~, obj.linModel.dofOrder] = sort(dofOrder);
                        
            % Get fixed or free-interface modes
            fixedSet = [bcs; boundaries]; 
            
            if strcmpi(method, 'matlab') % If MATLAB is selected
                for i = 1
                    % Partition system matrices to internal nodes
                    M = obj.linModel.M(internalDofAbs, internalDofAbs);
                    K = obj.linModel.K(internalDofAbs, internalDofAbs);
                    if (length(nModes) == 2) % Frequency range option
                        if ~exist('obj.linModel.omegaN')
                            warning('CBICE.makeCB does not support "MATLAB" option with frequency inputs for first call.');
                            warning('Calling Abaqus. "MATLAB" may be used for subsequent calls to CBICE.makeCB');
                            method = 'abaqus';
                            break;
                        else
                            nModes = length(obj.linModel.omegaN);
                        end
                    end
                    fprintf('Finding %i modes for %s\n', nModes, obj.abint.modelname);
                    [phiFI, fnFI] = eigs(K, M, nModes, 'SM');
                    fnFI = sqrt(diag(fnFI))/2/pi;
                    [fnFI, sortInd] = sort(fnFI);
                    phiFI = phiFI(:, sortInd);
                end
            end
       
            
            if strcmpi(method, 'abaqus') % If Abaqus is selected
                [phiFI, fnFI] = obj.abint.modal(fixedSet, nModes);
                if isempty(fnFI)
                    warning('No modes detected in range [%.0f, %.0f]; re-running with 5 FI modes', ...
                        nModes(1), nModes(2));
                    [phiFI, fnFI] = obj.abint.modal(fixedSet, 5);
                end
                phiFI = phiFI(internalDofAbs, :);
            end
            
            omegaN = fnFI*2*pi;
            
            % Get characteristic constraint modes [Use stiffness matrices]
            Kii = obj.linModel.K(internalDofAbs, internalDofAbs);
            Kib = obj.linModel.K(internalDofAbs, connDofAbs);
            psiIB = -Kii\Kib; % Very badly scaled!
            
            % Save modal matrices to model
            obj.normFactor = max(abs(phiFI), [], 1); % Normalization factor for each mode
            obj.linModel.phiFI = phiFI;
            obj.linModel.omegaN = omegaN;
            obj.linModel.psiIB = psiIB;
            
            % Build mass "hat" matrices
            Mii = obj.linModel.M(internalDofAbs, internalDofAbs);
            Mib = obj.linModel.M(internalDofAbs, connDofAbs);
            MkbHat = phiFI'*(Mii*psiIB + Mib);
            
            Mbb = obj.linModel.M(connDofAbs, connDofAbs);
            MbbHat = psiIB'*Mii*psiIB + Mib'*psiIB + psiIB'*Mib + Mbb;
            
            % Build stiffness "hat" matrix
            Kbb = obj.linModel.K(connDofAbs, connDofAbs);
            KbbHat = psiIB'*Kii*psiIB + Kib'*psiIB + psiIB'*Kib + Kbb;

            % Save mass and stiffness matrices
            obj.linModel.MkbHat = MkbHat;
            obj.linModel.MbbHat = MbbHat;
            obj.linModel.KbbHat = KbbHat;
            
        end
        
        
        function [obj, data] = qrFit(obj, qrInds, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            %    [obj, data] = qrFit(obj, qrInds, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            % 
            % Use a subset of FI and CC modes in order to form a nonlinear
            % fit of the component. FI and CC modes are indexed from one.
            % Include springs if necessary. Quads and triples default to 1
            % (active), "flist" is only used if supplied, otherwise the
            % curve fit routine will select it automatically.
            
            % (1) Form variables of convenience
            nFI = length(obj.linModel.omegaN); nCC = size(ccBasisFull, 2);
            nB = length(obj.linModel.connDofAbs);
            
            
            % (2) Obtain physical basis for each set of modes
            % First obtain the FI/Constraint basis
            fiBasis = [obj.linModel.phiFI; zeros(nB, nFI)];
            ccBasis = [obj.linModel.psiIB; 
                          eye(nB, nB)]*ccBasisFull((nFI + 1):end, :);
            basis = [fiBasis, ccBasis];
            
            % Get loaded DOF partition and associated system matrices
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            Kl = obj.linModel.K;
            for i = 1:length(springs.dof)
                dof = springs.dof(i); k = springs.k(i);
                Kl(dof, dof) = Kl(dof, dof) + k;
            end
            Kl = Kl(loadDof, loadDof);
            

            % (2) Perform a QR decomposition of the FI/CC basis to get a reduced set of
            % NL modes
            [Q, R] = qr(basis, 0);
            loadBasis = Q(:, qrInds);
            theta = R;
            
            tic;
            if isempty(fitData)
            
            % (4) Use icloading to generate loadings for the basis vectors
            [P, F] = icloading_v2(loadBasis, thick*obj.thick, eye(length(loadBasis)), Kl, 1, quads, triples, 1);
            % Zero out forces on boundary nodes
            boundaryNodes = (length(obj.linModel.internalDofAbs) + ...
                 1):(length(obj.linModel.internalDofAbs) + length(obj.linModel.connDofAbs));
%             F(boundaryNodes, :) = 0; warning('Boundary forces zeroed in
%             SVD fit');
            
            % (5) Send to Abaqus
            
            [displacement, cf] = obj.abint.static(obj.linModel.bcDofAbs, F, loadDof, 'nonlinear', springs);
            else
                F = fitData.F;
                P = fitData.P;
                displacement = fitData.displacement;
            end
            
            
            maxDisp = max(abs(displacement));
            for i = 1:length(qrInds)
                fprintf('Displacement ratio SVD %i = %.3f%%\n', i, maxDisp(i)/(thick*obj.thick)*100);
            end
            
                
            % (6) Perform fit
            [Knl, tlist, dispErr, forceErr] = noregress(loadBasis'*Kl*loadBasis, loadBasis, F, ... % Can we fit with a reduced set of U vectors?
                    displacement, quads, triples, flist);
            fitTime = toc;
%             keyboard

            % (7) Zero out all other NL terms
            emptyModes = (length(qrInds) + 1):(nFI + nCC);
            tlistEmpty = repmat(emptyModes', 1, 3);
            nEmpty = size(tlistEmpty, 1);
            obj.NLROM.tlist = [tlist; tlistEmpty];
            obj.NLROM.Knl = blkdiag(Knl, zeros(nEmpty));
            [obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist] = Knl2Nash(obj.NLROM.Knl, obj.NLROM.tlist);
            obj.NLROM.theta = theta;
            obj.NLROM.basis = basis;
            
            % Optional outputs
            data.basis = basis; data.F = F; data.P = P; data.displacement = displacement;
            data.loadBasis = loadBasis; data.theta = theta; data.fitTime = fitTime;
            data.dispErr = dispErr; data.forceErr = forceErr; data.Kl = Kl;
            
        end
        
        function [obj, data] = qrFitIS(obj, refLoc, ab, nodeset, bcset, qrInds, ccBasisFull, thick, quads, triples, flist, fitData)
            %     [obj, data] = qrFitIS(obj, refLoc, ab, nodeset, bcset, qrInds, ccBasisFull, thick, quads, triples, flist, fitData)
            % 
            % Use a subset of FI and CC modes in order to form a nonlinear
            % fit of the component. FI and CC modes are indexed from one.
            % Include springs if necessary. Quads and triples default to 1
            % (active), "flist" is only used if supplied, otherwise the
            % curve fit routine will select it automatically.
            %
            % This version uses the "in situ" FEA fit process to get a
            % curve fit.
            
            % (1) Form variables of convenience
            nFI = length(obj.linModel.omegaN); nCC = size(ccBasisFull, 2);
            nB = length(obj.linModel.connDofAbs);
            
            
            % (2a) Obtain physical basis for each set of modes
            % First obtain the FI/Constraint basis
            fiBasis = [obj.linModel.phiFI; zeros(nB, nFI)];
            ccBasis = [obj.linModel.psiIB; 
                          eye(nB, nB)]*ccBasisFull((nFI + 1):end, :);
            basis = [fiBasis, ccBasis];
            
            % Get loaded DOF partition and associated system matrices
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            
            Kl = obj.linModel.K(loadDof, loadDof);
            
            % (2b) Perform a QR decomposition of the FI/CC basis to get a reduced set of
            % NL modes
            [Q, R] = qr(basis, 0);
            loadBasis = Q(:, qrInds);
            theta = R;

            tic
            if isempty(fitData) % No fitData present
                % (3) Use icloading to generate loadings for the basis vectors
                [P, F] = icloading_v2(loadBasis, thick*obj.thick, eye(length(loadBasis)), Kl, 1, quads, triples, 1);
                % Zero out forces on boundary nodes
                boundaryNodes = (length(obj.linModel.internalDofAbs) + ...
                     1):(length(obj.linModel.internalDofAbs) + length(obj.linModel.connDofAbs));
                F(boundaryNodes, :) = 0;

                % (4a) Line up the nodeset and ABINT object supplied with
                % the desired load DOF vector using node locations.
                % - Get the node locations of the component load DOF
                compDof = obj.abint.getDofNumber(loadDof);
                compNodes = unique(floor(compDof), 'stable');
                compNodesAbs = obj.abint.getAbsNodeNumber(compNodes);
                compXYZ = obj.abint.nodes(compNodesAbs, 2:4);
                
                % - Get the node locations of the ABINT nodeset DOF
                abNodesAbs = ab.getAbsNodeNumber(nodeset);
                abXYZ = ab.nodes(abNodesAbs, 2:4);
                
                % Get the offset between coordinate systems
                offset = obj.abint.nodes(refLoc(1), 2:4) - ...
                         ab.nodes(refLoc(2), 2:4);
                abXYZOffset = bsxfun(@plus, abXYZ, offset(:)');
                
                % Now run through and match the DOF at each node
                loadDofAb = zeros(size(loadDof));
                for i = 1:length(compNodes)
                    compNode = compNodes(i);
                    compLoc = compXYZ(i, :);
                    d = bsxfun(@minus, abXYZOffset, compLoc);
                    dnorm = dot(d, d, 2);
                    [distance, abNodeAbs] = min(dnorm);
                    if (distance > 1E-6)
                        error('Comp Node %i out of tolerance; |d| = %.3E\n', ...
                            compNode, distance);
                    end
                    abNode = nodeset(abNodeAbs);
                    for j = 1:6 % This is written only for beams and shells
                        dofNum = (i - 1)*6 + j;
                        loadDofAb(dofNum) = abNode + j/10;
                    end
                end
                loadDofAbs = ab.getAbsDofNumber(loadDofAb); 
                
                % Figure out the BC set - not written generally at all
                % right now, just assumes each node is fully fixed.
                bcDOF = bsxfun(@plus, bcset, [1:6]/10)';
                bcDOF = bcDOF(:);
                bcDOFAbs = ab.getAbsDofNumber(bcDOF);
                
                % (4b) Run Abaqus job on the assembly nodeset of interest
                [displacement, cf] = ab.static(bcDOFAbs, F, loadDofAbs, 'nonlinear');
%                 [displacement, cf] = obj.abint.static(obj.linModel.connDofAbs, F, loadDof, 'nonlinear', springs);
                
            else
                F = fitData.F;
                P = fitData.P;
                displacement = fitData.displacement;
            end
            
            
            % (5) Perform fit
            lambda = loadBasis'*Kl*loadBasis;
%             [Knl, tlist, dispErr, forceErr] = noregress(lambda, basis, F, ...
%                     displacement, quads, triples, flist);
                
            if isequal(size(lambda), [1, 1])
                lambda = sqrt(lambda)/2/pi;
            end
            [Knl, tlist, dispErr, forceErr] = noregress(lambda, loadBasis, F, ...
                    displacement, quads, triples, flist);
            fitTime = toc;
            
            % (6) Add all non-included modes
            t = tlist(:); 
            discardModes = 1:(nFI + nCC);
            for i = 1:length(qrInds)
                ti = (t == i);
                t(ti) = qrInds(i);
                discardModes(i) = 0;
            end
            tlistActual = reshape(t, size(tlist));
            discardModes(discardModes == 0) = [];
            tlistEmpty = repmat(discardModes(:), 1, 3);

            
            KnlActual = zeros(nFI + nCC, size(Knl, 2));
            KnlActual(qrInds, :) = Knl;
            
            nEmpty = size(tlistEmpty, 1);
            obj.NLROM.tlist = [tlistActual; tlistEmpty];
            
            obj.NLROM.Knl = [KnlActual, zeros(nFI + nCC, nEmpty)];
            [obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist] = Knl2Nash(obj.NLROM.Knl, obj.NLROM.tlist);
            obj.NLROM.dispErr = dispErr; obj.NLROM.forceErr = forceErr;
            obj.NLROM.fitTime = fitTime;
            obj.NLROM.theta = theta;
            obj.NLROM.basis = basis;
            
            % Optional outputs
            data.basis = basis; data.F = F; data.P = P; data.displacement = displacement;
            data.fitTime = fitTime; data.Kl = Kl;
            data.dispErr = dispErr; data.forceErr = forceErr;
            
         end
        
        function [obj, data] = svdFit(obj, svdInd, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            %    [obj, data] = svdFit(obj, svdInd, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            % 
            % Use a subset of FI and CC modes in order to form a nonlinear
            % fit of the component. FI and CC modes are indexed from one.
            % Include springs if necessary. Quads and triples default to 1
            % (active), "flist" is only used if supplied, otherwise the
            % curve fit routine will select it automatically.
            
            % (1) Form variables of convenience
            nFI = length(obj.linModel.omegaN); nCC = size(ccBasisFull, 2);
            nB = length(obj.linModel.connDofAbs);
            nSvd = nFI + nCC;
            
            
            % (2) Obtain physical basis for each set of modes
            % First obtain the FI/Constraint basis
            fiBasis = [obj.linModel.phiFI; zeros(nB, nFI)];
            ccBasis = [obj.linModel.psiIB; 
                          eye(nB, nB)]*ccBasisFull((nFI + 1):end, :);
            basis = [fiBasis, ccBasis];
            
            % Get loaded DOF partition and associated system matrices
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            Kl = obj.linModel.K;
            for i = 1:length(springs.dof)
                dof = springs.dof(i); k = springs.k(i);
                Kl(dof, dof) = Kl(dof, dof) + k;
            end
            Kl = Kl(loadDof, loadDof);
            

            % (2) Perform an SVD of the FI/CC basis to get a reduced set of
            % NL modes
            [U, S, V] = svd(basis, 'econ');
            if isempty(svdInd) % This is the input option to examine the SVD results
                plot(diag(S));
                keyboard
            end
            loadBasis = U(:, svdInd);
            theta = S*V';
            
            
            tic
            if isempty(fitData) % No fitData present
                % (3) Use icloading to generate loadings for the basis vectors
                [P, F] = icloading_v2(loadBasis, thick*obj.thick, eye(length(loadBasis)), Kl, 1, quads, triples, 1);
                % Zero out forces on boundary nodes
%                 boundaryNodes = (length(obj.linModel.internalDofAbs) + ...
%                      1):(length(obj.linModel.internalDofAbs) + length(obj.linModel.connDofAbs));
    %             F(boundaryNodes, :) = 0;

                % (4) Send to Abaqus
                [displacement, cf] = obj.abint.static(obj.linModel.bcDofAbs, F, loadDof, 'nonlinear', springs);
                
            else
                F = fitData.F;
                P = fitData.P;
                displacement = fitData.displacement;
            end
            
            maxDisp = max(abs(displacement));
            for i = 1:length(svdInd)
                if length(thick) > 1; t = thick(i); 
                else t = thick;
                end
                fprintf('Displacement ratio SVD %i = %.3f%%\n', svdInd(i), maxDisp(i)/(t*obj.thick)*100);
            end
            
            % (6) Perform fit
            keyboard
            [Knl, tlist, dispErr, forceErr] = noregress(loadBasis'*Kl*loadBasis, loadBasis, F, ... % Can we fit with a reduced set of U vectors?
                    displacement, quads, triples, flist);
            fitTime = toc;
%             keyboard
            % (7) Zero out all other NL terms
            
            emptyModes = 1:nSvd;
            for i = emptyModes
                if sum(svdInd == i)
                    emptyModes(i) = nan;
                end
            end
            emptyModes(isnan(emptyModes)) = [];
            
            tlistEmpty = repmat(emptyModes', 1, 3);
            nEmpty = size(tlistEmpty, 1);
            obj.NLROM.tlist = [tlist; tlistEmpty];
            obj.NLROM.Knl = blkdiag(Knl, zeros(nEmpty));
            [obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist] = Knl2Nash(obj.NLROM.Knl, obj.NLROM.tlist);
            obj.NLROM.theta = theta;
            obj.NLROM.basis = basis;
            
            % Optional outputs
            data.basis = basis; data.F = F; data.P = P; data.displacement = displacement;
            data.loadBasis = loadBasis; data.theta = theta; data.fitTime = fitTime;
            data.dispErr = dispErr; data.forceErr = forceErr; data.Kl = Kl;
            
        end
        
         function [obj, data] = ficcFitIS(obj, refLoc, ab, nodeset, bcset, fiMind, ccMind, ccBasisFull, thick, quads, triples, flist, fitData)
            %     [obj, data] = ficcFitIS(obj, refLoc, ab, nodeset, bcset, fiMind, ccMind, ccBasisFull, thick, quads, triples, flist, fitData)
            % 
            % Use a subset of FI and CC modes in order to form a nonlinear
            % fit of the component. FI and CC modes are indexed from one.
            % Include springs if necessary. Quads and triples default to 1
            % (active), "flist" is only used if supplied, otherwise the
            % curve fit routine will select it automatically.
            %
            % This version uses the "in situ" FEA fit process to get a
            % curve fit.
            
            % (1) Form variables of convenience
            nFI = length(obj.linModel.omegaN); nCC = size(ccBasisFull, 2);
            nM = length(fiMind); nB = length(obj.linModel.connDofAbs);
            
            
            % (2) Obtain physical basis for each set of modes
            % First obtain the FI/Constraint basis
            fiBasis = [obj.linModel.phiFI(:, fiMind); zeros(nB, nM)];
            ccBasis = [obj.linModel.psiIB; 
                          eye(nB, nB)]*ccBasisFull((nFI + 1):end, ccMind);
            basis = [fiBasis, ccBasis];
            
            % Get loaded DOF partition and associated system matrices
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            
            Kl = obj.linModel.K(loadDof, loadDof);
            
%             keyboard
            tic
            if isempty(fitData) % No fitData present
                % (3) Use icloading to generate loadings for the basis vectors
                [P, F] = icloading_v2(basis, thick*obj.thick, Kl, Kl, 1, quads, triples, 1);
                % Zero out forces on boundary nodes
                boundaryNodes = (length(obj.linModel.internalDofAbs) + ...
                     1):(length(obj.linModel.internalDofAbs) + length(obj.linModel.connDofAbs));
                F(boundaryNodes, :) = 0;

                % (4a) Line up the nodeset and ABINT object supplied with
                % the desired load DOF vector using node locations.
                % - Get the node locations of the component load DOF
                compDof = obj.abint.getDofNumber(loadDof);
                compNodes = unique(floor(compDof), 'stable');
                compNodesAbs = obj.abint.getAbsNodeNumber(compNodes);
                compXYZ = obj.abint.nodes(compNodesAbs, 2:4);
                
                % - Get the node locations of the ABINT nodeset DOF
                abNodesAbs = ab.getAbsNodeNumber(nodeset);
                abXYZ = ab.nodes(abNodesAbs, 2:4);
                
                % Get the offset between coordinate systems
                offset = obj.abint.nodes(refLoc(1), 2:4) - ...
                         ab.nodes(refLoc(2), 2:4);
                abXYZOffset = bsxfun(@plus, abXYZ, offset(:)');
                
                % Now run through and match the DOF at each node
                loadDofAb = zeros(size(loadDof));
                for i = 1:length(compNodes)
                    compNode = compNodes(i);
                    compLoc = compXYZ(i, :);
                    d = bsxfun(@minus, abXYZOffset, compLoc);
                    dnorm = dot(d, d, 2);
                    [distance, abNodeAbs] = min(dnorm);
                    if (distance > 1E-6)
                        error('Comp Node %i out of tolerance; |d| = %.3E\n', ...
                            compNode, distance);
                    end
                    abNode = nodeset(abNodeAbs);
                    for j = 1:6 % This is written only for beams and shells
                        dofNum = (i - 1)*6 + j;
                        loadDofAb(dofNum) = abNode + j/10;
                    end
                end
                loadDofAbs = ab.getAbsDofNumber(loadDofAb); 
                
                % Figure out the BC set - not written generally at all
                % right now, just assumes each node is fully fixed.
                bcDOF = bsxfun(@plus, bcset, [1:6]/10)';
                bcDOF = bcDOF(:);
                bcDOFAbs = ab.getAbsDofNumber(bcDOF);
                
                % (4b) Run Abaqus job on the assembly nodeset of interest
                [displacement, cf] = ab.static(bcDOFAbs, F, loadDofAbs, 'nonlinear');
%                 [displacement, cf] = obj.abint.static(obj.linModel.connDofAbs, F, loadDof, 'nonlinear', springs);
                
            else
                F = fitData.F;
                P = fitData.P;
                displacement = fitData.displacement;
            end
            
            
            % (5) Perform fit
            lambda = basis'*Kl*basis;
%             [Knl, tlist, dispErr, forceErr] = noregress(lambda, basis, F, ...
%                     displacement, quads, triples, flist);
                
            phiInv = lambda\basis'*Kl;
            if isequal(size(lambda), [1, 1])
                lambda = sqrt(lambda)/2/pi;
            end
            [Knl, tlist, dispErr, forceErr] = noregress(lambda, basis, F, ...
                    displacement, quads, triples, flist, phiInv);
            fitTime = toc;
            
            % (6) Add all non-included modes
            t = tlist(:); 
            for i = ccMind
                tCC = (t == (nM + i));
                t(tCC) = nFI + i;
            end
            tlistActual = reshape(t, size(tlist));
            allFIModes = 1:nFI;
            allCCModes = (nFI + 1):(nFI + nCC);
            allFIModes(fiMind) = []; allCCModes(ccMind) = [];
            
            KnlActual = zeros(nFI + nCC, size(Knl, 2));
            KnlActual([fiMind(:); nFI + ccMind(:)], :) = Knl;

            tlistEmpty = repmat([allFIModes, allCCModes]', 1, 3);
            nEmpty = size(tlistEmpty, 1);
            obj.NLROM.tlist = [tlistActual; tlistEmpty];
            
            obj.NLROM.Knl = [KnlActual, zeros(nFI + nCC, nEmpty)];
            [obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist] = Knl2Nash(obj.NLROM.Knl, obj.NLROM.tlist);
            obj.NLROM.dispErr = dispErr; obj.NLROM.forceErr = forceErr;
            obj.NLROM.fitTime = fitTime;
            
            % Optional outputs
            data.basis = basis; data.F = F; data.P = P; data.displacement = displacement;
            data.fitTime = fitTime; data.Kl = Kl;
            data.dispErr = dispErr; data.forceErr = forceErr;
            
         end
        
        function [obj, data] = ficcFit(obj, fiMind, ccMind, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            %    [obj, data] = ficcFit(obj, fiMind, ccMind, ccBasisFull, springs, thick, quads, triples, flist, fitData)
            % 
            % Use a subset of FI and CC modes in order to form a nonlinear
            % fit of the component. FI and CC modes are indexed from one.
            % Include springs if necessary. Quads and triples default to 1
            % (active), "flist" is only used if supplied, otherwise the
            % curve fit routine will select it automatically.
            
            % (1) Form variables of convenience
            nFI = length(obj.linModel.omegaN); nCC = size(ccBasisFull, 2);
            nM = length(fiMind); nB = length(obj.linModel.connDofAbs);
            
            
            % (2) Obtain physical basis for each set of modes
            % First obtain the FI/Constraint basis
            fiBasis = [obj.linModel.phiFI(:, fiMind); zeros(nB, nM)];
            ccBasis = [obj.linModel.psiIB; 
                          eye(nB, nB)]*ccBasisFull((nFI + 1):end, ccMind);
            basis = [fiBasis, ccBasis];
            
            % Get loaded DOF partition and associated system matrices
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            
            Ksprings = blkdiag(zeros(length(obj.linModel.internalDofAbs)), ...
                               diag(springs.k));
            Kl = obj.linModel.K(loadDof, loadDof);
            
%             keyboard
            tic
            if isempty(fitData) % No fitData present
                % (3) Use icloading to generate loadings for the basis vectors
                [P, F] = icloading_v2(basis, thick*obj.thick, Kl, Kl + Ksprings, 1, quads, triples, 1);
                % Zero out forces on boundary nodes
                boundaryNodes = (length(obj.linModel.internalDofAbs) + ...
                     1):(length(obj.linModel.internalDofAbs) + length(obj.linModel.connDofAbs));
                F(boundaryNodes, :) = 0;

                % (4) Send to Abaqus
                [displacement, cf] = obj.abint.static(obj.linModel.bcDofAbs, F, loadDof, 'nonlinear', springs);
%                 [displacement, cf] = obj.abint.static(obj.linModel.connDofAbs, F, loadDof, 'nonlinear', springs);
                
            else
                F = fitData.F;
                P = fitData.P;
                displacement = fitData.displacement;
            end
            
            
            % (5) Perform fit
            lambda = basis'*Kl*basis;
            
%             [Knl, tlist, dispErr, forceErr] = noregress(lambda, basis, F, ...
%                     displacement, quads, triples, flist);
                
            phiInv = lambda\basis'*Kl;
            if isequal(size(lambda), [1, 1])
                lambda = sqrt(lambda)/2/pi;
            end
            [Knl, tlist, dispErr, forceErr] = noregress(lambda, basis, F, ...
                    displacement, quads, triples, flist, phiInv);
            fitTime = toc;
            
            % (6) Add all non-included modes
            t = tlist(:); 
            for i = ccMind
                tCC = (t == (nM + i));
                t(tCC) = nFI + i;
            end
            tlistActual = reshape(t, size(tlist));
            allFIModes = 1:nFI;
            allCCModes = (nFI + 1):(nFI + nCC);
            allFIModes(fiMind) = []; allCCModes(ccMind) = [];
            
            KnlActual = zeros(nFI + nCC, size(Knl, 2));
            KnlActual([fiMind(:); nFI + ccMind(:)], :) = Knl;

            tlistEmpty = repmat([allFIModes, allCCModes]', 1, 3);
            nEmpty = size(tlistEmpty, 1);
            obj.NLROM.tlist = [tlistActual; tlistEmpty];
            
            obj.NLROM.Knl = [KnlActual, zeros(nFI + nCC, nEmpty)];
            [obj.NLROM.N1, obj.NLROM.N2, obj.NLROM.Nlist] = Knl2Nash(obj.NLROM.Knl, obj.NLROM.tlist);
            obj.NLROM.dispErr = dispErr; obj.NLROM.forceErr = forceErr;
            obj.NLROM.fitTime = fitTime;
            
            % Optional outputs
            data.basis = basis; data.F = F; data.P = P; data.displacement = displacement;
            data.fitTime = fitTime; data.Kl = Kl;
            data.dispErr = dispErr; data.forceErr = forceErr;
            
        end
        
        function obj = genLoadCases(obj, basis, springs)
            % Develop cases for the model using the supplied
            % modal basis. Basis should be supplied in terms of the
            % fixed-interface and constraint mode amplitudes
            
            % Check the input lengths on "springs"
            if (length(springs.k) ~= length(springs.dof))
                error('Incompatible lengths on spring DOF and stiffnesses');
            end
            
            nM = length(obj.linModel.omegaN); nB = length(obj.linModel.connDofAbs);
            nC = size(basis, 2) - nM; nI = length(obj.linModel.internalDofAbs);
            % Convert CC basis into physical basis
            
            ccBasis = [ obj.linModel.psiIB; 
                          eye(nB, nB)]*basis((nM + 1):end, (nM + 1:end));
            loadDof = [obj.linModel.internalDofAbs; obj.linModel.connDofAbs];
            
            [dispTest, cf] = obj.abint.static(obj.linModel.bcDofAbs, ccBasis, loadDof, 'linear', springs);
       
            
            
            % Run preliminary load cases; linear [FI]
            phi = [obj.linModel.phiFI; zeros(nB, nM)];
            Kl = obj.linModel.K(loadDof, loadDof);
            Ml = obj.linModel.M(loadDof, loadDof);
            
            loadFI = bsxfun(@rdivide, obj.thick*Kl*phi, max(abs(phi)));
            loadFI((nI + 1):end, :) = 0;
            [dispFI, cf] = obj.abint.static(obj.linModel.bcDofAbs, loadFI, loadDof, 'nonlinear', springs);
            % This was all used for that "deflection estimate" stuff which
            % I don't like anymore
            
%             
%             maxFI = max(abs(dispFI))';
%             phiNl = diag(pinv(phi)*dispFI);
%             fprintf('FI Displacement Ratios; Initial:\n');
%             for i = 1:nM
%                 fprintf('\tMode %i: %.2f%% \n', i, (maxFI(i))/obj.thick);
%             end
%             % Now fit to a quadratic, cubic, or cubic + quadratic 1-mode
%             % model
%             frHat = diag(phi'*loadFI);
%             bEst = (diag(phi'*loadFI) - obj.linModel.omegaN(:).^2.*phiNl)./phiNl.^3;
%             
%             % Use these numbers to estimate modal force level to obtain +/-
%             % 15% of the linear response
%             % f = frHat/omega^2; c = b/omega^2
%             cubicEq = @(f, c) fzero(@(q) q + c*q^3 - f, 0);
%             deflect = zeros(nM, 1);
%             for i = 1:nM
%                 omega = obj.linModel.omegaN(i);
%                 c = bEst(i)/omega^2;
%                 startPoint = frHat(i)/omega^2;
%                 force = fzero(@(f) (f - cubicEq(f, c))^2/f^2 - 0.15^2, startPoint);
%                 deflect(i) = max(abs(phi(:, i)))*force;
%             end
            
            % Generate load cases using only + load cases for FI and CC
            % modes
            deflectCC = obj.thick*ones(size(ccBasis, 2), 1);
            freqCC = sqrt(diag(ccBasis'*Kl*ccBasis))/2/pi;
            loadCC = Kl*ccBasis;
            [dispCCLin, cf] = obj.abint.static(obj.linModel.bcDofAbs, loadCC, loadDof, 'linear', springs);
            maxCCLin = max(abs(dispCCLin));
            % The linear step above gives us an idea on how to scale the
            % loads to get deflections on the order of our base thickness
            loadCC = bsxfun(@times, obj.thick./maxCCLin, loadCC);

            % Now run a nonlinear step
            [dispCCNLin, cf] = obj.abint.static(obj.linModel.bcDofAbs, loadCC, loadDof, 'nonlinear', springs);

            disp = [dispFI, dispCCNLin]; F = [loadFI, loadCC];
            basis = [phi, ccBasis]; 
            freq = [obj.linModel.omegaN/2/pi; freqCC];
            flist = repmat((1:(nM + nC))', 1, 3); 
            
            % Save results if desired
            save('panelROM', 'disp', 'F', 'basis', 'freq', 'flist');
            

            % Single-mode or FI-only ROM Generation
            singleMode = false;
            if singleMode
                fitMode = 1; % Can do the FI modes also if setting a range on "fitMode"
                disp = disp(:, fitMode); F = F(:, fitMode); basis = basis(:, fitMode);
                freq = freq(fitMode); flist = flist(1:length(fitMode), :);
                
                [obj.phi_m, obj.nums_m] = expansion(disp, basis);
                [Knl, tlist] = noregress(freq, basis, F, ...
                    disp, 1, 1, flist);
                
                n = length(fitMode); N = nM + nC;
                tlist = [flist; repmat(((n + 1):N)', 1, 3)];
                Knl = [Knl, zeros(n, N - n); zeros(N - n, N)];
                
                util.F = F; util.disp = disp; util.basis = basis;
            end
            
                
            doubleBasis = false;
            if doubleBasis
                
                % Use FI mode 1 and CC mode 2; compute separately
                [U, S, V] = svd(basis, 'econ');
                
                obj.abint.static([], U(:, 1:end), loadDof, 'linear', springs);
                
                [obj.phi_m, obj.nums_m] = expansion(disp, basis);
                keyboard
                % FI Mode
                [Knl1, tlist1] = noregress(freq(1), basis(:, 1), F(:, 1), ...
                        disp(:, 1), 0, 1, []); 
                    
                % CC Mode - compute at a higher load level
                [dispCCNLin, cf] = obj.abint.static([], 2.5*loadCC(:, 1), loadDof, 'nonlinear', springs);
                [Knl2, tlist2] = noregress(freq(6), basis(:, 6), dispCCNLin, ...
                        disp(:, 6), 1, 1, []); 
                    
                Knl = diag([Knl1, zeros(1, 4), Knl2(1), zeros(1, 7)]);
                tlist = flist; tlist(6, 3) = 0;
                util.F = F; util.disp = disp; util.basis = basis;
            end
            
            diagBasis = true;
            if diagBasis
                
                % For this one, use all modes to make diagonal basis; using
                % both cubic and quadratic fits (NOPE, only cubic)
%                 dispNeg = obj.abint.static([], -F, loadDof, 'nonlinear', springs);
%                 keyboard
%                 F = [F, -F]; disp = [disp, dispNeg]; flist = [flist; flist];
%                 flist((end/2 + 1):end, 3) = 0;
%                 keyboard
                [obj.phi_m, obj.nums_m] = expansion(disp, basis);
                [Knl, tlist] = noregress(freq, basis, F, ...
                    disp, 1, 1, flist);
                
%                 Knl = zeros(nM + nC, nM + nC);
%                 tlist = flist;
%                 tlist((nM + 1):(nM + nC), 3) = 0;
%                 
%                 for i = 1:nM
%                     Knl(i, i) = noregress(freq(i), basis(:, i), F(:, i), ...
%                         disp(:, i), 0, 1, [1, 1, 1]);
%                 end
%                 for i = (nM + 1):(nM + nC)
%                     Knl(i, i) = noregress(freq(i), basis(:, i), F(:, i), ...
%                         disp(:, i), 1, 0, [1, 1, 0]);
%                 end
                
                util.F = F; util.disp = disp; util.basis = basis;
            end
            
            
            coupledBasis = false;
            if coupledBasis
                deflectCC = obj.thick*ones(size(ccBasis, 2), 1);
                freqCC = sqrt(diag(ccBasis'*Kl*ccBasis));
                [util.P, util.F] = icloading([obj.linModel.omegaN/2/pi; freqCC], ...
                [phi, ccBasis], [deflect; deflectCC], Kl, [], 1, 1, 1, ...
                obj.options.rf, 0, obj.options.cf);
                        
                util.F((nI + 1):end, :) = 0;
                [util.disp, cfCheck] = obj.abint.static([], util.F(:, 1:end), loadDof, 'nonlinear', springs);
                keyboard
                [obj.phi_m, obj.nums_m] = expansion(util.disp, phi); %#ok
                [Knl,tlist] = noregress_constrained([obj.linModel.omegaN/2/pi; freqCC], [phi, ccBasis], cfCheck, ...
                    util.disp, obj.options.quads, obj.options.triple, []);
            end
            
            % Need a nice way to add zeros to Knl and tlist for omitted
            % modes
            
            [NLROM.N1, NLROM.N2, NLROM.Nlist]= Knl2Nash(Knl,tlist);
            NLROM.Knl = Knl;    
            NLROM.tlist = tlist;
            
            % Must modify Khat and Mhat to accomodate the included CC basis
            % functions
            NLROM.Khat = blkdiag(diag(obj.linModel.omegaN.^2), ...
                         ccBasis'*Kl*ccBasis);
            NLROM.Mhat = blkdiag(eye(length(obj.linModel.omegaN)), ...
                         ccBasis'*Ml*ccBasis);    
                
            obj.NLROM = NLROM;
            obj.util = util;

        end
        
        function obj = zeroNLROM(obj, ncc)
            % Construct an NLROM with no nonlinearity
            N = length(obj.linModel.omegaN) + ncc;
            Knl = zeros(N);
            tlist = repmat((1:N)', 1, 3);
            [NLROM.N1, NLROM.N2, NLROM.Nlist]= Knl2Nash(Knl,tlist);
            NLROM.Knl = Knl;    
            NLROM.tlist = tlist;
            obj.NLROM = NLROM;
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
        
        function obj = fitCoefficients(obj)
            error('Not implemented');
            % Build NLROM Equation of motion
            % 1a: Hollkamps code
            % *Note that we need to use constrained with triple cubic terms in order
            % for N1(q) and N2(q,q) to be symmetric.

            F = obj.util.F; disp = obj.util.disp;
            quads = obj.options.quads; triple = obj.options.triple;
            constr = obj.options.constr;
            %keyboard
            if constr == 1; % constrained fit
                [Knl,tlist] = noregress_constrained(obj.freqsReduced, obj.phiReduced, F, ...
                    disp, quads, triple, []);
            elseif constr == 0;    % unconstrained fit
                [Knl,tlist] = noregress(obj.freqsReduced, obj.phiReduced, F, ...
                    disp, quads,triple, []);
            end
            
            % Pack NLROM for saving
            obj.NLROM = ICEROM.packNLROM(Knl, tlist, obj.mind, obj.x_rat, obj.model, ...
                obj.phiReduced, obj.linModel.DOF, obj.linModel.fn, ...
                obj.linModel.nodes, obj.linModel.eles, obj.phi_m, obj.nums_m,...
                obj.NNMoptions);
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
        
        
        function sys = defsys(obj, mode)
            %    sys = defsys(obj)
            %
            % Create the "sys" object expected by NNMcont using the
            % included "NNMoptions" structure
            
            % Linear System:
            sys.Klin = []; % obj.NLROM.Khat;
            sys.Mlin = []; %obj.NLROM.Mhat;
            
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
            sys.filename = sprintf('CMS_%s_NNM-%i', obj.abint.modelname, mode);
            sys.mode = mode;
            
        
        end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% "SET" METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = setLinModel(obj, linModel)
            obj.linModel = linModel;
        end
        
        function obj = setNLROM(obj, NLROM)
            obj.NLROM = NLROM;
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UTILITY METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function save(obj)
            cbice = obj;
            filename = [obj.abint.modelname];
            save(fullfile(obj.filepath, filename), 'cbice');
        end
        
        function plot(obj, amplitudes)
            %    plot(amplitudes)
            %
            % Plot the part using the given modal amplitudes. These should
            % be supplied as a column vector
            %
            % [qi; xb]
            %
            % Where qi are the fixed interface modes and xb are the
            % boundary DOF. Critically, the modal amplitudes must be given
            % in PHYSICAL COORDINATES to allow uniformity of plotting
            % between different parts in a substructure. If the modal
            % amplitude vector does not have the correct size, the function
            % will issue an error.
            
            % Check lengths on the input coordinates
            if ((size(amplitudes, 1) ~= (length(obj.linModel.connDofAbs) ...
                    + length(obj.linModel.omegaN))) || (size(amplitudes, 2) ~= 1))
                error('Wrong size on "amplitudes" in CBICE plot');
            end
            
            % Build X vector of displacements
            qi = amplitudes(1:length(obj.linModel.omegaN));
            xb = amplitudes((length(obj.linModel.omegaN) + 1):end);
            
            X = [obj.linModel.phiFI*qi + obj.linModel.psiIB*xb;
                 xb; zeros(size(obj.linModel.bcDofAbs))];
            
            obj.abint.deform(X(obj.linModel.dofOrder), 1);
            
        end
            
        function plotMode(obj, mode, scale)
            nM = length(obj.linModel.omegaN);
            nB = length(obj.linModel.connDofAbs);
            if (mode <= nM)
                amps = zeros(nM + nB, 1); 
                amps(mode) = scale;
            else
                Z = zeros(nB, 1); Z(mode - nM) = 1;
                amps = [zeros(nM, 1);
                        Z(:)]; 
            end
            obj.plot(amps);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% OPERATOR OVERLOADS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function disp(obj)
            % Display the object in a readable format
            fprintf('\nCBICE Object\n');
            fprintf('\t%s\n', obj.abint.modelname);
            if isempty(obj.NLROM)
                fprintf('\tNo ROM fit exists\n');
            else
                fprintf('\tROM fit exists');
            end
            byteSize = whos('obj'); byteSize = byteSize.bytes;
            fprintf('\tSize of %.3f mb\n\n', byteSize/1024^2);
        end
        
    end
    
    methods(Static)
        
        function linModel = packLinearModel(phi, fn, DOF, N, Nmds, Nnodes, M, ...
                nodes, eles)
            
            linModel.phi = phi;
            linModel.fn = fn;
            linModel.DOF = DOF;
            linModel.N = N;
            linModel.Nmds = Nmds;
            linModel.Nnodes = Nnodes;
            linModel.M = M;
            linModel.nodes = nodes;
            linModel.eles = eles;
            linModel.phiTM = phi'*M;
        end
        
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
                opt.TFS=100; % Time steps per period integration
            end
            if ~isfield(opt,'PrecS'),
                opt.PrecS = 1E-6; % 
            end
            if ~isfield(opt,'itoptS'),
                opt.itoptS = 3; % No idea what this does
            end
            if ~isfield(opt,'itmaxS'),
                opt.itmaxS = 10; % "h_max???"
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
            if ~isfield(opt,'alpha_f'),
                opt.alpha_f = 0.01; % Don't know what this does
            end
            
            if ~isfield(opt,'h'); opt.h = 1E-5; end
            if ~isfield(opt, 'hMax'); opt.hMax = 1E-2; end
            if ~isfield(opt, 'hMin'); opt.hMin = 1E-9; end
            if ~isfield(opt, 'betamax'); opt.betamax = 90; end
            if ~isfield(opt, 'betamin'); opt.betamin = 0; end
            if ~isfield(opt, 'wstop'); opt.wstop = []; end
            if ~isfield(opt, 'shootingPeriod'); opt.shootingPeriod = 1; end
            if strcmpi(opt.shootingPeriod, 'HALF'); opt.shootingPeriod = 1; end
            if strcmpi(opt.shootingPeriod, 'FULL'); opt.shootingPeriod = 2; end
            if (opt.shootingPeriod ~= 1) && (opt.shootingPeriod ~= 2)
                error('Incorrect specification for shootingPeriod')
            end
        end
                
        function NLROM = packNLROM(Knl, tlist, mind, modelname, ...
                phiReduced, fn, phi_m, nums_m, NNMoptions, linDispLC, thick)
            
            N = length(mind);
            [NLROM.N1, NLROM.N2, NLROM.Nlist]= Knl2Nash(Knl,tlist);
            NLROM.Knl = Knl;    
            NLROM.tlist = tlist;
            NLROM.mind = mind;  
            NLROM.modelname = modelname;
            NLROM.phi = phiReduced;     % This is phi(:, mind);
            NLROM.fn = fn(mind);    
            NLROM.phi_m = phi_m;
            NLROM.nums_m = nums_m;
            NLROM.NNMoptions = NNMoptions;
            NLROM.Khat = diag((2*pi*fn(mind)).^2);
            NLROM.Mhat = eye(N);
            NLROM.linDispLC = linDispLC;
            NLROM.thick = thick;
        end
    end
end




