% CBSS(filepath, filename) 
% 
% Craig-Bampton substructuring routine for CBICE models. 
%
% 
 

classdef CBSS

    properties(GetAccess = 'public', SetAccess = 'public')
        
        % Administrative
        matDir;     % Starting MATLAB directory for any function call (always returns here)
        filepath;   % Filepath for save of .mat file
        filename;   % Filename for save of .mat file
        
        % Connectivity cell array
        connect = cell(0);
        
        % Structure array of parts
        parts = struct('part', [], 'refNode', [], 'location', [], ...
                'connectors', [], 'nBoundaryDOF', [], 'nModalDOF', [], ...
                'amplitudes', [], 'internalModes', []);
        npart = 0;
        
        % DOF Handling
        nBoundaryDOF = 0;
        nModalDOF = 0;
        nConstr = 0;
        boundaryRelations = containers.Map('KeyType', 'char', 'ValueType', 'double');
        B; % Signed Boolean B-matrix
        L; % Unsigned Boolean assembly matrix
        
        % CB Solutions
        Msys; % CB system mass matrix
        Ksys; % CB System stiffness matrix
        phi;  % CB Modal matrix
        fn;   % CB natural frequencies [Hz]
        
        % Full order model handling
        fullOrder = struct('nodes', [], ... % Full order nodes
                           'dof',   [], ... % Full order DOF ordering
                           'phi',   [], ... % Full order modeshapes
                           'fn',    []);    % Full order frequencies
        phiFull; % Full system phi matrix for comparison with FEA
        mac;     % MAC of the system
    end
    
    
    methods
        % Constructor
        function obj = CBSS(filepath, filename)
            
            obj.matDir = cd;
            
            % Assign working directory, filepath, and filename
            if ~ischar(filepath); error('Invalid filepath'); end
            if ~ischar(filename); error('Invalid filename'); end
            
            % Perform input checks. This will fail if directories do not
            % exist, it will not create them.
            try 
                cd(filepath); obj.filepath = cd; cd(obj.matDir);
            catch err
                error('%s does not exist', filepath);
            end
            
            obj.filename = filename;
            
        end
        
        function [Ki, ki] = interfaceReduction(obj, partNumber)
            %    [Ki, ki] = interfaceReduction(obj, partNumber)
            %
            % Find the effective interface stiffness as seen by
            % "partNumber."
            
            % Build a partition of the boundary DOF only for the model
            boundaryDOF = [];
            currDOF = 1;
            for i = 1:obj.npart
                currDOF = currDOF + obj.parts(i).nModalDOF;
                endDOF = currDOF + obj.parts(i).nBoundaryDOF - 1;
                boundaryDOF = [boundaryDOF, currDOF:endDOF];
                currDOF = endDOF + 1;
            end
            if (length(boundaryDOF) ~= obj.nBoundaryDOF)
                error('Wrong size of Boundary DOF Vector; check code');
            end
            
            % Get global boundary DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(partNumber - 1)
                endDOF = currDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = partNumber;
            nB = obj.parts(i).nBoundaryDOF;
            
            % Obtain the appropriate set of interface DOF to be partitioned
            bSet = currDOF:(currDOF + nB - 1); % Boundary set
            iSet = [1:(currDOF - 1), (currDOF + nB):obj.nBoundaryDOF];
            
            % Make a "reduced" L matrix
            Lb = obj.L(boundaryDOF, sum(obj.L) > 1);
            Lr = Lb; Lr(bSet, :) = 0;
%             Br = obj.B;(:, iSet);
%             Lr2 = null(Br, 'r');
            % Assemble new global stiffness matrix
            Kr = Lr'*full(obj.Ksys(boundaryDOF, boundaryDOF))*Lr;
            
            % Find the assembled interface sets
            iSet = 1:size(Lb, 2);
            for i = 1:length(bSet)
                ind = find(Lb(bSet(i), :), 1);
                bSet(i) = ind;
                iSet(ind) = 0;
            end
            iSet(iSet == 0) = [];
            
            % Partition global stiffness matrix
            Krii = Kr(iSet, iSet);
            Krib = Kr(iSet, bSet);
            Krbi = Krib.';
            Krbb = Kr(bSet, bSet);
            
            Ki = Krbb - Krbi*(Krii\Krib);
            
            % Get fully reduced stiffness at each DOF
            ki = zeros(nB, 1);
            
            for i = 1:nB
                bSet = i;
                iSet = [1:(i - 1), (i + 1):nB];
                Kii = Ki(iSet, iSet);
                Kib = Ki(iSet, bSet);
                Kbi = Kib.';
                Kbb = Ki(bSet, bSet);
            
                ki(i) = Kbb - Kbi*(Kii\Kib);
            end
        end
        
        function obj = removeConnection(obj, part1, part2)
            %    obj = removeConnection(part1, part2)
            %
            % For each element of part1, remove connections to all elements
            % of part2
            for i = part1
                for j = part2
%                     keyboard
                    if (i == j); continue; end
                    if isempty(obj.connect{i, j}); continue; end
                    obj.connect{i, j} = []; obj.connect{j, i} = [];
                end
            end
        end
        
        function obj = buildCBSystem(obj, nModes)
            %    obj = buildCBSystem(obj, nModes)
            %
            % Assembly the system using the primal CB method. 
            
            obj = obj.buildBL();
            Mlocal = obj.L'*full(obj.Msys)*obj.L;
            Klocal = obj.L'*full(obj.Ksys)*obj.L;
            [phiLocal, lambda] = eigs(Klocal, Mlocal, nModes, 'SM');
            obj.phi = obj.L*real(phiLocal);
            obj.fn = real(sqrt(diag(lambda))/2/pi);
            
            % Assign each part its set of modal amplitudes
            obj = assignAmplitudes(obj);
            
            
        end
        
        function obj = assignAmplitudes(obj)
            %    obj = assignAmplitudes()
            %
            % Assign each part its set of modal amplitudes. Requires an
            % existing phi matrix.
            currDOF = 1;
            for i = 1:obj.npart
                nDof = obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                ampVec = currDOF:(currDOF + nDof - 1);
                try
                    obj.parts(i).amplitudes = obj.phi(ampVec, :);
                catch err
                    keyboard
                end
                currDOF = currDOF + nDof;
            end
        end
        
        function obj = buildBL(obj)
            % Perform preliminary system substructuring for both the CB and 
            % the CC substructuring techniques. Constructs a
            % signed binary B matrix of DOF constraints then uses SVD to
            % find the null space of B, yielding the boolean assembly
            % matrix L. L is then used as a transformation from a block
            % diagonal matrix of each CB component matrix to the final
            % system CB matrix.
            
            % 1) Build the B matrix. First, we need a list of the system
            % DOF number of each boundary DOF in the system. This also
            % determines the total number of boundary DOF in the system.
            % (We don't know otherwise).
            % Reset "boundaryRelations" map
            obj.boundaryRelations = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.nBoundaryDOF = 0;
            currDOF = 0; % Running count of "current" DOF in loop
            for i = 1:obj.npart % Run through each part
                currDOF = currDOF + obj.parts(i).nModalDOF;
                connSet = []; % Need to get a UNIQUE list of DOF that are actually attached
                obj.parts(i).nBoundaryDOF = 0;
                for j = 1:obj.npart % Run through each connection set
                    % Need to get a UNIQUE list of nodes that are actually
                    % being attached to.
                        if (i == j); continue; end % Skip boundary sets
                        setA = obj.connect{i, j};
                        if isempty(setA); continue; end % Skip empty sets
                        % Change to node.dof format (note transpose)
                        
                        setA = bsxfun(@plus, setA(:, 1), bsxfun(@times, setA(:, 2:7),(1:6)/10))';
                        setA = setA(:);
                        setA((setA - floor(setA) == 0)) = [];
                        
                        connSet = [connSet; setA];
                        % Also handle the connectivity, noting that a node
                        % may have more than 1 attaching DOF

                end

                connSet = unique(connSet, 'stable'); % Only unique DOF
                nBDof = length(connSet);
                % Now, assign to each of these DOF an absolute system DOF
                % number
                for j = 1:length(connSet);
                    % DOF on different parts need unique ID; use a string
                    dofStr = sprintf('part(%i).%.1f', i, connSet(j));
                    obj.boundaryRelations(dofStr) = currDOF + j; 
                end
                currDOF = currDOF + nBDof;
                obj.parts(i).nBoundaryDOF = obj.parts(i).nBoundaryDOF + nBDof;
                obj.nBoundaryDOF = obj.nBoundaryDOF + nBDof;
            end
            
            

            % Now we know the absolute DOF number of all boundary DOF in
            % the model. We have to do the same thing, but this time go
            % through and keep track of DOF connectivity.
            currConstr = 1; % Current constraint number
            obj.B = sparse(1, obj.nModalDOF + obj.nBoundaryDOF);
            for i = 1:obj.npart % Run through each part
                for j = (i + 1):obj.npart % Run through each connection set
                    
                    % Now we can go through again, but since we actually
                    % know where the DOF are in the model, we can build the
                    % B matrix.
                    setA = obj.connect{i, j};
                    setB = obj.connect{j, i};
                    if isempty(setA); continue; end

                    % Change to node.dof format
                    setA = bsxfun(@plus, setA(:, 1), bsxfun(@times, setA(:, 2:7),(1:6)/10))';
                    setA = setA(:);
                    setA((setA - floor(setA) == 0)) = [];
                    setB = bsxfun(@plus, setB(:, 1), bsxfun(@times, setB(:, 2:7),(1:6)/10))';
                    setB = setB(:);
                    setB((setB - floor(setB) == 0)) = [];

                    if (length(setA) ~= length(setB))
                        error('Incompatible connection sets for parts %i and %j', i, j);
                    end
                    
                    % NOW we can build the B-matrix
                    for k = 1:length(setA)
                        dofStrI = sprintf('part(%i).%.1f', i, setA(k));
                        dofStrJ = sprintf('part(%i).%.1f', j, setB(k));
                        try
                        constrI = obj.boundaryRelations(dofStrI);
                        constrJ = obj.boundaryRelations(dofStrJ);
                        catch err
                            keyboard
                        end
                        obj.B(currConstr, constrI) = 1;
                        obj.B(currConstr, constrJ) = -1;
                        currConstr = currConstr + 1;
                    end
                end
            end
            assert(size(obj.B, 2) == (obj.nBoundaryDOF + obj.nModalDOF), ...
                'Incompatible size of B Matrix');
            % Get total number of constraints
            obj.nConstr = currConstr - 1;
            % Finally, the B matrix is here! Now use build the L matrix.
%             spLeft = spspaces(obj.B.', 1);
%             obj.L = spLeft{1}(spLeft{3}, :).';
            warning('Using MATLAB null due to out-of-order results in spspaces');
            obj.L = null(obj.B, 'r');
            
            % Finally, assemble the full system and transform it
            obj.Msys = sparse([]); obj.Ksys = sparse([]); 
            for i = 1:obj.npart
                linModel = obj.parts(i).part.linModel;
                Mcb = [eye(length(linModel.omegaN)), linModel.MkbHat;
                       linModel.MkbHat', linModel.MbbHat];
                Z = zeros(size(linModel.MkbHat));
                Kcb = [diag(linModel.omegaN.^2), Z; 
                       Z', linModel.KbbHat];
                obj.Msys = blkdiag(obj.Msys, Mcb);
                obj.Ksys = blkdiag(obj.Ksys, Kcb);
            end
        end
            
        
        function obj = buildCBModels(obj, method)
            %    obj = buildCBModels(obj)
            %
            % Build all CB models of the components. Make sure to call
            % "makeConnections" first; no checking on the attachments is
            % performed.
            %
            
            if (nargin < 2); method = 'abaqus'; end
            
            % Run through each part and make its CB model
            obj.nModalDOF = 0;
            for i = 1:obj.npart
                % This call will assign:
                % omegaN
                % MkbHat
                % MbbHat
                % KbbHat
                %
                % To the structure obj.parts(i).part.linModel, a field of
                % the CBICE object. 
                goFlag = 0;
                obj.parts(i).nModalDOF = 0;
                % Make the connection set
                connSet = [];
                for j = 1:obj.npart
                    if (i == j); continue; end
                    connSet = [connSet; obj.connect{i, j}];
                end
                % Sort to unique connections; stable sort
                [~, ia, ~] = unique(connSet(:, 1), 'stable');
                connSet = connSet(ia, :);
                
                % See if we can use models from any existing parts
                for j = 1:(i - 1)
                    if ~(obj.parts(i).part.abint.isequal(obj.parts(j).part.abint));
                        continue;
                    end
                    compareSet = [];
                    for k = 1:obj.npart
                        if (j == k); continue; end
                        compareSet = [compareSet; obj.connect{j, k}];
                    end
                    % Sort to unique connections; stable sort
                    [~, ja, ~] = unique(compareSet(:, 1), 'stable');
                    compareSet = compareSet(ja, :);
                    if ~(isequal(size(compareSet), size(connSet))); continue; end
                    if (norm(compareSet - connSet) == 0) % If comparison is true
                        obj.parts(i).part = obj.parts(i).part.setLinModel( ...
                            obj.parts(j).part.linModel);
                        goFlag = 1;
                        break
                    end
                end
                   
                if (goFlag); 
                    obj.nModalDOF = obj.nModalDOF + obj.parts(j).nModalDOF;
                    obj.parts(i).nModalDOF = obj.parts(j).nModalDOF;
                    continue; 
                end
                
                % Call the CB routine
                [obj.parts(i).part, omegaN] = obj.parts(i).part.makeCB(...
                    obj.parts(i).internalModes, obj.connect{i, i}, connSet, ...
                    method);
                
                % Add number of modes to total number of modal DOF in SS
                obj.nModalDOF = obj.nModalDOF + length(omegaN);
                obj.parts(i).nModalDOF = length(omegaN);
            end
            
        end
        
        function obj = makeConnections(obj, tol)
            %    obj = makeConnections(obj, tol)
            %
            % Find close nodes for eligible connections between each part.
            %
            
            if (nargin < 2); tol = 1E-3; end
            obj.nBoundaryDOF = 0; % Reset this counter
            for i = 1:obj.npart; obj.parts(i).nBoundaryDOF = 0; end % And this one
            
            % Run through and make connections between all parts. 
            for i = 1:obj.npart % Part A connection index
                for j = (i + 1):obj.npart % Part B connection index
                    
                partA = obj.parts(i).part.abint; % Part A (ABINT)
                partB = obj.parts(j).part.abint; % Part B (ABINT)
                
                connSetA = obj.parts(i).connectors; % Eligible nodes on part A for connection
                ncandA = size(connSetA, 1); % Number of part A candidate nodes
                connSetB = obj.parts(j).connectors; % Elgibile nodes on part B for connection
                ncandB = size(connSetB, 1); % Number of part B candidate nodes
                
                
                % Get node locations and find distances between nodes
                distances = zeros(ncandA, ncandB);
                partANodeLocs = partA.nodes(partA.getAbsNodeNumber(connSetA(:, 1)), 2:4);
                partBNodeLocs = partB.nodes(partB.getAbsNodeNumber(connSetB(:, 1)), 2:4);

                for k = 1:ncandB
                    % Find the distance between all nodes in part A and the
                    % ith node in part B
                    positions = bsxfun(@minus, partANodeLocs, partBNodeLocs(k, :));
                    distances(:, k) = sqrt(dot(positions, positions, 2));
                end

                filteredNodes = (distances < tol); % Filter nodes below a certain tolerance

                % Find indices of connected nodes
                % I: Connected nodes of part 1
                % J: Connected nodes from part 2
                [I, J] = find(filteredNodes, max(ncandA, ncandB));
                
                % Verify DOF compatibility
                connSetA = connSetA(I, :); connSetB = connSetB(J, :);
                for k = 1:length(J)
                    % Assign only DOF that
                    % are active for both parts to the connection
                    dofVec = connSetA(k, 2:end).*connSetB(k, 2:end);
                    connSetA(k, 2:end) = dofVec;
                    connSetB(k, 2:end) = dofVec;
                end

                % Assign to "connect" cell array
                obj.connect{i, j} = connSetA;
                obj.connect{j, i} = connSetB;
                

                end
            end
        end
   
        
        % Add a part to the model
        function obj = addPart(obj, part, refNode, location, boundary, connectors, internalModes)
            %    obj = addPart(part, refNode, location, boundary, connectors, internalModes)
            %
            % Add a component (CBICE object) to the model.
            %
            % INPUTS
            %
            % part: CBICE model of the part to be added
            %
            % location: [3 x 1] Position of the refNode in the model
            %
            % refNode: Node number for use as reference for positioning
            %
            % boundary: [nBound x 7] Set of fixed boundary DOF for the
            %           model; leave empty for none
            %
            % connectors: [nConn x 7] Set of connector DOF for the model
            %
            % internalModes: Either the number of FI modes to retain, or
            %                the frequency range of FI modes to retain
            %
            
            %%% Input checks %%%
            
            % Part check
            if ~strcmpi(class(part), 'CBICE'); error('Require CBICE input for "part"'); end
            
            % RefNode check
            try 
                absRefNode = part.abint.getAbsNodeNumber(refNode); 
            catch err
                error('Node %i not found in %s', refNode, part.abint.modelname);
            end
            
            % Location check
            if (sum(size(location(:)) ~= [3, 1]) ~= 0)
                sz = size(location);
                error('"Location" has size [%i, %i], requires size [3, 1] or [1, 3]', ...
                    sz(1), sz(2));
            end
            
            % Boundary check
            if isempty(boundary)
                nBC = 0;
            elseif (size(boundary, 2) ~= 7); 
                error('Need 7 columns in boundary for part %s', part.abint.modelname); 
            else
                nBC = size(boundary, 1);
            end
            if (nBC > part.abint.n); error('More boundary nodes than model nodes'); end
            if nBC ~= 0
                try 
                    part.abint.getAbsNodeNumber(boundary(:, 1));
                catch err
                    error('Bad boundary nodeset for part %s', part.abint.modelname);
                end
            end
            
            % Connector check
            if (size(connectors, 2) ~= 7); 
                error('Need 7 columns in connector for part %s', part.abint.modelname); 
            else
                nConn = size(connectors, 1);
            end
            if (nConn > part.abint.n); error('More connector nodes than model nodes for part %s', part.abint.modelname); end
            try 
                part.abint.getAbsNodeNumber(connectors(:, 1));
            catch err
                error('Bad connector nodeset for part %s', part.abint.modelname);
            end
            
            % internalModes check
            if (nargin < 7); internalModes = 10; end
            
            % Add part and location to internal structure array
            % Change location to use an absolute offset
            location = location(:)' - part.abint.nodes(absRefNode, 2:4);
            part.abint = part.abint.offsetNodes(location);
            obj.npart = obj.npart + 1;
            obj.parts(obj.npart) = struct('part', part, ...
                'refNode', refNode, 'location', location, ...
                'connectors', connectors, 'nBoundaryDOF', 0, ...
                'nModalDOF', 0, 'amplitudes', [], 'internalModes', internalModes);
            
            
            % Add boundary set to cell array 
            obj.connect{obj.npart, obj.npart} = boundary;
             
        end
        
        
        function obj = clearParts(obj)
            %    obj = clearParts(obj)
            %
            % Remove all parts from the model and reset. This KEEPS the full
            % order model information. (That's the point).

            % Connectivity cell array
            obj.connect = cell(0);

            % Structure array of parts
            obj.parts = struct('part', [], 'refNode', [], 'location', [], ...
                    'connectors', [], 'nBoundaryDOF', [], 'nModalDOF', [], ...
                    'amplitudes', [], 'internalModes', []);
            obj.npart = 0;

            % DOF Handling
            obj.nBoundaryDOF = 0;
            obj.nModalDOF = 0;
            obj.nConstr = 0;
            obj.boundaryRelations = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.B; % Signed Boolean B-matrix
            obj.L; % Unsigned Boolean assembly matrix

            % CB Solutions
            obj.Msys; % CB system mass matrix
            obj.Ksys; % CB System stiffness matrix
            obj.phi;  % CB Modal matrix
            obj.fn;   % CB natural frequencies [Hz]

            obj.phiFull = []; % Full system phi matrix for comparison with FEA
            obj.mac = [];     % MAC of the system
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% FULL ORDER COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function compareTable(obj, nModes, type, outFile)
            %    compareTable(obj, nModes, type, outFile)
            %
            % Write a comparison table between full-order and CB results.
            % Specify number of modes and "type" as "text" "html" or "tex".
            % "outFile," if provided, will save the table output to the
            % specified filepath.
            
            if nargin < 2; nModes = length(obj.fn); end
            if nargin < 3; type = 'text'; end
            if nargin < 4; 
                fid = 1;
            else
                fid = fopen(outFile, 'wt');
                if (fid == -1); error('Invalid File'); end
            end
            
            % Use the CB modes as the reference. For each mode, find the
            % "matching" MAC value as well as the 2nd-highest MAC value on
            % each column.
            [macVals, matchInds] = max(obj.mac, [], 1); % Column-wise maximum
            macHold = obj.mac; 
            for i = 1:length(matchInds); macHold(matchInds(i), i) = 0; end
            offVals = max(macHold, [], 1);
            fullFreqs = obj.fullOrder.fn(matchInds);
            
            freqCompare = @(i) [obj.fn(i), fullFreqs(i), (obj.fn(i) - fullFreqs(i))/ ...
                                fullFreqs(i)*100]; % Frequency error comparison
            macCompare = @(i) [macVals(i), offVals(i)];
                        
            % Write out tables
            if strcmpi(type, 'text');
                fprintf(fid, 'Mode\tMatch Mode\tCB Freq\tAb Freq\tError [%%]\tMAC\tMAC "Off-Diag"\n');
                for i = 1:nModes
                    fprintf(fid, '%i\t\t%i\t\t\t%.2f\t%.2f\t%.4f\t\t%.2f\t%.3f\n', i, ...
                        matchInds(i), freqCompare(i), macCompare(i));
                end
            end
            
            % LaTeX - not done yet
            % Writes "\begin{tabular} ..." through "\end{tabular}
            if strcmpi(type, 'tex');
                fprintf(fid, '\\begin{tabular}{ c || c | c | c | c | c | c}\n');
                fprintf(fid, 'Mode & Match Mode & SS Freq & Ab Freq & Error [\\%%] & MAC & Off-Mode MAC \\\\ \n');  
                fprintf(fid, '\\hline \\hline\n');
                for i = 1:nModes
%                     fprintf(fid, '\\hline\n');
                    fprintf(fid, '%i & %i & %.2f & %.2f & %.2f & %.2f & %.3f \\\\ \n', i, ...
                        matchInds(i), freqCompare(i), macCompare(i));
                end
                fprintf(fid, '\\end{tabular}\n');
            end
            
            % HTML
            if strcmpi(type, 'html');
                fprintf(fid, '<table>\n');
                fprintf(fid, '<tr><td>Mode</td><td>Match Mode</td><td>CB Freq</td><td>Ab Freq</td><td>Error [%%]</td><td>MAC</td><td>Off-Mode MAC</td></tr>\n');
                for i = 1:nModes
                    fprintf(fid, '<tr><td>%i</td><td>%i</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.3f</td></tr>\n', i, ...
                        matchInds(i), freqCompare(i), macCompare(i));
                end
                fprintf(fid, '</table>\n');
            end
            try 
                fclose(fid);
            catch err
            end
        end
        
        
        function fig = comparePlot(obj, maxErr, nModes, cmap)
            %    fig = comparePlot(obj, maxErr, nModes, cmap)
            %
            % Display a comparison summary image demonstrating both the MAC
            % and the frequency error with mode comparisons based on the
            % MAC.
            
            % Use the CB modes as the reference. For each mode, find the
            % "matching" MAC value as well as the 2nd-highest MAC value on
            % each column.
            if nargin < 2; maxErr = 2.5; end
            if nargin < 3; nModes = size(obj.phiFull, 2); end
            if nargin < 4; cmap = colormap('jet'); end
            close(gcf);
            
            [~, matchInds] = max(obj.mac, [], 1); % Column-wise maximum
            macHold = obj.mac; 
            for i = 1:length(matchInds); macHold(matchInds(i), i) = 0; end
            fullFreqs = obj.fullOrder.fn(matchInds);
            freqError = (obj.fn - fullFreqs)./fullFreqs*100;
            
            % Create figure/axes
            fig = figure('units', 'inches', 'position', [4, 2, 8, 8]);
            axMac = subplot(3, 1, 1:2); axFreq = subplot(3, 1, 3);
            
            
            % Handle MAC plot (this is copied from "MAC_plot" so that I can
            % put it on the axis where I want it)
            [a, b] = size(obj.mac); xs = [1:a]; ys = [1:b];
            nC = length(cmap) - 1;
            axes(axMac);
            for ii = 1:1:a;
                for jj = 1:1:b;
                    xcoord = xs(ii) + 0.5*[-1  1 1 -1];
                    ycoord = ys(jj) + 0.5*[-1 -1 1  1];
                    col_vec = cmap(round(obj.mac(ii,jj)*nC)+1,:);
                    patch(xcoord, ycoord, col_vec);
                end
            end
            set(axMac, 'XLim', [0.5, nModes + 0.5],'YLim', [0.5, nModes + 0.5],...
                'YTick', ys,'XTick', xs);
            ylabel('\bfFEA Modes');
            filestr = strrep(obj.filename, '_', '\_');
            titlestr = sprintf('%s; %i Modal DOF; %i Boundary DOF', filestr, ...
                obj.nModalDOF, obj.nBoundaryDOF);
            title(titlestr);
            colormap(cmap);
            colorbar('southoutside');
            xLim = xlim(axMac);
            
            
            % Handle frequency error plot
            axes(axFreq);
            for i = 1:nModes
                cmapInd = min(round(abs(freqError(i))/maxErr*nC) + 1, nC + 1);
                bar(axFreq, i, freqError(i), ...
                    'faceColor', cmap(cmapInd, :));
                hold(axFreq, 'on');
                % Add text denoting actual frequency
                offSet = -sign(freqError(i))*maxErr/10;
                text(i, offSet, sprintf('%.0f', obj.fullOrder.fn(matchInds(i))), ...
                    'HorizontalAlignment', 'center', 'fontSize', 10);
                % Add text denoting SS frequency
                loc1 = abs(freqError(i)) + maxErr/10;
                loc2 = maxErr - maxErr/10;
                loc3 = abs(freqError(i)) - maxErr/10; loc3 = loc3 + maxErr*(loc3 < 0);
                loc = sign(freqError(i))*min([loc1, loc2, loc3]);
                text(i, loc, sprintf('%.0f', obj.fn(i)), 'HorizontalAlignment', ...
                    'center', 'fontSize', 10);
            end
            set(axFreq, 'xlim', xLim, 'xTick', 1:nModes, 'yGrid', 'on', ...
                'ylim', maxErr*[-1, 1]);
            grid on
            ylabel(axFreq, '\bfFrequency Error [%]');
            xlabel(axFreq, '\bfCB Mode Number');
            title(axFreq, '\bfFrequency Errors; Annotations in Hz');
            
        end
        
        
        
        function [obj, mac] = makeMac(obj, tol)
            %    [obj, mac] = makeMac(obj, tol)
            %
            % Perform a MAC of the model and Abaqus part; obtain resulting
            % matrix. The mac is called as 
            %
            % mac = MAC(phiCB, phiFE)
            %
            
            if (nargin < 2); tol = 1E-3; end
            % (1) For each part, loop through and find the matching node
            % index. Assign to global phi matrix as we move along.
            nModes = size(obj.phi, 2); nFeDof = size(obj.fullOrder.phi, 1);
            obj.phiFull = zeros(nFeDof, nModes); % SS modal matrix to be filled
            feNodes = obj.fullOrder.nodes;
            
            for i = 1:obj.npart    
                part = obj.parts(i).part;
                fprintf('Matching part %i of %i (%i DOF):', i, obj.npart, ...
                    part.abint.N);
                % Make the part modal displacement matrix
                E = [part.linModel.phiFI, part.linModel.psiIB;
                     zeros(obj.parts(i).nBoundaryDOF, obj.parts(i).nModalDOF), ...
                     eye(obj.parts(i).nBoundaryDOF); ...
                     zeros(length(obj.connect{i, i})*6, length(obj.parts(i).amplitudes))];
                 
                % "phiPart" is the modal matrix for each component
                phiPart = E*obj.parts(i).amplitudes;
                absDofVector = [part.linModel.internalDofAbs;
                                part.linModel.connDofAbs;
                                part.linModel.bcDofAbs];
                            
                nodes = part.abint.nodes; 
                n = length(nodes);
                matchedNodes = 0; unmatchedNodes = 0;
                for j = 1:n % For each node
                    % Get distance from all FE nodes
                    positions = bsxfun(@minus, nodes(j, 2:4), feNodes(:, 2:4));
                    distances = sqrt(dot(positions, positions, 2));
                    [val, I] = min(distances); % I is the absolute FE node number
                    
                    if (val > tol); % If out of tolerance, continue
                        unmatchedNodes = unmatchedNodes + 1;
                        keyboard
                        continue;
                    end 
                    matchedNodes = matchedNodes + 1;
                    feNode = feNodes(I, 1);
                    
                    % J1:J2 is the absolute FE DOF range for this noe
                    J1 = find(floor(obj.fullOrder.dof) == feNode, 1, 'first');
                    J2 = J1 + find(floor(obj.fullOrder.dof((J1 + 1):end)) == feNode, 1, 'last');
                    
                    % Now get the DOF number in the Craig-Bampton ordering
                    cbNode = part.abint.getNodeNumber(j);
                    cbDof = part.abint.getAbsDofNumber(cbNode + 0.1);
                    cbDof = cbDof:(cbDof + (J2 - J1)); % List of DOF present in full model
                      
                    
                    % Now, find where each DOF is in "absDofVector;"
                    % assign that row from "phiPart" to the correct row in
                    % "phi"
                    for k = 1:length(cbDof)
                        try 
                            pullIdx = find(absDofVector == cbDof(k), 1, 'first');
                            obj.phiFull(J1 + k - 1, :) = phiPart(pullIdx, :);
                        catch err
                            keyboard
                            rethrow(err);
                        end
                    end
                end
                fprintf('%i Matched Nodes; %i Unmatched Nodes\n', ...
                    matchedNodes, unmatchedNodes);
            end 
            
            % Do the MAC
            
            obj.mac = MAC(obj.phiFull, obj.fullOrder.phi(:, 1:nModes));
            mac = obj.mac;
        end
            
        
        
        function obj = loadOdb(obj, odbPath, odbFile, pythonPath)
            %    obj = loadOdb(obj, odbPath, odbFile, pythonPath)
            %
            % Load an ODB corresponding to a modal step and store the
            % results for later comparison. 
            
            obj.matDir = cd; % Hold starting directory
            % Get absolute python path
            cd(pythonPath); pythonPath = cd; cd(obj.matDir);
           
            % First, get the nodes from the model
            cd(odbPath);
            execStr = sprintf('abaqus python "%s" "%s" model', ...
                fullfile(pythonPath, 'odb2mat_v2.py'), odbFile);
            status = system(execStr);
            loadName = sprintf('%sFromAbaqus', odbFile);
            
            try
                load(loadName);
                delete([loadName '.mat']);
                % Assign ODB properties to ABINT
                obj.fullOrder.nodes = FE.nodes;    
            catch err
                cd(obj.matDir);
                error('Python read error');
            end
            
            
            % Now get modes, freqs, and DOF
            execStr = sprintf('abaqus python "%s" "%s" shapes', ...
                fullfile(pythonPath, 'odb2mat_v2.py'), odbFile);
            status = system(execStr);
            if status ~= 0; cd(obj.matDir), error('Python read error'); end
            
            loadName = sprintf('%sFromAbaqus', odbFile);
            load(loadName);
            
            % Clean file directory and leave this place
            delete([loadName '.mat']);
            cd(obj.matDir);
            
            % The object is called "mode" in the workspace; has fields
            % "shape" and "freqs"
            obj.fullOrder.phi = mode.shapes;
            obj.fullOrder.fn = mode.freqs;
            obj.fullOrder.dof = mode.dof;
        
        end
        
        function plotmode(obj, mode, varargin)
            %    pltmode(obj, mode, varargin)
            % Plot the deformation of a given mode for all parts
            if (nargin > 2)
                scale = varargin{1};
            else
                scale = 1;
            end
            for i = 1:obj.npart
                % Call the CBICE "deform" function with the appropriate set
                % of modal amplitudes
                obj.parts(i).part.plot(scale*obj.parts(i).amplitudes(:, mode));
            end
            title(sprintf('\\bfMode %i; f_n = %.3f Hz', mode, obj.fn(mode)));
        end
            
        function deform(obj, amplitudes, varargin)
            %    pltmode(obj, mode, varargin)
            % Plot a given deformation for all parts
            
            if ((size(amplitudes, 1) ~= (obj.nBoundaryDOF + obj.nModalDOF)) || (size(amplitudes, 2) ~= 1))
                error('Wrong size on "amplitudes" in CBSS plot');
            end
            
            currDof = 1;
            for i = 1:obj.npart
                % Call the CBICE "deform" function with the appropriate set
                % of modal amplitudes
                endDof = currDof + obj.parts(i).nBoundaryDOF + obj.parts(i).nModalDOF;
                partAmp = amplitudes(currDof:endDof - 1);
%                 keyboard
                obj.parts(i).part.plot(partAmp);
                currDof = endDof;
            end
            
        end
        
        function plot(obj, varargin)
            % Overloaded call to plot function. Additional name/value plot
            % property pairs will be carried forward.
            
            for i = 1:obj.npart
                % Patch plot of each part
                plotPart = obj.parts(i).part.abint;
                
                plotPart.plot(varargin);
                
                % Highlight boundary set
                if ~isempty(obj.connect{i, i})
                    plotPart.highlightNodes(obj.connect{i, i}(:, 1), ...
                        'MarkerEdgeColor', 'k', ...
                        'marker', '^', ...
                        'linestyle', 'none', ...
                        'MarkerFaceColor', 0.2*[1, 1, 1], ...
                        'MarkerSize', 7);
                    set(gca, 'nextplot', 'add');
                end
                
                % Highlight potential connector set
                plotPart.highlightNodes(obj.parts(i).connectors(:, 1), ...
                    'MarkerEdgeColor', [0.1, 0.3, 0.1], ...
                    'MarkerFaceColor', [0.1, 0.6, 0.1], ...
                    'Marker', 'o', ...
                    'linestyle', 'None', ...
                    'MarkerSize', 4);
                set(gca, 'nextplot', 'add');
                
                % Highlight attached connector set
                for j = 1:obj.npart
                    if i == j; continue; end
                    if isempty(obj.connect{i, j}); continue; end
                        plotPart.highlightNodes(obj.connect{i, j}(:, 1), ...
                            'MarkerEdgeColor', [0.3, 0.1, 0.1], ...
                            'MarkerFaceColor', [0.6, 0.1, 0.1], ...
                            'Marker', 's', ...
                            'Linestyle', 'none', ...
                            'MarkerSize', 5);
                        set(gca, 'nextplot', 'add');
                end
            end
            
        end
        
        
        function save(obj)
            cbss = obj; %#ok
            save(fullfile(obj.filepath, obj.filename), 'cbss');
        end
        
        function disp(obj)
            % Display the object in a readable format
            fprintf('\nCBSS Object\n');
            if obj.npart == 0;
                fprintf('No included parts\n');
            else
                fprintf('\t%s\n', obj.filename);
                fprintf('\t%i parts\n', obj.npart);
                fprintf('\t%i Modal DOF; %i Boundary DOF\n', ...
                    obj.nModalDOF, obj.nBoundaryDOF);
            end
            
            byteSize = whos('obj'); byteSize = byteSize.bytes;
            fprintf('\tSize of %.3f mb\n\n', byteSize/1024^2);
        end
        
        function summaryTable(obj, parts, type)
            % Display a readable table with a summary of desired parts
            % along with the total model
            if nargin < 2; parts = 1:obj.npart; end
            if nargin < 3; type = 'text'; end
            fid = 1; % Modify later to write files
            
            % Get total DOF count
            dofCount = 0;
            for i = 1:obj.npart
                dofCount = dofCount + obj.parts(i).part.abint.N;
            end
            
            if strcmpi(type, 'text');
                fprintf(fid, 'Part\t\tFreq. Range\tFI Modes\tBoundary DOF\n');
                % Part info
                for i = parts
                    if (length(obj.parts(i).internalModes) == 1)
                        fRange = '--';
                    else
                        fRange = sprintf('[%i, %i]', obj.parts(i).internalModes(1), ...
                            obj.parts(i).internalModes(2));
                    end
                    fprintf(fid, '%s\t%s\t%i\t%i\n', obj.parts(i).part.abint.modelname, ...
                        fRange, obj.parts(i).nModalDOF, obj.parts(i).nBoundaryDOF);
                end
                % Overall Model info
                fprintf(fid, 'Nconstraints: %i\t\tN Total DOF: %i\n',...
                    size(obj.B, 1), size(obj.B, 2));
                fprintf(fid, 'Physical DOF: %i\t\tAssembled DOF: %i\n', ...
                    dofCount, size(obj.L, 2));
            end
            
            if strcmpi(type, 'html');
                fprintf(fid, '<table><tr><td><b>Part</b></td><td><b>Freq. Range</b></td><td><b>FI Modes</b></td><td><b>Boundary DOF</b></td></tr>\n');
                % Part info
                for i = parts
                    if (length(obj.parts(i).internalModes) == 1)
                        fRange = '--';
                    else
                        fRange = sprintf('[%i, %i]', obj.parts(i).internalModes(1), ...
                            obj.parts(i).internalModes(2));
                    end
                    fprintf(fid, '<tr><td>%s</td><td>%s</td><td>%i</td><td>%i</td></tr>\n', obj.parts(i).part.abint.modelname, ...
                        fRange, obj.parts(i).nModalDOF, obj.parts(i).nBoundaryDOF);
                end
                % Overall Model info
                fprintf(fid, '<tr><td><b>Nconstraints</b></td><td>%i</td><td><b>N Total DOF</b></td><td>%i</td>\n', ...
                    size(obj.B, 1), size(obj.B, 2));
                fprintf(fid, '<tr><td><b>Physical DOF</b></td><td>%i</td><td><b>N Assembled DOF</b></td><td>%i</td></tr></table>\n', ...
                    dofCount, size(obj.L, 2));
            end
        end
            
    end
end