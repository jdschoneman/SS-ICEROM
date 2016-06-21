% ABINT(workingDir, filepath, filename, pythonpath) 
% Abaqus interface 
% class. Used to store Abaqus models as MATLAB data for
% rapid manipulation. Can write model input files and step 
% input files for selected analysis procedures. Reads 
% results and stores as necessary.
%
% Dependencies: odb2mat_v2.py in "pythonpath"
 

classdef ABINT

    properties(GetAccess = 'public', SetAccess = 'public')
        
        % Administrative
        matDir;     % Starting MATLAB directory for any function call (always returns here)
        workingDir; % Working directory for computation
        workingDirRel; % Relative (original) specification of working directory
        filepath;   % Filepath for save of .mat file
        filepathRel; % Relative (original) specification of file save path
        filename;   % Filename for save of .mat file
        pythonpath; % Path for odb2mat_v2 python file
        pythonpathRel; % Relative (original) specification of python path
        
        modelname;  % Abaqus model name (same as .INP filename)
        updateTime; % Time Abaqus model was read
        
        % Properties obtained from .odb file
        n;        % Number of nodes
        N;        % Number of DOF in model
        M;        % Number of elements
       
        nodes;     % [n x 4] matrix of element numbers and locations
        absNodes;  % Container.map relating node numbers to absolute node numbers
        elements;  % [M x 10] matrix of element connectivity
        nodeSets;  % Structure of model node sets
        eleSets;   % Structure of model element sets
        materials; % Structure of any material models
        
        % DOF properties
        dof;       % [M x 1] vector of dof in node.dof format
        absDof;    % Container.map relating node.dof to absolute dof number
        
        % Properties obtained from .inp file
        sections = containers.Map('KeyType', 'char', 'ValueType', 'any');  % Section assigments
        ties = containers.Map('KeyType', 'char', 'ValueType', 'any'); % Tie assignments
        surfaces = containers.Map('KeyType', 'char', 'ValueType', 'any'); % Surface assignments
        
        % Properties assigned manually
        boundaries; % Boundary condition vector
        feLoad;     % Load case vector [N x 1] for static or [N x t] for time
        
        % Stored properties
        response;    % Finite element response
        h;           % Patch handle
        lengthScale; % Length scale for plotting
    end
    
    
    methods
        % Constructor
        function obj = ABINT(workingDir, filepath, pythonpath)
            % Where does the help go?
            obj.matDir = cd;
            % Assign working directory, filepath, and filename
            if ~ischar(workingDir); error('Invalid workingDir'); end
            if ~ischar(filepath); error('Invalid filepath'); end
            if ~ischar(pythonpath); error('Invalid pythonpath'); end
            
            % Perform input checks. This will fail if directories do not
            % exist, it will not create them.
            try 
                cd(workingDir); obj.workingDir = cd; cd(obj.matDir);
                obj.workingDirRel = workingDir;
            catch err
                error('Working directory specification %s does not exist', workingDir);
            end
            try 
                cd(filepath); obj.filepath = cd; cd(obj.matDir);
                obj.filepathRel = filepath;
            catch err
                error('Filepath designation %s does not exist', filepath);
            end
            try 
                cd(pythonpath); obj.pythonpath = cd; cd(obj.matDir);
                obj.pythonpathRel = pythonpath;
            catch err
                error('Python path %s does not exist', pythonpath);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% "GET" METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function nodeNum = getNodeNumber(obj, absNodeNumber)
            %    nodeNum = getNodeNumber(absNodeNumber)
            %
            % Given absolute node numbers for this part, return actual node
            % numbers.
            if (size(absNodeNumber, 2) > 1); error('Wrong size on absNodeNumber'); end
            nodeNum = obj.nodes(absNodeNumber, 1);
        end
        
        function absNodeNum = getAbsNodeNumber(obj, nodeNumber)
            %    absNodeNum = getAbsNodeNumber(nodeNumber)
            %
            % Given actual node numbers for this part, return absolute node
            % numbers.
            if (size(nodeNumber, 2) > 1); error('Wrong size on nodeNumber');end
            absNodeNum = zeros(size(nodeNumber, 1), 1);
            for i = 1:length(absNodeNum)
                absNodeNum(i) = obj.absNodes(nodeNumber(i));
            end
        end
        
        function dofNum = getDofNumber(obj, absDofNumber)
            %    dofNum = getDofNumber(absDofNumber)
            %
            % Given absolute DOF numbers for this part, return actual DOF
            % numbers.
            if (size(absDofNumber, 2) > 1); error('Wrong size on absDofNumber'); end
            dofNum = obj.dof(absDofNumber);
        end
        
        function absDofNum = getAbsDofNumber(obj, dofNumber)
            %    absDofNum = getAbsDofNumber(dofNumber)
            %
            % Given actual DOF numbers for this part, return absolute DOF
            % numbers.
            if (size(dofNumber, 2) > 1); error('Wrong size on dofNumber'); end
            absDofNum = zeros(size(dofNumber, 1), 1);
            for i = 1:length(absDofNum)
                absDofNum(i) = obj.absDof(dofNumber(i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% OTHER METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = offsetNodes(obj, location)
            %    obj = offSetNodes(location)
            %
            % Offset all of the nodes in the model by the given location.
            
            % Check validity of input
            if (sum(size(location(:)) ~= [3, 1]) ~= 0)
                sz = size(location);
                error('"Location" has size [%i, %i], requires size [3, 1] or [1, 3]', ...
                    sz(1), sz(2));
            end
            
            % Modify nodes
            obj.nodes(:, 2:4) = bsxfun(@plus, obj.nodes(:, 2:4), location(:)');
            
        end
        function highlightNodes(obj, nodes, varargin)
            %    highlightNodes(obj, nodes, varargin)
            %
            % Highlight a subset of node numbers in a (hopefully existing)
            % patch plot. Use the string 'highlightSpec' to specify the
            % plot type.
            
            plotNodes = obj.getAbsNodeNumber(nodes);
            nodeLocs = obj.nodes(plotNodes, 2:4);
            set(gca, 'nextplot', 'add');
            tempHandle = plot3(nodeLocs(:, 1), nodeLocs(:, 2), nodeLocs(:, 3));
            
            for i = 1:length(varargin)/2
                set(tempHandle, varargin{2*i - 1}, varargin{2*i});
            end
            
        end
            
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% ANALYSIS STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        function [x, t] = dynamicFree(obj, boundary, x0, T, dt)
            %    [x, t] = dynamicFree(obj, boundary, x0, T, dt)
            %  
            % Perform linear/nonlinear static analysis of an Abaqus part.
            %
            % INPUT
            %
            % boundary: Vector of NODES to be fully fixed.
            %
            % x0:       [nDof x 1] Vector of initial conditions for model.
            %           Must include all DOF in structure.
            %
            % T:        [scalar] Time period of integration.
            %
            % dt:       [scalar] Starting/maximum time increment.
            %
            %
            % OUTPUT
            %
            % x:        [nDof x nT] Part displacements vs time.
            %
            % t:        [1 x nT] Time increments.
            
            obj.matDir = cd;
            % Input handling
            nBC = size(boundary, 1);
            
            [nDof, ~] = size(x0);
            if(nDof ~= obj.N) % Check on number of DOF in model
                error('Incorrect DOF number supplied');
            end
            
            % Pre-allocate variables
            stepname = [obj.modelname, '_dynamicFree']; % Name of this step
            icNodes = floor(obj.dof); icDof = round((obj.dof - icNodes)*10);

            
            % Write the base input file.
            obj.writeInp();
            
            % Open a step input file for writing in the working directory
            fid = fopen(fullfile(obj.workingDir, [stepname, '.inp']), 'w');
            
            % Write step
            try % Wrap in a try block so we can close file on error
                fprintf(fid, '*Include, INPUT=%s.inp\n', obj.modelname);
                fprintf(fid, '**NL DYNAMIC FREE STEP FROM ABINT.M\n');
                
                % First write a static step with the IC's to deflect model
                % to initial state
                fprintf(fid, '*STEP, NLGEOM=YES, NAME=INITIALCONDS\n');
                fprintf(fid, '*STATIC\n');
                fprintf(fid, '**INITIAL INCREMENT, TIME PERIOD, MINIMUM INCREMENT, MAXIMUM INCREMENT\n');
                fprintf(fid, '0.5, 1, 0.1, 1\n');
                
                % Write initial conditions
                fprintf(fid, '**INITIAL CONDITIONS\n');
                fprintf(fid, '*Boundary, TYPE=DISPLACEMENT, OP=NEW\n');
                for j = 1:obj.N
                    if sum(icNodes(j) == boundary); x0(j) = 0; end
                    fprintf(fid, '%i, %i, , %.6E\n', icNodes(j), icDof(j), x0(j));
                end
                fprintf(fid, '*END STEP\n');

                % Now write a free static step to get the response
                fprintf(fid, '*STEP, NLGEOM=YES, INC=%i, NAME=DYNAMICFREE\n', ...
                    round(1.1*T/dt));

                fprintf(fid, '*DYNAMIC, ALPHA=0, DIRECT\n');
%                 fprintf(fid, '**INITIAL INCREMENT, TIME PERIOD, MINIMUM INCREMENT, MAXIMUM INCREMENT\n');
%                 fprintf(fid, '%.6E, %.6E, %.6E, %.6E\n', dt, T, dt/1000, dt);
                fprintf(fid, '%.6E, %.6E\n', dt, T);

                fprintf(fid, '*OUTPUT, FIELD, FREQUENCY=1\n');
                fprintf(fid, '*NODE OUTPUT\n');
                fprintf(fid, 'U\n');
                fprintf(fid, '*ELEMENT OUTPUT\n');
                fprintf(fid, 'S\n');

                % Re-write initial conditions
                fprintf(fid, '**BOUNDARY CONDITIONS\n');
                fprintf(fid, '*Boundary, TYPE=DISPLACEMENT, OP=NEW\n');
                for j = 1:nBC
                    fprintf(fid, '%i, 1, 6, 0\n', boundary(j));
                end
                
                
                fprintf(fid, '*END STEP\n');

            catch err
                % Close file and exit on error
                fclose(fid);
                warning('Write error on static step; closing file and rethrowing');
                rethrow(err);
            end
            
            fclose(fid); % Close step file
            fprintf('Write complete to %s\\%s.inp\n', obj.workingDir, stepname);
            
            % Run job
            execStr = sprintf('abaqus j="%s" cpus=4 inter double ask_delete=off', ...
                stepname);
            obj.runJob(execStr, obj.workingDir);
            
            % Extract results
            cd(obj.workingDir);
            
            execStr = sprintf('abaqus python "%s" "%s" fieldOutput DYNAMICFREE all', ...
                fullfile(obj.pythonpath, 'odb2mat_v2.py'), stepname);
            status = system(execStr);
            if status ~= 0; cd(obj.matDir), error('Python read error'); end
            
            loadName = sprintf('%sFromAbaqus', stepname);
            dat = load(loadName);
            outputs = dat.fieldOutputs.DYNAMICFREE;
            
            % Clean file directory and leave this place
            delete([loadName '.mat']);
            % Perform file cleanup
            extensions = {'.com', '.dat', '.msg', '.prt', '.sim', ...
                '.sta'};
            for i = 1:length(extensions);
                delete([stepname, extensions{i}]);
            end
            cd(obj.matDir);
            
            % Assign results
            % The field names are not sorted at all, which is awesome.
            framenames = fieldnames(outputs);
            nframes = length(framenames);
            
            x = zeros(obj.N, nframes);
            t = zeros(1, nframes);
            
            for i = 1:nframes;
                framename = framenames{i};
                field = outputs.(framename); 
                U = field.U;
                UR = field.UR;
                
                frameEnd = framename(7:end);
                
                if strcmpi(frameEnd, 'last')
                    fn = nframes;
                else
                    fn = str2double(frameEnd) + 1;
                end
                t(fn) = field.TIME;

                % Place into full DOF vector
                for j = 1:obj.n % For each node
                    x(6*(j - 1) + (1:6), fn) = [U.data(j, :)'; UR.data(j, :)'];
                end
            end
            
        end
        
        
        function [disp, cf] = static(obj, boundary, loads, dof, method, springs)
            %    [disp, cf] = static(obj, boundary, loads, dof, method, springs)
            %  
            % Perform linear/nonlinear static analysis of an Abaqus part.
            %
            % INPUT
            %
            % boundary: [nBC x 7] Matrix of nodes (first column) and
            %           associated boundary conditions (remaining columns.)
            %           Non-zero entries in columns 2-7 correspond to a
            %           fixed degree of freedom at that node. Pass empty
            %           variable to ignore boundary conditions.
            %
            % load:     [nDof x nL] Matrix of load cases for simulation,
            %           arranged as column vectors. Does NOT have to 
            %           include all DOF in the model.
            %
            % dof:      [nDof x 1] Vector of absolute DOF numbers
            %           corresponding to the load DOF.
            %
            % method:   [string] 'linear' or 'nonlinear' analysis.
            %           Default nonlinear.
            %
            % springs: 
            % OUTPUT
            %
            % disp:     [nDof n nL] Part displacements, returned in the
            %           same DOF order as the supplied loading.
            
            obj.matDir = cd;
            % Input handling
            nBC = size(boundary, 1);
            
            [nDof, nL] = size(loads);
            if(nDof > obj.N) % Check on number of DOF in model
                error('Too many DOF supplied for load case');
            end
            
            dof = dof(:);
            if (length(dof) ~= nDof) % Check on loads vs. DOF size
                error('Size of loads matrix does not match size of DOF vector');
            end
            
            if (nargin < 5); method = 'nonlinear'; end
            if (nargin < 6); springs = []; end
            
            % Pre-allocate variables
            stepname = [obj.modelname, '_static']; % Name of this step
            dofAb = obj.getDofNumber(dof); bcAb = obj.getDofNumber(boundary);
            nodes = floor(dofAb); nodeDof = round((dofAb - nodes)*10);
            bcNodes = floor(bcAb); bcDof = round((bcAb - bcNodes)*10);
            
            % Write the base input file.
            obj.writeInp(springs);
            
            % Open a step input file for writing in the working directory
            fid = fopen(fullfile(obj.workingDir, [stepname, '.inp']), 'w');
            
            % Write
            method1 = true;
            try % Wrap in a try block so we can close file on error
                fprintf(fid, '*Include, INPUT=%s.inp\n', obj.modelname);
                fprintf(fid, '**NL STATIC STEP FROM ABINT.M\n**NUMBER OF LOAD CASES = %i\n', nL);
                fprintf(fid, '**LOADS APPLIED TO %i DOF\n', nDof);
                % Write boundary conditions
                fprintf(fid, '**BOUNDARY CONDITIONS\n');
                fprintf(fid, '*Boundary, TYPE=DISPLACEMENT, OP=NEW\n');
                for j = 1:nBC
                    fprintf(fid, '%i, %i, , 0\n', bcNodes(j), bcDof(j));
                end
                
                
                if (method1) % Old single-step method
                    for i = 1:nL
                        
                        
                        
                        if (strcmpi(method, 'nonlinear'));
                            fprintf(fid, '*STEP, NLGEOM=YES, NAME=STATIC-%05i\n', i);
                        elseif strcmpi(method, 'linear');
                            fprintf(fid, '*STEP, PERTURBATION, NAME=STATIC-%05i\n', i);
                        else
                            error('Invalid analysis method specified');
                        end
                        
                        fprintf(fid, '*STATIC\n');
                        fprintf(fid, '**INITIAL INCREMENT, TIME PERIOD, MINIMUM INCREMENT, MAXIMUM INCREMENT\n');
                        fprintf(fid, '0.5, 1, 1E-6, 1\n');
                        fprintf(fid, '*OUTPUT, FIELD, FREQUENCY=99999\n');
                        fprintf(fid, '*NODE OUTPUT\n');
                        fprintf(fid, 'U, CF\n');
%                         fprintf(fid, '*ELEMENT OUTPUT\n');
%                         fprintf(fid, 'S\n');
                        
                        % Write loads
                        fprintf(fid, '*CLOAD, OP=new\n');
                        for j = 1:nDof
                            fprintf(fid, '%i, %i, %.6E\n', nodes(j), nodeDof(j), loads(j, i));
                        end
                        

                        fprintf(fid, '*END STEP\n');
                        
  
                    end
                 end
                
            catch err
                % Close file and exit on error
                fclose(fid);
                warning('Write error on static step; closing file and rethrowing');
                rethrow(err);
            end
            
            fclose(fid); % Close step file
            fprintf('Write complete to %s\\%s.inp\n', obj.workingDir, stepname);
            
            % Run job
            execStr = sprintf('abaqus j="%s" cpus=4 inter double ask_delete=off', ...
                stepname);
            obj.runJob(execStr, obj.workingDir);
            
            % Extract results
            cd(obj.workingDir);
            execStr = sprintf('abaqus python "%s" "%s" fieldOutput', ...
                fullfile(obj.pythonpath, 'odb2mat_v2.py'), stepname);
            status = system(execStr);
            if status ~= 0; cd(obj.matDir), error('Python read error'); end
            
            loadName = sprintf('%sFromAbaqus', stepname);
            dat = load(loadName);
            
            % Clean file directory and leave this place
            delete([loadName '.mat']);
            % Perform file cleanup
            extensions = {'.com', '.dat', '.msg', '.prt', '.sim', ...
                '.sta'};
            for i = 1:length(extensions);
                delete([stepname, extensions{i}]);
            end
            cd(obj.matDir);
            
            % Assign results
            stepnames = sort(fieldnames(dat.fieldOutputs));
            disp = zeros(size(loads));
            cf = zeros(size(loads));
            
            for i = 1:length(stepnames);
                stepname = stepnames{i};
                field = getfield(dat.fieldOutputs, stepname);
                U = field.frame_last.U;
                UR = field.frame_last.UR;
                CF = field.frame_last.CF;
                CM = field.frame_last.CM;
                % Place into full DOF vector
                
                for j = 1:obj.n % For each node
                    disp(6*(j - 1) + (1:6), i) = [U.data(j, :)'; UR.data(j, :)'];
                    cf(6*(j - 1) + (1:6), i) = [CF.data(j, :)'; CM.data(j, :)'];
                end
            end
            
            % Put into correct order
            disp = disp(dof, :);
            cf = cf(dof, :);
        end
        
        
        function [M, K] = matrix(obj, boundary, springs)
            %    [M, K] = matrix(obj, boundary, springs)
            %  
            % Get the M and K matrices along with the DOF vector of an
            % Abaqus part.
            %
            % INPUT
            %
            % boundary: [nBC x 7] Matrix of nodes (first column) and
            %           associated boundary conditions (remaining columns.)
            %           Non-zero entries in columns 2-7 correspond to a
            %           fixed degree of freedom at that node.
            %
            % OUTPUT
            %
            % M:   [nDOF x nDOF] Sparse mass matrix
            %
            % K:   [nDOF x nDOF] Sparse stiffness matrix
            %
            % DOF: [nDOF x 1] DOF matrix
            
            obj.matDir = cd;
            % Input handling
            if isempty(boundary)
                nBC = 0;
            elseif (size(boundary, 2) ~= 7); 
                error('Need 7 columns in boundary'); 
            else
                nBC = size(boundary, 1);
            end
            
            if (nBC > obj.n); error('More boundary nodes than model nodes'); end
            stepname = [obj.modelname, '_matrix']; % Name of this step
            
            if (nargin < 3)
                springs = [];
            end
            
            % Write the base input file.
            obj.writeInp(springs);
            
            % Open a step input file for writing in the working directory
            fid = fopen(fullfile(obj.workingDir, [stepname, '.inp']), 'w');
            
            % Write
            try % Wrap in a try block so we can close file on error
                fprintf(fid, '*Include, INPUT=%s.inp\n', obj.modelname);
                fprintf(fid, '*Step\n');
                fprintf(fid, '**MATRIX EXTRACTION STEP FROM ABINT.M\n');
                fprintf(fid, 'SUBTITLE=%s\n', obj.modelname);
                fprintf(fid, '*MATRIX GENERATE, MASS, STIF\n');
                fprintf(fid, '*MATRIX OUTPUT, MASS, STIF\n');
                fprintf(fid, '**BOUNDARY CONDITIONS\n');
                fprintf(fid, '*Boundary, TYPE=DISPLACEMENT\n');
                
                for i = 1:nBC
                    bcRow = find(boundary(i, 2:end));
                    if any(diff(bcRow) > 1) % If we have separate indices
                        % Write out separately
                        for j = 1:length(bcRow)
                            fprintf(fid, '%i, %i, ,\n', boundary(i, 1), ...
                                bcRow(j));
                        end
                    else
                        % Write out all at once
                        fprintf(fid, '%i, %i, %i, \n', boundary(i, 1), ...
                            bcRow(1), bcRow(end));
                    end
                end
                
                fprintf(fid, '*END STEP\n');
                
            catch err
                % Close file and exit on error
                fclose(fid);
                warning('Write error on matrix step; closing file and rethrowing');
                rethrow(err);
            end
            
            fclose(fid); % Close step file
            fprintf('Write complete to %s\\%s.inp\n', obj.workingDir, stepname);
            
            % Run job
            execStr = sprintf('abaqus j="%s" cpus=4 inter ask_delete=off double', ...
                stepname);
            obj.runJob(execStr, obj.workingDir);
            
            % Extract results
            cd(obj.workingDir);
            
            % file=extension(file,'mtx');
            fid = fopen([stepname, '_MASS1.mtx'], 'r');
            xMass = fscanf(fid, '%g,%g,%g,%g,%g', [5 inf])';
            fclose(fid);
            fid = fopen([stepname, '_STIF1.mtx'], 'r');
            xStiff = fscanf(fid, '%g,%g,%g,%g,%g', [5 inf])';
            fclose(fid);
            
            % Clean file directory and leave this place
            delete([stepname, '_MASS1.mtx']);
            delete([stepname, '_STIF1.mtx']);
            extensions = {'.com', '.dat', '.msg', '.odb', '.prt', '.sim', ...
                '.sta'};
            for i = 1:length(extensions);
                delete([stepname, extensions{i}]);
            end
            cd(obj.matDir);
            
            
            % Get matrices
            M = obj.sparseMake(xMass);
            K = obj.sparseMake(xStiff);
            
        end
        
        function spMat = sparseMake(obj, x)
            massVec = x(:,5);
            dof_key = x(:,1:4);
            dofPerNode = max(max(x(:,[2 4])));
            
            % Make a dummy vector of sorted nodes
            sortedNodes = sort(unique(dof_key(:, 1)));
            tic
            row_id = (obj.ifind(sortedNodes, dof_key(:, 1)) - 1)*dofPerNode + dof_key(:, 2);
            clm_id = (obj.ifind(sortedNodes, dof_key(:, 3)) - 1)*dofPerNode + dof_key(:, 4);
            spMat = sparse(row_id, clm_id, massVec, dofPerNode*obj.n, dofPerNode*obj.n);
            spMat = spMat + spMat' - diag(diag(spMat));
        end
        
        
        function [phi, frequencies, dof] = modal(obj, boundary, nModes, springs)
            %    [phi, frequencies, dof] = modal(obj, boundary, nModes, springs)
            %  
            % Get the modeshapes of an Abaqus part. 
            %
            % INPUT
            %
            % boundary: [nBC x 7] Matrix of nodes (first column) and
            %           associated boundary conditions (remaining columns.)
            %           Non-zero entries in columns 2-7 correspond to a
            %           fixed degree of freedom at that node. Pass empty
            %           variable to ignore boundary conditions.
            %
            % nModes:   Integer, number of modes to recover OR 2 x 1 vector
            %           of a frequency range to use for recovery.
            %
            % OUTPUT
            %
            % phi:         [nDOF x nModes] Modal matrix
            %
            % frequencies: [nModes x 1] Natural frequency (Hz) of each
            %              mode.
            % 
            % dof:         [nDof x 1] DOF numbering scheme for the part
            
            obj.matDir = cd;
            % Input handling
            if isempty(boundary)
                nBC = 0;
            elseif (size(boundary, 2) ~= 7); 
                error('Need 7 columns in boundary'); 
            else
                nBC = size(boundary, 1);
            end
            if(length(nModes) == 1) % Number of modes
                fLower = []; fUpper = [];
            elseif (length(nModes) == 2) % Frequency range
                if (nModes(1) < 1E-6); fLower = [];
                    if (nModes(1) ~= 0)warning('Adjusting fLower to 0 Hz'); end
                else
                    fLower = nModes(1); 
                end
                fUpper = nModes(2); nModes = [];
            else
                error('"nModes" input too long, need either 1 or 2 values');
            end
            
            if (nBC > obj.n); error('More boundary nodes than model nodes'); end
            stepname = [obj.modelname, '_modal']; % Name of this step
            
            if (nargin < 4)
                springs = [];
            end
            
            % Write the base input file.
            obj.writeInp(springs);
            
            % Open a step input file for writing in the working directory
            fid = fopen(fullfile(obj.workingDir, [stepname, '.inp']), 'w');
            
            % Write
            try % Wrap in a try block so we can close file on error
                fprintf(fid, '*Include, INPUT=%s.inp\n', obj.modelname);
                fprintf(fid, '*Step\n');
                fprintf(fid, '**MODAL STEP FROM ABINT.M\n**NMODES = %i\n', nModes);
                fprintf(fid, 'SUBTITLE=%s\n', obj.modelname);
                fprintf(fid, '*Frequency, EIGENSOLVER=LANCZOS, NORMALIZATION=MASS\n');
                fprintf(fid, '**num_EV, F_lower, F_upper, shift_point, block_size\n');
                fprintf(fid, '%i, %.6E, %.6E, , ,\n', nModes, fLower, fUpper);
                fprintf(fid, '**BOUNDARY CONDITIONS\n');
                fprintf(fid, '*Boundary, TYPE=DISPLACEMENT\n');
                for i = 1:nBC
                    bcRow = find(boundary(i, 2:end));
                    if any(diff(bcRow) > 1) % If we have separate indices
                        % Write out separately
                        for j = 1:length(bcRow)
                            fprintf(fid, '%i, %i, ,\n', boundary(i, 1), ...
                                bcRow(j));
                        end
                    else
                        % Write out all at once
                        fprintf(fid, '%i, %i, %i, \n', boundary(i, 1), ...
                            bcRow(1), bcRow(end));
                    end
                end
                
                fprintf(fid, '*END STEP\n');
                
            catch err
                % Close file and exit on error
                fclose(fid);
                warning('Write error on modal step; closing file and rethrowing');
                rethrow(err);
            end
            
            fclose(fid); % Close step file
            fprintf('Write complete to %s\\%s.inp\n', obj.workingDir, stepname);
            
            % Run job
            execStr = sprintf('abaqus j="%s" cpus=4 inter ask_delete=off', ...
                stepname);
            obj.runJob(execStr, obj.workingDir);
            
            % Extract results
            cd(obj.workingDir);
            execStr = sprintf('abaqus python "%s" "%s" shapes', ...
                fullfile(obj.pythonpath, 'odb2mat_v2.py'), stepname);
            status = system(execStr);
            if status ~= 0; cd(obj.matDir), error('Python read error'); end
            
            loadName = sprintf('%sFromAbaqus', stepname);
            load(loadName);
            
            % Clean file directory and leave this place
            delete([loadName '.mat']);
            % Perform file cleanup
            extensions = {'.com', '.dat', '.msg', '.prt', '.sim', ...
                '.sta'};
            for i = 1:length(extensions);
                delete([stepname, extensions{i}]);
            end
            cd(obj.matDir);
            
            % The object is called "mode" in the workspace; has fields
            % "shape" and "freqs"
            phi = mode.shapes;
            frequencies = mode.freqs;
            dof = mode.dof;
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% READ/WRITE INPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        function obj = getDof(obj)
            %    getDof(obj)
            %
            % Obtain the node.dof numbering scheme of the model through a
            % call to the modal analysis step.
            
            % Get node.dof vector
            [~, ~, obj.dof] = obj.modal([], 3);
            obj.N = length(obj.dof);
            
            % Construct absolute DOF number
            obj.absDof = containers.Map(obj.dof, 1:obj.N);
        end
        
        function writeInp(obj, springs)
            %    writeInp(obj, springs)
            %
            % Write a base model input file using a stored Abaqus model.
            % Note that no loads or boundary conditions are written in this
            % step; those must be included in analysis steps.
            
            % Spring argument is optional
            if (nargin < 2); springs = []; end
            
            % Check if a model actually exists
            assert(~isempty(obj.modelname), 'Model is empty');
            
            % Open a file in the working directory
            fid = fopen(fullfile(obj.workingDir, [obj.modelname, '.inp']), 'w');
            
            try % Wrap in a try block so we can close file on error
                % Start writing
                fprintf(fid, '** AUTO-GENERATED INPUT FILE FROM ABINT.M **\n');
                fprintf(fid, '** ABAQUS VERSION 6.12 **\n');
                fprintf(fid, '*Heading\n** Job name: %s Model name: %s\n', ...
                    obj.modelname, obj.modelname);
                
                % Write nodes
                fprintf(fid, '*Node\n');
                for i = 1:obj.n
                    fprintf(fid, '%i, %.6E, %.6E, %.6E\n', obj.nodes(i, 1), ...
                        obj.nodes(i, 2), obj.nodes(i, 3), obj.nodes(i, 4));
                end
                
                % Write elements
                eleTypes = unique(obj.elements.eleTypes);
                nTypes = length(eleTypes);
                
                for i = 1:nTypes
                    type = eleTypes{i};
                    fprintf(fid, '*Element, type=%s\n', type);
                    for j = 1:obj.M
                        if strcmpi(obj.elements.eleTypes{j}, type)
                            connect = obj.elements.eles(j, 3:end);
                            connect = connect(connect ~= 0);
                            fprintf(fid, '%i', obj.elements.eles(j, 1));
                            for k = 1:length(connect)
                                fprintf(fid, ', %i', connect(k));
                            end
                            fprintf(fid, '\n');
                        end
                    end
                end
                
                % Write element sets
                setNames = fieldnames(obj.eleSets);
                for i = 1:length(setNames)
                    name = setNames{i};
                    set = obj.eleSets.(name);
                    fprintf(fid, '*Elset, ELSET=%s\n', name);
                    for j = 1:length(set);
                        fprintf(fid, '%i\n', set(j));
                    end
                end
                
                % Write node sets
                setNames = fieldnames(obj.nodeSets);
                for i = 1:length(setNames)
                    name = setNames{i};
                    set = obj.nodeSets.(name);
                    fprintf(fid, '*Nset, NSET=%s\n', name);
                    for j = 1:length(set);
                        fprintf(fid, '%i\n', set(j));
                    end
                end
                
                % Write sections
                keys = obj.sections.keys;
                for i = 1:length(keys)
                    key = keys{i};
                    text = obj.sections(key);
                    fprintf(fid, '%s\n', key);
                    for j = 1:length(text);
                        fprintf(fid, '%s\n', text{j});
                    end
                end
                
                % Write surfaces
                keys = obj.surfaces.keys;
                for i = 1:length(keys)
                    key = keys{i};
                    text = obj.surfaces(key);
                    fprintf(fid, '%s\n', key);
                    for j = 1:length(text);
                        fprintf(fid, '%s\n', text{j});
                    end
                end
                
                % Write ties
                keys = obj.ties.keys;
                for i = 1:length(keys)
                    key = keys{i};
                    text = obj.ties(key);
                    fprintf(fid, '%s\n', key);
                    for j = 1:length(text);
                        fprintf(fid, '%s\n', text{j});
                    end
                end
                    
                % Write material
                matNames = fieldnames(obj.materials);
                for i = 1:length(matNames);
                    matName = matNames{i};
                    fprintf(fid, '*Material, NAME=%s\n', matName);
                    mat = obj.materials.(matName);
                    props = fieldnames(mat);
                    for j = 1:length(props)
                        prop = props{j};
                        values = mat.(prop);
                        fprintf(fid, '*%s\n', prop);
                        for k = 1:length(values);
                            fprintf(fid, '%.6E, ', values(k));
                        end
                        fprintf(fid, '\n');
                    end
                end
                
                % OPTIONAL - Spring section
                if ~isempty(springs)
                    % First write 1 spring per DOF. We have the desired DOF
                    % and stiffness associated with each spring in the
                    % "springs" structure array
                    
                    nSprings = length(springs.k);
                    lastElement = max(obj.elements.eles(:, 1));
                    for i = 1:nSprings
                        % Get attached node and DOF
                        fprintf(fid, '*ELEMENT, type=SPRING1, ELSET=spring-%i\n', i);
                        springDof = obj.getDofNumber(springs.dof(i));
                        node = floor(springDof);
                        springDof = round(10*(springDof - node));
                        
                        % Write node
                        fprintf(fid, '%i, %i\n', lastElement + i, node);
                        
                        % Assign spring property
                        fprintf(fid, '*Spring, ELSET=spring-%i\n', i);
                        fprintf(fid, '%i\n', springDof);
                        fprintf(fid, '%.6E\n', springs.k(i));
                    end
                    
                end
                
            catch err
                fclose(fid);
                warning('Write error; closing file and rethrowing');
                rethrow(err);
            end
            fclose(fid);
            fprintf('Write complete to %s\\%s.inp\n', obj.workingDir, obj.modelname);
        end
        
        function obj = readInp(obj, modelpath, modelname)
            obj.matDir = cd; addpath(obj.matDir);
            if strcmpi(modelname(end-3:end), '.inp');
                modelname = modelname(1:end-4);
            end
            obj.modelname = modelname;
            
            % Run a datacheck of the given input file to generate a base
            % .odb
            execStr = sprintf('abaqus j="%s" inter ask_delete=off', obj.modelname);
            status = obj.runJob(execStr, modelpath); % Run the desired job
            if status ~= 0; return; end
            
            % Now use Python to read the .ODB. Call the "odb2mat_v2.py". 
            cd(modelpath);
            execStr = sprintf('abaqus python "%s" "%s" model', ...
                fullfile(obj.pythonpath, 'odb2mat_v2.py'), obj.modelname);
            status = system(execStr);
            if status ~= 0; cd(obj.matDir), error('Python read error'); end
            
            loadName = sprintf('%sFromAbaqus', obj.modelname);
            load(loadName);
            delete([loadName '.mat']);
            
            % Assign ODB properties to ABINT
            obj.nodes = FE.nodes;    
            obj.elements = FE.elements;     
            obj.nodeSets = FE.nodeSets; 
            obj.eleSets = FE.eleSets;  
            obj.materials = FE.materials;
            
            obj.n = size(obj.nodes, 1);  
            obj.M = size(obj.elements.eles, 1);  
            
            % Build the absolute DOF container.map object
            obj.absNodes = containers.Map(obj.nodes(:, 1), 1:obj.n);
            
            % Fix the element types
            eleTypes = cell(obj.M, 1);
            for i = 1:obj.M
                eleTypes{i} = obj.elements.eleTypes(i, :);
            end
            obj.elements.eleTypes = eleTypes;
            
            % Open .inp file and read section data
            fid = fopen([obj.modelname, '.inp'], 'r');
            line = fgetl(fid);
            
            while ischar(line)
                if length(line) < 4; line = fgetl(fid); continue; end
                lineCheck = strsplit(line, ',');
                
                % *TIE checks
                if length(lineCheck{1}) < 4; line = fgetl(fid); continue; end
                if strcmpi(lineCheck{1}(1:4), '*TIE');
                    tieKey = line; tieCell = cell(0); flag = true;
                    while flag;
                        line = fgetl(fid);
                        if strcmpi(line(1), '*')
                            flag = false;
                        else
                            tieCell{end + 1} = line;
                        end
                    end
                    obj.ties(tieKey) = tieCell;
                    continue
                end
                
                
                % *SURFACE checks
                if length(lineCheck{1}) < 8; line = fgetl(fid); continue; end
                if strcmpi(lineCheck{1}(1:8), '*SURFACE');
                    surfKey = line; surfCell = cell(0); flag = true;
                    while flag;
                        line = fgetl(fid);
                        if strcmpi(line(1), '*')
                            flag = false;
                        else
                            surfCell{end + 1} = line;
                        end
                    end
                    obj.surfaces(surfKey) = surfCell;
                    continue
                end
                
                % Section checks
                if length(lineCheck{1}) < 7; line = fgetl(fid); continue; end
                if strcmpi(lineCheck{1}(end-7:end), ' section');
                    sectionKey = line; sectionCell = cell(0); flag = true;
                    while flag;
                        line = fgetl(fid);
                        if strcmpi(line(1), '*')
                            flag = false;
                        else
                            sectionCell{end + 1} = line;
                        end
                    end
                    obj.sections(sectionKey) = sectionCell;
                    continue
                end
                
                % Get next line and continue otherwise
                line = fgetl(fid);
            end
            
            fclose(fid);
        
            % Perform file cleanup
            extensions = {'.com', '.dat', '.msg', '.odb', '.prt', '.sim', ...
                '.sta'};
            for i = 1:length(extensions);
                delete([obj.modelname, extensions{i}]);
            end
            cd(obj.matDir);
        end
        
        
        function status = runJob(obj, execStr, workingDir)
            
            cd(workingDir);
            status = system(execStr);
            
            % Error handling
            if (status ~= 0)
                warning('Analysis error; scanning %s.dat for information...', obj.modelname);
                fprintf('Echoing ***ERROR lines:\n');
                fid = fopen(sprintf('%s.dat', obj.modelname));
                line = fgetl(fid);
                k = 1;
                while ischar(line)
                    lineCheck = strrep(line, ' ', '');
                    if length(lineCheck) < 8; line = fgetl(fid); k = k + 1; continue; end
                    if strcmpi(lineCheck(1:8), '***ERROR');
                        fprintf('Line %i: %s\n', k, line);
                    end
                    line = fgetl(fid);
                    k = k + 1;
                end
                fclose(fid);
                cd(obj.matDir);
            end
            cd(obj.matDir);
        end
            
        
        
        function deform(obj, X, maxNorm)
            %    deform(obj, X, maxNorm)
            %
            % Deform each DOF of the object by the values given in X. X
            % must be a column vector of size N x 1 where N is the number
            % of DOF in the model. If a patch handle "h" does not already
            % exist for the object, it will be created in the current
            % figure window. "maxNorm" is used for coloring purposes, to
            % normalize all deflections relative to the largest in the
            % model.
            
            % Check size of x
            if ((size(X, 1) ~= obj.N) || (size(X, 2) ~= 1))
                error('Incorrect size of deform input X');
            end
            
            % If no patch handle exists, call plot
            if isempty(obj.h); obj = obj.plot(); end
            
            obj.h.Vertices = [obj.nodes(:, 2) + X(1:6:end), ...
                              obj.nodes(:, 3) + X(2:6:end), ...
                              obj.nodes(:, 4) + X(3:6:end)]; % This only works for 6-DOF nodes
                          
            % Set colors
            X = X(:);
            mag = dot([X(1:6:end), X(2:6:end), X(3:6:end)], [X(1:6:end), X(2:6:end), X(3:6:end)], 2);
            mag = mag/maxNorm;
            set(obj.h, 'FaceVertexCData', mag, 'FaceColor', 'Interp');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% OPERATOR OVERLOADS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = isequal(obj, rhs)
            % Check for same name
            if ~strcmpi(obj.modelname, rhs.modelname); out = 0; return; end
            % Check for same nodes, up to an offset
            offset1 = obj.nodes(1, 2:4); offset2 = rhs.nodes(1, 2:4);
            nodesLhs = bsxfun(@minus, obj.nodes(:, 2:4), offset1);
            nodesRhs = bsxfun(@minus, rhs.nodes(:, 2:4), offset2);
            if (norm(nodesLhs - nodesRhs) ~= 0); out = 0; return; end
            
            % Check same elements
            if (norm(obj.elements.eles - rhs.elements.eles) > 0)
                out = 0;
                return
            end
            
            % Check for same materials
            if ~(isequal(obj.materials, rhs.materials))
                out = 0;
                return
            end
            out = 1;
        end
        
        function save(obj)
            abint = obj;
            save(fullfile(obj.filepath, obj.modelname), 'abint');
        end
        
        function disp(obj)
            % Display the object in a readable format
            fprintf('\nABINT Object\n');
            if isempty(obj.modelname)
                fprintf('Model is empty\n');
            else
                fprintf('\t%s\n', obj.modelname);
                fprintf('\t%i nodes\n', obj.n);
                if ~isempty(obj.N); fprintf('\t%i DOF\n', obj.N); end
                fprintf('\t%i elements\n', obj.M);
                fprintf('\t%i section assignments\n', length(obj.sections.keys));
            end
            
            byteSize = whos('obj'); byteSize = byteSize.bytes;
            fprintf('\tSize of %.3f mb\n\n', byteSize/1024^2);
        end
        
        
        
        function obj = plot(obj, varargin)
            % Overloaded call to plot function. Additional name/value plot
            % property pairs will be carried forward.
            
            % Use "patch" to make the plot. Key tasks are stripping out all
            % empty values from the element connectivity matrix, and
            % re-ordering the values in the element connectivity matrix to
            % correspond to a sorted order of nodes.
            
            % Make a dummy vector of sorted nodes
            [sortedNodes, sortedNodesIdx] = sort(unique(obj.nodes(:, 1)));
            
            eles = obj.ifind(sortedNodes, obj.elements.eles(:, 3:end));
            
            nanSum = sum(~isnan(eles), 1); % Which columns have only NaN?
            eles = eles(:, nanSum ~= 0); % Keep only those columns
            nodeLocs = obj.nodes(sortedNodesIdx, 2:4);
            
            % Get a length scale for plotting purposes
            obj.h = patch('Vertices', nodeLocs,'Faces', eles);
            
            % Formatting
            set(obj.h, 'MarkerEdgeColor',[0 0 0], ...
                'Marker','None', ...
                'FaceColor',[0.6 0.6 1], ...
                'FaceAlpha', 0.5, ...
                'EdgeColor',[0 0 0.3]);
            view([-45, 45])
            axis equal
            xlabel('\bfX'); ylabel('\bfY'); zlabel('\bfZ');

        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% STATIC METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        
        function idx = ifind(array, values)
            % Return the indices in "array" corresponding to occurence of
            % "values." Detects only the first occurence of each value
            idx = zeros(size(values));
            for i = 1:numel(values)
                thisIndex = find((abs(array - values(i)) < eps), 1, 'first');
                if isempty(thisIndex); thisIndex = NaN; end
                idx(i) = thisIndex;
            end
        end
        
        function obj = loadobj(obj)
            % Have to check on our filepaths in case we're loading on a
            % different computer
            obj.matDir = cd;
            try 
                cd(obj.workingDir);
            catch err % If we can't change to the absolute directory, try changing to a relative directory
                warning('Absolute working directory specification %s does not exist; checking relative directory %s', ...
                    obj.workingDir, obj.workingDirRel);
                try 
                    cd(obj.workingDirRel);
                    obj.workingDir=cd;
                catch err;
                    warning('Relative working directory specification %s does not exist', obj.workingDirRel);
                    newdir = uigetdir(cd, 'Select ABAQUS Working Directory');
                    if (newdir == 0); error('No working directory selected; ABINT will not function correctly'); end
                    obj.workingDir = newdir;
                end 
            end
            cd(obj.matDir);
            
            try 
                cd(obj.filepath);
            catch err
                warning('Filepath designation %s does not exist, checking relative directory %s', ...
                    obj.filepath, obj.filepathRel);
                try 
                    cd(obj.filepathRel)
                    obj.filepath=cd;
                catch err
                    warning('Relative filepath directory specification %s does not exist', obj.filepathRel);
                    newdir = uigetdir(cd, 'Select ABINT Save Directory');
                    if (newdir == 0); error('No save directory selected; ABINT will not function correctly'); end
                    obj.filepath = newdir;
                end
            end
            cd(obj.matDir);
            
            try 
                cd(obj.pythonpath);
            catch err
                warning('Filepath designation python path %s does not exist; checking relative directory %s', ...
                    obj.pythonpath, obj.pythonpathRel);
                try 
                    cd(obj.pythonpathRel);
                    obj.pythonpath=cd;
                catch err
                    warning('Relative python path directory specification %s does not exist', obj.pythonpathRel);
                    newdir = uigetdir(cd, 'Select Python Directory');
                    if (newdir == 0); error('No python path selected; ABINT will not function correctly'); end
                    obj.pythonpath = newdir;
                end
            end
            cd(obj.matDir);
        end
    end
                
        
    
end