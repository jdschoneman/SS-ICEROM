% CCSS(filepath, filename) 
% 
% Craig-Bampton with characteristic constraint mode substructuring 
% routine for CBICE models. 
%
% 
 

classdef CCSS < CBSS

    properties(GetAccess = 'public', SetAccess = 'public')
        
        Mcc;
        Kcc;
        phiCB;      % Original Craig-Bampton modal matrix
        psiCC;      % Matrix of characteristic constraint modes
        fnCC;       % "Frequencies" of the CC modes
        phiCC;      % Resulting modal matrix using characteristic constraint modes
        Tcc;        % CC transformation matrix
        modalCB;    % Modal DOF indices for CB matrix
        boundaryCB; % Boundary DOF indices for CB matrix
        nCharDOF;   % Number of characteristic constraint modes retained
        ccBasis;
        
        % Logical values
        cbModel = false; % Last SS was a CB model
        ccModel = false; % Last SS was a CC model
        
       
    end
    
    
    methods
        % Constructor
        function obj = CCSS(filepath, filename)
            
            % Call CBSS constructor
            obj@CBSS(filepath, filename);
            
        end
        
        function obj = initNLSystem(obj)
            %
            % Initialize an NL system by setting all part nonlinearities to
            % "0."
            for i = 1:obj.npart
                obj.parts(i).part = obj.parts(i).part.zeroNLROM(obj.nCharDOF);
            end
        end
        
        function obj = copyNLPart(obj, basePart, copyPart)
            %    obj = copyNLPart(obj, basePart, copyPart)
            %
            % Copy the NL model of "basePart" to that of "copyPart."
            obj.parts(copyPart).part.NLROM = obj.parts(basePart).part.NLROM;
        end
        
        function [obj, fitData] = addNLPartIS(obj, ab, nodeset, bcset, refNode, basePart, fiMind, ccMind, thick, quads, triples, flist, fitData)
            %       [obj, fitData] = addNLPartIS(obj, ab, nodeset, bcset, refNode, basePart, fiMind, ccMind, thick, quads, triples, flist, fitData)
            %
            % Build nonlinear model of a part using the "in situ" FEA
            % modeling procedure.
            
            % Input handling
            if (nargin < 10); quads = obj.parts(basePart).part.options.quads; end
            if (nargin < 11); triples = obj.parts(basePart).part.options.triple; end
            if (nargin < 12); flist = []; end
            if (nargin < 13); fitData = []; end
                
            % Get the "part-ordered" set of CC modes
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            nC = obj.nCharDOF;
            if (length(ccMind) > nC); error('More CC modes requested than available'); end
            TccOrdered = obj.L*obj.Tcc(cc2cbConvert, (end - nC + 1):end);
            
            % Get global DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(basePart - 1)
                endDOF = currDOF + obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = basePart;
            nM = obj.parts(i).nModalDOF; nB = obj.parts(i).nBoundaryDOF;
            if (length(fiMind) > nM); error('More FI modes requested than available'); end
            endDOF = currDOF + nM + nB;
            
            % Create a basis in terms of the part's FI and CC modes
            psiChar = TccOrdered(((currDOF + nM):endDOF - 1), :);
%             fiBasis = [eye(nM); zeros(nB, nM)];
            ccBasis = [zeros(nM, nC); psiChar];
%             basis = [fiBasis(:, fiMind), ccBasis(:, fiMind)];
            rPart = obj.parts(i).part.abint.getAbsNodeNumber(obj.parts(i).refNode);
            rAb = ab.getAbsNodeNumber(refNode);
            refLoc = [rPart, rAb];
            % "refLoc" is [component reference node, ABINT reference node]
            [obj.parts(i).part, fitData] = obj.parts(i).part.ficcFitIS(refLoc, ...
                ab, nodeset, bcset, fiMind, ccMind, ccBasis, thick, quads, triples, flist, fitData);            
                    
        end
        
         function [obj, fitData] = addNLPartQRIS(obj, ab, nodeset, bcset, refNode, basePart, qrInds, thick, quads, triples, flist, fitData)
            %       [obj, fitData] = addNLPartQRIS(obj, ab, nodeset, bcset, refNode, basePart, qrInds, thick, quads, triples, flist, fitData)
            %
            % Build nonlinear model of a part using the "in situ" FEA
            % modeling procedure.
            
            % Input handling
            if (nargin < 9); quads = obj.parts(basePart).part.options.quads; end
            if (nargin < 10); triples = obj.parts(basePart).part.options.triple; end
            if (nargin < 11); flist = []; end
            if (nargin < 12); fitData = []; end
                
            % Get the "part-ordered" set of CC modes
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            nC = obj.nCharDOF;
            TccOrdered = obj.L*obj.Tcc(cc2cbConvert, (end - nC + 1):end);
            
            % Get global DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(basePart - 1)
                endDOF = currDOF + obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = basePart;
            nM = obj.parts(i).nModalDOF; nB = obj.parts(i).nBoundaryDOF;
            
            endDOF = currDOF + nM + nB;
            
            % Create a basis in terms of the part's FI and CC modes
            psiChar = TccOrdered(((currDOF + nM):endDOF - 1), :);
%             fiBasis = [eye(nM); zeros(nB, nM)];
            ccBasis = [zeros(nM, nC); psiChar];
%             basis = [fiBasis(:, fiMind), ccBasis(:, fiMind)];
            rPart = obj.parts(i).part.abint.getAbsNodeNumber(obj.parts(i).refNode);
            rAb = ab.getAbsNodeNumber(refNode);
            refLoc = [rPart, rAb];
            % "refLoc" is [component reference node, ABINT reference node]
            [obj.parts(i).part, fitData] = obj.parts(i).part.qrFitIS(refLoc, ...
                ab, nodeset, bcset, qrInds, ccBasis, thick, quads, triples, flist, fitData);            
                    
        end
        
        function [obj, fitData] = addNLPart(obj, basePart, fiMind, ccMind, springs, thick, quads, triples, flist, fitData)
            %    [obj, fitData] = addNLPart(obj, basePart, fiMind, ccMind, springs, thick, quads, triples, flist, fitData)
            %
            % Build nonlinear model of a part
            
            % Input handling
            if (nargin < 7); quads = obj.parts(basePart).part.options.quads; end
            if (nargin < 8); triples = obj.parts(basePart).part.options.triple; end
            if (nargin < 9); flist = []; end
            if (nargin < 10); fitData = []; end
                
            % Get the "part-ordered" set of CC modes
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            nC = obj.nCharDOF;
            if (length(ccMind) > nC); error('More CC modes requested than available'); end
            TccOrdered = obj.L*obj.Tcc(cc2cbConvert, (end - nC + 1):end);
            
            % Get global DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(basePart - 1)
                endDOF = currDOF + obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = basePart;
            nM = obj.parts(i).nModalDOF; nB = obj.parts(i).nBoundaryDOF;
            if (length(fiMind) > nM); error('More FI modes requested than available'); end
            endDOF = currDOF + nM + nB;
            
            % Get the static equivalent stiffness of the part if springs
            % are not supplied
%             if (nargin < 3) % None of this works right now
%                 xCBDOF = (currDOF + nM):(endDOF - 1);
%                 xADOF = [];
%                 for j = xCBDOF
%                     xADOF = [xADOF; find(obj.L(j, :))];
%                 end
%                 xEDOF = (1:(size(obj.L, 2)))';
%                 for j = xADOF
%                     xEDOF(j) = 0;
%                 end
%                 xEDOF(xEDOF == 0) = [];
% 
%                 % Partition matrices
%                 KA = obj.L'*obj.Ksys*obj.L;
%                 cDof = obj.parts(i).part.linModel.connDofAbs;
%                 Kaa = KA(xADOF, xADOF) - obj.parts(i).part.linModel.K(cDof, cDof); % Kae = KA(xADOF, xEDOF);
%     %             Kea = Kae'; Kee = KA(xEDOF, xEDOF);
%     %             Te = -Kee\Kea;
%     %             Ka = [Kaa, Kae; Kea, Kee]*[eye(length(xADOF), size(Te, 2)); Te];
%     %             
%                 springs.dof = obj.parts(i).part.linModel.connDofAbs;
%                 nNodes = length(springs.dof)/6;
%                 springs.k = full(diag(Kaa));
%             end
            
            % Create a basis in terms of the part's FI and CC modes
            psiChar = TccOrdered(((currDOF + nM):endDOF - 1), :);
%             fiBasis = [eye(nM); zeros(nB, nM)];
            ccBasis = [zeros(nM, nC); psiChar];
%             basis = [fiBasis(:, fiMind), ccBasis(:, fiMind)];
            [obj.parts(i).part, fitData] = obj.parts(i).part.ficcFit(fiMind, ccMind, ccBasis, springs, thick, quads, triples, flist, fitData);            
            
            
        end
        
          function [obj, fitData] = addNLPartSVD(obj, basePart, nSvd, springs, thick, quads, triples, flist, fitData)
            %      [obj, fitData] = addNLPartSVD(obj, basePart, nSvd, springs, thick, quads, triples, flist, fitData)
            %
            % Build nonlinear model of a part
            
            % Input handling
            if (nargin < 6); quads = []; end
            if (nargin < 7); triples = []; end
            if (nargin < 8); flist = []; end
            if (nargin < 9); fitData = []; end
                
            % Get the "part-ordered" set of CC modes
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            nC = obj.nCharDOF;
            
            TccOrdered = obj.L*obj.Tcc(cc2cbConvert, (end - nC + 1):end);
            
            % Get global DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(basePart - 1)
                endDOF = currDOF + obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = basePart;
            nM = obj.parts(i).nModalDOF; nB = obj.parts(i).nBoundaryDOF;

            endDOF = currDOF + nM + nB;
            
            
            
            % Get the static equivalent stiffness of the part if springs
            % are not supplied
            if (nargin < 3) % None of this works right now
                xCBDOF = (currDOF + nM):(endDOF - 1);
                xADOF = [];
                for j = xCBDOF
                    xADOF = [xADOF; find(obj.L(j, :))];
                end
                xEDOF = (1:(size(obj.L, 2)))';
                for j = xADOF
                    xEDOF(j) = 0;
                end
                xEDOF(xEDOF == 0) = [];

                % Partition matrices
                KA = obj.L'*obj.Ksys*obj.L;
                cDof = obj.parts(i).part.linModel.connDofAbs;
                Kaa = KA(xADOF, xADOF) - obj.parts(i).part.linModel.K(cDof, cDof); % Kae = KA(xADOF, xEDOF);
    %             Kea = Kae'; Kee = KA(xEDOF, xEDOF);
    %             Te = -Kee\Kea;
    %             Ka = [Kaa, Kae; Kea, Kee]*[eye(length(xADOF), size(Te, 2)); Te];
    %             
                springs.dof = obj.parts(i).part.linModel.connDofAbs;
                nNodes = length(springs.dof)/6;
                springs.k = full(diag(Kaa));
            end
            
            % Create a basis in terms of the part's FI and CC modes
            psiChar = TccOrdered(((currDOF + nM):endDOF - 1), :);
%             fiBasis = [eye(nM); zeros(nB, nM)];
            ccBasis = [zeros(nM, nC); psiChar];
%             basis = [fiBasis(:, fiMind), ccBasis(:, fiMind)];
            [obj.parts(i).part, fitData] = obj.parts(i).part.svdFit(nSvd, ccBasis, springs, thick, quads, triples, flist, fitData);            
            
            
          end
        
          
          function [obj, fitData] = addNLPartQR(obj, basePart, qrInds, springs, thick, quads, triples, flist, fitData)
            %      [obj, fitData] = addNLPartQR(obj, basePart, qrInds, springs, thick, quads, triples, flist, fitData)
            %
            % Build nonlinear model of a part
            
            % Input handling
            if (nargin < 6); quads = []; end
            if (nargin < 7); triples = []; end
            if (nargin < 8); flist = []; end
            if (nargin < 9); fitData = []; end
                
            % Get the "part-ordered" set of CC modes
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            nC = obj.nCharDOF;
            
            TccOrdered = obj.L*obj.Tcc(cc2cbConvert, (end - nC + 1):end);
            
            % Get global DOF of the part in order to obtain the correct CC
            % mode partition
            currDOF = 1;
            for i = 1:(basePart - 1)
                endDOF = currDOF + obj.parts(i).nModalDOF + obj.parts(i).nBoundaryDOF;
                currDOF = endDOF;
            end
            
            i = basePart;
            nM = obj.parts(i).nModalDOF; nB = obj.parts(i).nBoundaryDOF;

            endDOF = currDOF + nM + nB;
            
            
            
            % Get the static equivalent stiffness of the part if springs
            % are not supplied
            if (nargin < 3) % None of this works right now
                xCBDOF = (currDOF + nM):(endDOF - 1);
                xADOF = [];
                for j = xCBDOF
                    xADOF = [xADOF; find(obj.L(j, :))];
                end
                xEDOF = (1:(size(obj.L, 2)))';
                for j = xADOF
                    xEDOF(j) = 0;
                end
                xEDOF(xEDOF == 0) = [];

                % Partition matrices
                KA = obj.L'*obj.Ksys*obj.L;
                cDof = obj.parts(i).part.linModel.connDofAbs;
                Kaa = KA(xADOF, xADOF) - obj.parts(i).part.linModel.K(cDof, cDof); % Kae = KA(xADOF, xEDOF);
    %             Kea = Kae'; Kee = KA(xEDOF, xEDOF);
    %             Te = -Kee\Kea;
    %             Ka = [Kaa, Kae; Kea, Kee]*[eye(length(xADOF), size(Te, 2)); Te];
    %             
                springs.dof = obj.parts(i).part.linModel.connDofAbs;
                nNodes = length(springs.dof)/6;
                springs.k = full(diag(Kaa));
            end
            
            % Create a basis in terms of the part's FI and CC modes
            psiChar = TccOrdered(((currDOF + nM):endDOF - 1), :);
%             fiBasis = [eye(nM); zeros(nB, nM)];
            ccBasis = [zeros(nM, nC); psiChar];
%             basis = [fiBasis(:, fiMind), ccBasis(:, fiMind)];
            [obj.parts(i).part, fitData] = obj.parts(i).part.qrFit(qrInds, ccBasis, springs, thick, quads, triples, flist, fitData);            
            
            
        end
        
        
        
        % Build CC model
        function obj = buildCCSystem(obj, nModes, nCModes)
            %    obj = buildCCSystem(nModes, nCModes)
            %
            % Build a system using characteristic constraint modes. Inputs
            % are the number of modes to compute and the number of CC modes
            % to retain.
            
            % Build B, L, matrix etc.
            obj = obj.buildBL();
            
            % Partition mass and stiffness matrices to the boundary
            % set.
            obj.nCharDOF = nCModes;
            dofCount = sum(obj.L, 1);
            obj.modalCB = find(dofCount == 1);   % Modal DOF in the model
            obj.boundaryCB = find(dofCount > 1); % Connection DOF in the model
            [~, cc2cbConvert] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            
            % Get assembled mass and stiffness matrices
            Ma = obj.L'*obj.Msys*obj.L;
            Ka = obj.L'*obj.Ksys*obj.L;
            
            % Partition them to boundary DOF
            Mii = Ma(obj.modalCB, obj.modalCB);
            Kii = Ka(obj.modalCB, obj.modalCB);
            Mib = Ma(obj.modalCB, obj.boundaryCB);
            Kib = Ka(obj.modalCB, obj.boundaryCB);
            Mbb = Ma(obj.boundaryCB, obj.boundaryCB);
            Kbb = Ka(obj.boundaryCB, obj.boundaryCB);
            
            % Perform eigenvalue analysis
            [obj.psiCC, obj.fnCC] = eigs(Kbb, Mbb, nCModes, 'SM');
            obj.fnCC = sqrt(diag(obj.fnCC))/2/pi;
            
            % Construct CC transformation matrix
            I = eye(obj.nModalDOF); 
            Zi = zeros(obj.nModalDOF, nCModes);
            Zb = zeros(length(obj.boundaryCB), obj.nModalDOF);
            obj.Tcc = [I, Zi; Zb, obj.psiCC];
            
            % Further transform CB matrix; perform modal
            obj.Mcc = obj.Tcc'*[Mii, Mib; Mib', Mbb]*obj.Tcc;
            obj.Kcc = obj.Tcc'*[Kii, Kib; Kib', Kbb]*obj.Tcc;
            [obj.phiCC, fn] = eigs(obj.Kcc, obj.Mcc, nModes, 'SM');
            obj.fn = sqrt(diag(fn))/2/pi;
            
            % Sort modes to ascending order
            [obj.fn, modeOrder] = sort(obj.fn);
            obj.phiCC = obj.phiCC(:, modeOrder);
            
            
            % Convert the CC modal matrix back to Craig-Bampton "phi" for
            % plotting and MAC comparison
            obj.phi = obj.Tcc*obj.phiCC; % Transform
            obj.phi = obj.L*obj.phi(cc2cbConvert, :); % Re-order
            
            % Assign modal amplitudes to each part
            obj = obj.assignAmplitudes();
            
            obj.cbModel = false;
            obj.ccModel = true;
        end
        
        % Slightly altered "comparePlot" routine
        function fig = comparePlot(obj, maxErr, nModes, cmap)
            if nargin < 2; maxErr = 2.5; end
            if nargin < 3; nModes = size(obj.phiFull, 2); end
            if nargin < 4; cmap = colormap('jet'); end
            
            fig = comparePlot@CBSS(obj, maxErr, nModes, cmap);
            
            % Modify title
            if (obj.ccModel)
                children = get(fig, 'children');
                ax = children(3);
                titlestr = get(ax, 'title'); titlestr = titlestr.String;
                titlestr = sprintf('%s Reduced to %i CC Modes', titlestr, ...
                    obj.nCharDOF);
                title(ax, titlestr);
            end
        end
        
        % Disallow calling "buildCBSystem" from this class, so there's no
        function obj = buildCBSystem(obj, nModes)
%             error('CB Substructuring Disabled for CCSS Class');
            obj.cbModel = true;
            obj.ccModel = false;
            obj = buildCBSystem@CBSS(obj, nModes);
        end
        
        function plotConstraintmode(obj, mode, factor, varargin)
            %    plotConstraintmode(obj, mode, factor, varargin)
            % Plot the deformation of a given mode for all parts
            if (nargin < 3); factor = 1; end
            
            cMode = obj.Tcc(:, obj.nModalDOF + mode);
            [~, cc2cb] = sort([obj.modalCB(:); obj.boundaryCB(:)]);
            cMode = cMode(cc2cb);
            cMode = obj.L*cMode;
            currDOF = 1;
            for i = 1:obj.npart
                % Call the CBICE "deform" function with the appropriate set
                % of modal amplitudes
                endDOF = currDOF + obj.parts(i).nBoundaryDOF + obj.parts(i).nModalDOF;
                amps = cMode(currDOF:endDOF - 1);
                obj.parts(i).part.plot(amps*factor);
                currDOF = endDOF;
            end
            title(sprintf('\\bfConstraint Mode %i', mode));
        end
        
        
        function save(obj, filename)
            obj.filepath = 'CCSSs';
            if (nargin == 2)
                obj.filename = filename;
            end
            ccss = obj; %#ok
            save(fullfile(obj.filepath, obj.filename), 'ccss');
        end
        
        function disp(obj)
            % Display the object in a readable format
            fprintf('\nCCSS Object\n');
            if obj.npart == 0;
                fprintf('No included parts\n');
            else
                fprintf('\t%s\n', obj.filename);
                fprintf('\t%i parts\n', obj.npart);
                fprintf('\t%i Modal DOF; %i Boundary DOF\n', ...
                    obj.nModalDOF, obj.nBoundaryDOF);
                fprintf('\t%i Characteristic Constraint Modes\n', ...
                    obj.nCharDOF);
            end
            
            byteSize = whos('obj'); byteSize = byteSize.bytes;
            fprintf('\tSize of %.3f mb\n\n', byteSize/1024^2);
        end
            
    end
end