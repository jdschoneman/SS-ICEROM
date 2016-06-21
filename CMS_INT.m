

classdef CMS_INT
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        m; % Number of parts
        nCC; % Number of CC DOF
        ccBasis; % Reduced basis for the global CC modes
        fiDofs;  % Vector of FI DOF counts for each part
        nFIDof;  % Count of total FI DOF
        NLROMs; % Nonlinear models of each part (structure array)
        Mhat;   % Linear model mass matrix
        Khat;   % Linear model stiffness matrix
        L;      % Substructure assembly matrix
        psiCC;
        modalCB;  % Indices of FI modal DOF
        boundaryCB; % Indices of constraint DOF
        q0;
        N1hat;
        N2hat;
    end
    
    methods
        
        function obj = CMS_INT(ccss)
            
            obj.m = length(ccss.parts);
            obj.Mhat = ccss.Mcc;
            obj.Khat = ccss.Kcc;
            
            obj.nCC = ccss.nCharDOF;
            obj.nFIDof = size(obj.Mhat, 1) - obj.nCC;
            
            obj.psiCC = ccss.psiCC;
            obj.modalCB = ccss.modalCB;
            obj.boundaryCB = ccss.boundaryCB;
            obj.ccBasis = eye(obj.nCC);
            
            obj.q0 = zeros(obj.nCC + obj.nFIDof, 1);
            
            obj.fiDofs = zeros(obj.m, 1);
            obj.NLROMs = struct('N1', [], 'N2', [], 'Nlist', [], 'theta', [], ...
                                'thetaInv', []);
            
            % Build up the NLROMs and L matrix
            obj.L = zeros(obj.nFIDof + obj.nCC*obj.m, ...
                          obj.nFIDof + obj.nCC);
            nEq = 0; nDof = 0;
            for i = 1:obj.m
                obj.fiDofs(i) = ccss.parts(i).nModalDOF;    
                for j = 1:obj.fiDofs(i)
                    obj.L(nEq + j, nDof + j) = 1;
                end
                nEq = nEq + j; nDof = nDof + j;
                for j = 1:obj.nCC
                    obj.L(nEq + j, obj.nFIDof + j) = 1;
                end
                nEq = nEq + j;
                    
                hold.N1 = ccss.parts(i).part.NLROM.N1;
                hold.N2 = ccss.parts(i).part.NLROM.N2;
                hold.Nlist = ccss.parts(i).part.NLROM.Nlist;
                
                if ~isfield(ccss.parts(i).part.NLROM, 'theta');
                    hold.theta = eye(size(hold.N1(:, :, 1)));
                else
                    hold.theta = ccss.parts(i).part.NLROM.theta;
%                     [U, S, V] = svd(ccss.parts(i).part.NLROM.basis, 'econ');
%                     hold.theta = V';
                end
                hold.thetaInv = inv(hold.theta);
                
                
                obj.NLROMs(i) = hold;
                
            end
            
            [obj.N1hat, obj.N2hat] = update(obj, [obj.q0; obj.q0]);
        end
            
            
        function fnl = fint_nl(obj, z)
            
            N = size(z, 1)/2;
            if (checkState(obj, z))
                [obj.N1hat, obj.N2hat] = update(obj, z);
                obj.q0 = z(1:N);
            end
            
            fnl = [1/2*obj.N1hat + 1/3*obj.N2hat]*z(1:N);
            
        end
        
        function J = dfint_nl(obj, z)
            N = size(z, 1)/2;
            if (checkState(obj, z))
                [obj.N1hat, obj.N2hat] = update(obj, z);
                obj.q0 = z(1:N);
            end
            
            J = [obj.N1hat + obj.N2hat];
            J = [J, zeros(size(J))];
            
        end
        
        function E = Energy_nl(obj, z)
            N = size(z, 1)/2;
            if (checkState(obj, z))
                [obj.N1hat, obj.N2hat] = update(obj, z);
                obj.q0 = z(1:N);
            end
            
            E = z(1:N)'*(1/6*obj.N1hat + 1/12*obj.N2hat)*z(1:N);
        end
        
        function [N1, N2] = update(obj, z)
            N = size(z, 1)/2;
            q = z(1:N);
%             try
            qu = obj.L*q; % Unconstrained modal coordinates

%             etaCC = pinv(obj.ccBasis)*qCC;
            
            N1 = zeros(obj.m*obj.nCC + obj.nFIDof);
            N2 = N1;
            
            % Build up a vector of forces in the fixed interface modes
            
            for j = 1:obj.m
                
                if (j == 1)
                    indFI = 1:obj.fiDofs(j);
                else
                    indFI = (indj(end) + 1):(indj(end) + obj.fiDofs(j));
                end
                indCC = (indFI(end) + 1):(indFI(end) + obj.nCC);
                indj = [indFI, indCC];
                nDof = size(obj.NLROMs(j).theta, 1);
                N1j = zeros(nDof);
                N2j = zeros(nDof);
                
                Nlistj = obj.NLROMs(j).Nlist;
                
                qj = qu([indFI, indCC]); % Unconstrained modal coordinates
                % Convert to the SVD coordinates
                eta = obj.NLROMs(j).theta*qj;
                
                eta2 = eta(Nlistj(:,1)).*eta(Nlistj(:,2));

                for k = 1:nDof;
                    N1j(:, k)= obj.NLROMs(j).N1(:, :, k)*eta;
                    N2j(:, k)= obj.NLROMs(j).N2(:, :, k)*eta2;
                end
                

                N1j = obj.NLROMs(j).theta'*N1j*obj.NLROMs(j).theta;
                N2j = obj.NLROMs(j).theta'*N2j*obj.NLROMs(j).theta;
                
                
                N1(indj, indj) = N1j;
                N2(indj, indj) = N2j;
%                 keyboard
            end
            
            N1 = obj.L'*N1*obj.L;
            N2 = obj.L'*N2*obj.L;
            
        end
        
        function new = checkState(obj, z)
            N = size(z, 1)/2;
            qCheck = z(1:N);
            new = true;
            return;
            if (norm(qCheck - obj.q0) > eps)
                new = true;
            else
                new = false;
            end
            
        end
        
    end
end