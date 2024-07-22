classdef postprocessor 
% POSTPROCESSOR Postprocessing steps to the controller matrices originating from the standard problem

% This file is part of gssshinfcd.
% Copyright (c) 2024, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% gssshinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% gssshinfcd is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
% License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with gssshinfcd. If not, see <https://www.gnu.org/licenses/>.
    
    properties
        % prob - The hinfcd problem
        % The extended H-infinity controller synthesis problem.
        prob
        
        % K - The controller resulting from the standard problem
        K
        
        % perf - Achieved performance
        perf
    end
    
    methods
        %% Constructor
        function obj = postprocessor(problem, K)
        % POSTPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'gssshinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            obj.K = K; 
            assert(~all(structfun(@isempty,K)),'Solver did not return a solution. The problem is infeasible, ill-posed or ill-conditioned.'); 
            
        end
        
        %% Compensation for the feedthrough component of the generalized plant
        function obj = compensateDyu(obj)
        % COMPENSATEDYU Compensates for direct feedthrough in the generalized plant, if present
            
            % helper function:  check for constant splines
            function [b,val] = isconstant(spline)
                if isa(spline,'splines.Function')
                    c = spline.coeff.data;
                    if ndims(c)<4; c = reshape(c,[1 size(c)]); end
                    for k=1:size(c,1)
                        for l=2:size(c,2)
                            dc(k,l-1,:,:) = c(k,l,:,:)-c(k,l-1,:,:);
                        end
                    end
                    b = ~(any(dc(:) > sqrt(eps)));
                    if b
                        cog = cellfun(@(x) mean(x.domain.data), spline.tensor_basis.bases);
                        val = spline.eval(cog);
                    else
                        val = spline;
                    end
                else
                    b = true;
                    val = spline;
                end
            end
            
            % first check if Dyu is constant
            [constant,Dyu] = isconstant(obj.prob.P.D((obj.prob.nz()+1):end,(obj.prob.nw()+1):end));
            
            % if not, then make the algebraic loop if Dyu ~= 0
            if ~constant && any(Dyu(:))
                for kk=1:length(obj.K)
                    obj.K(kk).K = obj.K(kk).K([1:obj.prob.nu() 1:end],[1:obj.prob.ny() 1:end]); 
                    obj.K(kk).K.D = obj.K(kk).K.D * blkkdiag(eye(obj.prob.ny()),-Dyu,eye(2*obj.prob.n())); 
                    obj.K(kk).fb = blkkdiag(eye(obj.prob.nu()),obj.K(kk).fb);
                end
            end
            
        end
           
    end
    
end