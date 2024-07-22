function obj = setoptions(obj,s)
% SETOPTIONS Sets the options for the control problem and its handling
  
% This file is part of gssshinfcd.
% Copyright (c) 2024, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% sshinfcd is free software: you can redistribute it and/or modify it under 
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

    % get default options
    defops = options(obj);
    
    % recursively merge options
    if nargin>1
        % make sure the problem has been defined
        obj.watchdog.assertCond(~isempty(obj.problem), 'First set the problem, then the options.'); 
        
        % assign options
        obj.opts = gssshinfcd.util.mergestruct(s,defops);
        
        % check gssshinfcd specific options
        obj.watchdog.assertCond(ischar(obj.opts.synthesis.gammasolver) && ismember(obj.opts.synthesis.gammasolver,{'lmilab','cvx','yalmip'}), 'gammasolver can only take the values ''lmilab'', ''cvx'' or ''yalmip''.'); 
        obj.watchdog.assertCond(isscalar(obj.opts.synthesis.knotinsertions) && isnumeric(obj.opts.synthesis.knotinsertions) && mod(obj.opts.synthesis.knotinsertions,1)==0 && obj.opts.synthesis.knotinsertions >= 0, 'Number of knot insertions must be a positive integer.'); 
        obj.watchdog.assertCond(isscalar(obj.opts.synthesis.degreeelevations) && isnumeric(obj.opts.synthesis.degreeelevations) && mod(obj.opts.synthesis.degreeelevations,1)==0 && obj.opts.synthesis.degreeelevations >= 0, 'Number of degree elevations must be a positive integer.'); 
        obj.watchdog.assertCond(isscalar(obj.opts.synthesis.zerotol) && isnumeric(obj.opts.synthesis.zerotol) && obj.opts.synthesis.zerotol >= 0, 'Feasibility tolerance for LMIs must be positive.'); 
        obj.watchdog.assertCond(isscalar(obj.opts.synthesis.lyapunovshape) && isnumeric(obj.opts.synthesis.lyapunovshape) && ismember(obj.opts.synthesis.lyapunovshape,[1 2 3 4]),'The Lyapunov shape can only be of type 1, 2, 3 or 4.');
        
        if ~isempty(obj.opts.synthesis.lyapunovbasis) % check the tensor basis
            obj.watchdog.assertCond(isa(obj.opts.synthesis.lyapunovbasis,'splines.TensorBasis'), 'If nonempty, the B-spline tensor basis of the Lyapunov matrix must be of type splines.TensorBasis.'); 
            
            % get the tensor basis of the plant
            TBP = tensor_basis([obj.problem.P.A obj.problem.P.B ; obj.problem.P.C obj.problem.P.D]);
            
            % compare tensor basis of Lyapunov parametrization with the tensor basis of the plant 
            argP = num2cell(TBP.arguments);
            argB = num2cell(obj.opts.synthesis.lyapunovbasis.arguments);
            obj.watchdog.assertCond(~any(cellfun(@(x) contains(x,' '), argB)), 'Do not use whitespaces in parameter names. They cause troubles in the OptiSpline toolbox.'); 
            obj.watchdog.warnCond(all(ismember(argP,argB)), 'Not all scheduling parameters of the plant are in the Lyapunov matrix tensor basis you provided. As a result, the Lyapunov matrix will be constant along these dimensions.');
            iB = ismember(argB,argP); 
            b = bases(obj.opts.synthesis.lyapunovbasis);
            if ~any(iB)
                obj.watchdog.warning('There are scheduling parameters in the Lyapunov matrix tensor basis you provided that are not in the plant. These will be omitted.'); 
                obj.opts.synthesis.lyapunovbasis = [];
            elseif ~all(iB) 
                obj.watchdog.warning('There are scheduling parameters in the Lyapunov matrix tensor basis you provided that are not in the plant. These will be omitted.'); 
                obj.opts.synthesis.lyapunovbasis = splines.TensorBasis(b(iB),argB(iB)); 
                argB = num2cell(obj.opts.synthesis.lyapunovbasis.arguments);
                b = b(iB);
            end
            
            % finally, check the intervals
            bP = bases(TBP); 
            for k=1:sum(iB)
                obj.watchdog.assertCond(all(data(domain(b{k})) == data(domain(bP{strcmp(argB{k},argP)}))), 'The domain of the parameters in the Lyapunov matrix tensor basis do not match the domain of the corresponding parameters in the plant.'); 
            end
        end
        
        if ~isempty(obj.opts.synthesis.gammabasis) % check the tensor basis
            obj.watchdog.assertCond(isa(obj.opts.synthesis.gammabasis,'splines.TensorBasis'), 'If nonempty, the B-spline tensor basis of the performance variables (gamma) must be of type splines.TensorBasis.'); 
            
            % get the tensor basis of the plant
            TBP = tensor_basis([obj.problem.P.A obj.problem.P.B ; obj.problem.P.C obj.problem.P.D]);
            
            % compare tensor basis of Lyapunov parametrization with the tensor basis of the plant 
            argP = num2cell(TBP.arguments);
            argB = num2cell(obj.opts.synthesis.gammabasis.arguments);
            obj.watchdog.assertCond(~any(cellfun(@(x) contains(x,' '), argB)), 'Do not use whitespaces in parameter names. They cause troubles in the OptiSpline toolbox.'); 
            obj.watchdog.warnCond(all(ismember(argP,argB)), 'Not all scheduling parameters of the plant are in the parameterization of the performance variables (gamma). If deliberate, you can ignore this warning.');
            iB = ismember(argB,argP); 
            b = bases(obj.opts.synthesis.gammabasis);
            if ~any(iB)
                obj.watchdog.warning('There are scheduling parameters in the parameterization of the performance variables (gamma) that are not in the plant. These will be omitted.'); 
                obj.opts.synthesis.gammabasis = [];
            elseif ~all(iB) 
                obj.watchdog.warning('There are scheduling parameters in the parameterization of the performance variables (gamma) that are not in the plant. These will be omitted.'); 
                obj.opts.synthesis.gammabasis = splines.TensorBasis(b(iB),argB(iB)); 
                argB = num2cell(obj.opts.synthesis.gammabasis.arguments);
                b = b(iB);
            end
            
            % finally, check the intervals
            bP = bases(TBP); 
            for k=1:sum(iB)
                obj.watchdog.assertCond(all(data(domain(b{k})) == data(domain(bP{strcmp(argB{k},argP)}))), 'The domain of the parameters in the parameterization of the performance variables (gamma) do not match the domain of the corresponding parameters in the plant.'); 
            end
            
        end
        
        if ~isempty(obj.opts.synthesis.controllerbasis) % check the tensor basis
            obj.watchdog.assertCond(isa(obj.opts.synthesis.controllerbasis,'splines.TensorBasis'), 'If nonempty, the B-spline tensor basis of the transformed controller variables must be of type splines.TensorBasis.'); 
            
            % get the tensor basis of the plant
            TBP = tensor_basis([obj.problem.P.A obj.problem.P.B ; obj.problem.P.C obj.problem.P.D]);
            
            % compare tensor basis of Lyapunov parametrization with the tensor basis of the plant 
            argP = num2cell(TBP.arguments);
            argB = num2cell(obj.opts.synthesis.controllerbasis.arguments);
            obj.watchdog.assertCond(~any(cellfun(@(x) contains(x,' '), argB)), 'Do not use whitespaces in parameter names. They cause troubles in the OptiSpline toolbox.'); 
            obj.watchdog.warnCond(all(ismember(argP,argB)), 'Not all scheduling parameters of the plant are in the tensor basis you provided for the transformed controller variables. This will lead to conservative results.');
            iB = ismember(argB,argP); 
            b = bases(obj.opts.synthesis.controllerbasis);
            if ~any(iB)
                obj.watchdog.warning('None of the scheduling parameters in the tensor basis you provided for the transformed controller variables are in the plant. Using the default tensor basis instead.'); 
                obj.opts.synthesis.controllerbasis = [];
            elseif ~all(iB) 
                obj.watchdog.warning('There are scheduling parameters in the tensor basis you provided for the transformed controller variables that are not in the plant. These will be omitted.'); 
                obj.opts.synthesis.controllerbasis = splines.TensorBasis(b(iB),argB(iB)); 
                argB = num2cell(obj.opts.synthesis.controllerbasis.arguments);
                b = b(iB);
            end
            
            % finally, check the intervals
            bP = bases(TBP); 
            for k=1:sum(iB)
                obj.watchdog.assertCond(all(data(domain(b{k})) == data(domain(bP{strcmp(argB{k},argP)}))), 'The domain of the parameters in the tensor basis for the transformed controller variables do not match the domain of the corresponding parameters in the plant.'); 
            end
        end
    else
        obj.opts = defops;
    end
    
end