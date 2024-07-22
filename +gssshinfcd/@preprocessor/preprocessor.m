classdef preprocessor 
% PREPROCESSOR Processor formulating a problem into an SDP

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
        % prob - The gssshinfcd problem
        % The H-infinity controller synthesis problem.
        prob
        
        % opts - The options for controller synthesis
        opts
    end
    
    methods
        %% Constructor
        function obj = preprocessor(problem, options)
        % PREPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'gssshinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            
            assert(isstruct(options), 'options should be a structure.');
            obj.opts = options;
            
        end
        
        %% Solution procedure
        function K = solve(obj)
            % SOLVE Executes the solution routines
                
                % solve with the standard method
               K = gssshinfcd.standard.nlcov.nlcov(obj.prob, obj.opts); 
             
        end
        
    end
    
end