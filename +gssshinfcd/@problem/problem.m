classdef problem 
% PROBLEM H-infinity controller synthesis problem

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
        % P - Generalized plant (ny measured outputs, nu control inputs)
        % A state-space realization of the augmented (weighted) generalized
        % plant, in the form of a gsss model. 
        P
        ny
        nu

        % specs - Desired performance specifications
        % Structure array indicating the performance channels and their relative weights.
        specs
        
        % region - Desired closed-loop D-stability specificions
        % Structure array indicating the closed-loop pole region constraints.
        region

        % watchdog - A watchdog for monitoring
        % The watchdog keeps track of all processes that happen on objects of this class. 
        watchdog
        
    end
    
    methods
        %% Constructor, input checking, standard form
        
        function obj = problem(watchdog, P, ny, nu, specs, region)
        % PROBLEM Constructor
        
            % construct the object
            assert(isa(watchdog,'gssshinfcd.watchdog.watchdog'), 'Invalid watchdog.');
            obj.watchdog = watchdog; 
            obj.P = P; 
            obj.ny = ny;
            obj.nu = nu; 
            obj.specs = specs;
            obj.region = region; 
            
            % check the validity of the input
            obj.checkInput(); 
            
            % preprocess the object
            obj = obj.sortchannels(); 
            
        end

        function checkInput(obj)
        % CHECKINPUT Checks the validity of the input arguments
        
            % check member types
            obj.watchdog.assertType(obj.watchdog, 'gssshinfcd.watchdog.watchdog', 'Invalid watchdog type.');
            obj.watchdog.assertType(obj.P, 'gssshinfcd.gsss', 'P is not gsss (gain-scheduled state-space) model.'); 
            obj.watchdog.assertCond(isnumeric(obj.ny) && isscalar(obj.ny) && obj.ny > 0 && mod(obj.ny,1)==0, 'ny must be a strictly positive integer.'); 
            obj.watchdog.assertCond(isnumeric(obj.nu) && isscalar(obj.nu) && obj.nu > 0, 'nu must be a strictly positive integer.');
            obj.watchdog.assertCond(obj.nu < size(obj.P,2), 'The plant does not have any exogenous (performance) inputs.');
            obj.watchdog.assertCond(obj.ny < size(obj.P,1), 'The plant does not have any performance outputs.');
            obj.watchdog.assertType(obj.specs, 'struct', 'specs is not a structure.'); 
            obj.watchdog.assertType(obj.region, 'struct', 'region is not a structure.'); 
            
            % check specifications
            obj.watchdog.assertField(obj.specs, 'in', 'specs requires a field ''in''.');
            obj.watchdog.assertField(obj.specs, 'out', 'specs requires a field ''out''.'); 
            obj.watchdog.assertField(obj.specs, 'weight', 'specs requires a field ''weight''.'); 
            cellfun(@(x) obj.watchdog.assertNonEmpty(x, 'A specification cannot have zero performance inputs.'), {obj.specs.in}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.in can only contain numeric values.'), {obj.specs.in}); 
            cellfun(@(x) obj.watchdog.assertCond(all(mod(x,1)==0) && all(x>0), 'specs.in can only contain integer values.'), {obj.specs.in});
            cellfun(@(x) obj.watchdog.assertNonEmpty(x, 'A specification cannot have zero performance outputs.'), {obj.specs.out}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.out can only contain numeric values.'), {obj.specs.out});
            cellfun(@(x) obj.watchdog.assertCond(all(mod(x,1)==0) && all(x>0), 'specs.out can only contain integer values.'), {obj.specs.out});
            cellfun(@(x) obj.watchdog.assertNonEmpty(x,'A specification requires a weight. Set to 0 to handle it as a normalized constraint.'), {obj.specs.weight}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.weight can only contain numeric values.'), {obj.specs.weight}); 
            cellfun(@(x) obj.watchdog.assertCond(all(x>=0) && isscalar(x), 'specs.weight can only positive scalars or 0.'), {obj.specs.weight});
            obj.watchdog.assertField(obj.region, 'M', 'region requires a field ''M''.'); 
            obj.watchdog.assertField(obj.region, 'L', 'region requires a field ''L''.'); 
            cellfun(@(x) obj.watchdog.assertCond(isnumeric(x) || isa(x,'splines.Function'), 'region.M can only contain numeric values or spline functions.'), {obj.region.M});
            cellfun(@(x) obj.watchdog.assertCond(isnumeric(x) || isa(x,'splines.Function'), 'region.L can only contain numeric values or spline functions.'), {obj.region.L});
            cellfun(@(x,y) obj.watchdog.assertCond(size(x,1)==size(y,1), 'region.M and region.L must have the same number of rows.'), {obj.region.M}, {obj.region.L});
            
            obj.watchdog.assertCond(max(vertcat(obj.specs.in))<=size(obj.P,2)-obj.nu, 'A performance input exceeding the number of performance inputs was requested.'); 
            obj.watchdog.assertCond(max(vertcat(obj.specs.out))<=size(obj.P,1)-obj.ny, 'A performance output exceeding the number of performance outputs was requested.'); 
            obj.watchdog.warnCond(all(ismember(1:(size(obj.P,2)-obj.nu),vertcat(obj.specs.in))),'Some performance inputs do not appear in the specifications. They will therefore be ignored.');
            obj.watchdog.warnCond(all(ismember(1:(size(obj.P,1)-obj.ny),vertcat(obj.specs.out))),'Some performance outputs do not appear in the specifications. They will therefore be ignored.');
        end
        
        function obj = sortchannels(obj)
        % SORTCHANNELS Sorts the specifications and its performance inputs and outputs
            obj.P = obj.P([vertcat(obj.specs.out); reshape((end-obj.ny+1):end,obj.ny,1)],[vertcat(obj.specs.in); reshape((end-obj.nu+1):end,obj.nu,1)]);
            chi = 0; cho = 0;
            for k=1:length(obj.specs)
                obj.specs(k).in = (chi+1):(chi+length(obj.specs(k).in)); chi = chi+length(obj.specs(k).in);
                obj.specs(k).out = (cho+1):(cho+length(obj.specs(k).out)); cho = cho+length(obj.specs(k).out); 
            end
        end
        
        %% Sampling time
        
        function Ts = Ts(obj)
        % TS Returns the sampling time of the plant 
            Ts = obj.P.Ts;
        end
        
        %% Dimensions
        
        function nw = nw(obj)
        % NW Returns number of exogenous (performance) inputs
            nw = size(obj.P,2)-obj.nu;
        end
        
        function nz = nz(obj)
        % NZ Returns number of performance outputs
            nz = size(obj.P,1)-obj.ny;
        end
      
        function n = n(obj)
        % N Returns order of the plant
            n = size(obj.P.A,1);
        end
        
        %% Matrices of P        
        
        function A = A(obj)
        % A Returns the A-matrix of the generalized plant
            A = obj.P.A;
        end
        
        function Bw = Bw(obj)
        % Bw Returns the Bw-matrix of the generalized plant
            Bw = obj.P.B(:,1:obj.nw());
        end
        
        function Bu = Bu(obj)
        % BU Returns the Bu-matrix of the generalized plant
            Bu = obj.P.B(:,(obj.nw()+1):end);
        end
        
        function Cz = Cz(obj)
        % Cz Returns the Cz-matrix of the generalized plant
            Cz = obj.P.C(1:obj.nz(),:);
        end
        
        function Cy = Cy(obj)
        % CY Returns the Cy-matrix of the generalized plant
            Cy = obj.P.C((obj.nz()+1):end,:);
        end
        
        function Dzw = Dzw(obj)
        % DEV Returns the Dev-matrix of the generalized plant
            Dzw = obj.P.D(1:obj.nz(),1:obj.nw());
        end
        
        function Dzu = Dzu(obj)
        % DEU Returns the Deu-matrix of the generalized plant
            Dzu = obj.P.D(1:obj.nz(),(obj.nw()+1):end);
        end
        
        function Dyw = Dyw(obj)
        % DYV Returns the Dyv-matrix of the generalized plant
            Dyw = obj.P.D((obj.nz()+1):end,1:obj.nw());
        end
        
        function Dyu = Dyu(obj)
        % DYU Returns the Dyu-matrix of the generalized plant
            Dyu = obj.P.D((obj.nz()+1):end,(obj.nw()+1):end);
        end
        
    end
    
end