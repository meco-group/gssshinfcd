classdef gsss
% GSSS Gain-scheduled state-space model 

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
        A           % A-matrix
        B           % B-matrix
        C           % C-matrix
        D           % D-matrix
        Ts          % sampling time
        parameters  % structure containing the scheduling parameters
                    %  parameters(i).name:  name of the scheduling
                    %  parameters(i).range: [m M], m minimal value of the parameter, M its maximal value (must be finite)
                    %  parameters(i).rate: [m M], m minimal value of the parameter rate, M its maximal value (can be infinite)
    end
    
    methods
        function obj = gsss(A,B,C,D,Ts,parameters)
        % Constructor
        
            % check types and get parameters
            sys = {A,B,C,D}; p = {}; r = {}; 
            assert(any(cellfun(@(x) isa(x,'splines.Function'), sys)), 'The state-space matrices are not parameter-dependent. There are better alternatives for LTI systems than gssshinfcd - use them instead.'); 
            for i=1:4
                if isa(sys{i},'splines.Function')
                    p = [p ; num2cell(sys{i}.tensor_basis.arguments)];
                    r = [r ; cellfun(@(x) x.data, sys{i}.domain.domains', 'un', 0)];
                else
                    assert(isnumeric(sys{i}), 'State-space matrices must be numeric or a B-spline function.'); 
                end
            end
            
            % check whether discrete-time or not
            if isnumeric(Ts)
                assert(isscalar(Ts) && Ts>=0, 'Sampling time must be a real positive scalar.');
            else
                parameters = Ts; 
                Ts = 0;
            end
            
            % check all parameters and match with parameters structure
            assert(isstruct(parameters) && ~isempty(parameters), 'The parameters structure must be a nonempty structure.'); 
            assert(all(ismember({'name','range'}, fieldnames(parameters))), 'The parameters structure must have (at least) the fields ''name'' and ''range''.');
            assert(all(cellfun(@(x) ischar(x) && ~isempty(x), {parameters.name})), 'The field ''name'' of the parameter structure can only contain nonempty strings.'); 
            assert(length(unique({parameters.name})) == length({parameters.name}),'All parameters must have unique names.');
            assert(all(cellfun(@(x) isnumeric(x) && isvector(x) && length(x)==2 && x(2)>x(1), {parameters.range})), 'The field ''range'' of the parameter structure can only contain numeric vectors of length 2. In addition, they should be sorted in monotonically increasing order.'); 
            assert(~any(cellfun(@(x) contains(x,' '), {parameters.name})), 'The name of parameters cannot contain whitespaces. They cause troubles in the OptiSpline toolbox.'); 
            if isfield(parameters,'rate')
                for i=1:length(parameters)
                    if isempty(parameters(i).rate) % assume constant/frozen
                        parameters(i).rate = [0 0];
                    else
                        assert(isnumeric(parameters(i).rate) && isvector(parameters(i).rate) && length(parameters(i).rate)==2 && (~any(parameters(i).rate) || parameters(i).rate(2)>=parameters(i).rate(1)), 'The field ''rate'' of the parameter structure can only contain numeric vectors of length 2. In addition, with the exception of [0 0], they should be sorted in monotonically increasing order.');
                    end
                end
            end
            
            p = cellfun(@strip, p, 'un', 0); % remove leading and trailing whitespaces
            for i=1:length(p)
                assert(ismember(p{i},{parameters.name}), ['Parameter ''' p{i} ''' occurs in the state-space matrices but is not included in the parameter list you provided.']); 
                assert(all(r{i}(:) == parameters(find(strcmp(p{i},{parameters.name}))).range(:)), 'The field ''range'' of the parameter structure can only contain numeric vectors of length 2. In addition, they should be sorted in monotonically increasing order.'); 
            end
            
            % check whether all parameters are in the plant
            ip = ismember({parameters.name},p); 
            if ~all(ip)
                warning('Not all parameters defined in the parameter structure occur in the state-space matrices of the plant. These will be omitted.'); 
                parameters(~ip) = [];
            end
            
            % check matrix dimensions
            assert(size(A,1) == size(A,2), 'A must be square.');
            assert(size(A,1) == size(B,1), 'B must have as many rows as A.');
            assert(size(B,2) == size(D,2), 'B and D must have the same number of columns.'); 
            assert(size(C,1) == size(D,1), 'C and D must have the same number of rows.');
            assert(size(C,2) == size(A,2), 'C must have as many columns as A.'); 
            
            % assign members
            obj.A = A;
            obj.B = B; 
            obj.C = C;
            obj.D = D; 
            obj.Ts = Ts; 
            obj.parameters = parameters;
            
        end
        
        function s = size(obj,varargin)
        % SIZE Returns the size of the model
            switch nargin
                case 1
                    s = [size(obj.D,1) size(obj.D,2)];
                otherwise
                    if varargin{1}==1
                        s = size(obj.D,1);
                    elseif varargin{1}==2
                        s = size(obj.D,2);
                    else
                        error('Gain-scheduled state-space models are only two-dimensional (output x input).'); 
                    end
            end
        end
        
        function varargout = subsref(obj,strct)
        % SUBSREF Overloads indexing of gsss objects
            if strcmp(strct(1).type,'()')
                siz = size(obj);
                if strct(1).subs{1}==':'; strct(1).subs{1}=1:siz(1);end
                if strct(1).subs{2}==':'; strct(1).subs{2}=1:siz(2);end
                assert(all(cellfun(@isnumeric,strct(1).subs)),'Only numeric indices are allowed.');
                assert(all(strct(1).subs{1}<=siz(1)) && all(strct(1).subs{2}<=siz(2)),'Index exceeds system dimensions.');
                varargout = {submodel(obj,strct(1).subs{1},strct(1).subs{2})};
                if length(strct)>1
                    varargout = {builtin('subsref',varargout{1},strct(2:end))};
                end
            else
                varargout = {builtin('subsref',obj,strct(:))};
            end
        end
        
        function sm = submodel(obj,idxout,idxin)
        % SUBMODEL Helper for subsref
            sm = gssshinfcd.gsss(obj.A,obj.B(:,idxin),obj.C(idxout,:),obj.D(idxout,idxin),obj.Ts,obj.parameters);
        end
        
        function idx = end(obj,k,~)
        % END Overloads the end operator of MATLAB
            switch k
                case 1; idx = size(obj.D,1);
                case 2; idx = size(obj.D,2); 
                otherwise; error('Gain-scheduled state-space models are only two-dimensional (output x input).');  
            end
        end
        
        function model = tolpvdss(obj)
        % TOLPVDSS To LCToolbox LPVDSSmod
            params = cell(1,length(obj.parameters));
            for k=1:length(obj.parameters)
                params{k} = SchedulingParameter(obj.parameters(k).name,obj.parameters(k).range,obj.parameters(k).rate);
            end
            model = LPVDSSmod(obj.A,obj.B,obj.C,obj.D,eye(size(obj.A)),params,obj.Ts);
        end
        
    end
    
    methods (Static)
        
        function obj = fromlpvdss(model)
            params = model.parameters;
            parameters = struct();
            for k=1:length(params)
                parameters(k).name = params{k}.argument;
                parameters(k).range = params{k}.range;
                parameters(k).rate = params{k}.rate;
            end
            obj = gssshinfcd.gsss(model.A,model.B,model.C,model.D,model.Ts,parameters);
        end
        
    end
end

