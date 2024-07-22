classdef  OptiSplineSol < casadi.OptiSol
    %OPTISPLINESOL 
    %
    %
    %
  methods
    function varargout = value(self,varargin)
    %VALUE Obtain value of expression at the current value
    %
    %  double = VALUE(self, DM x, {MX} values)
    %  double = VALUE(self, SX x, {MX} values)
    %  double = VALUE(self, MX x, {MX} values)
    %  Coefficient = VALUE(self, Coefficient c)
    %    
    %  Function = VALUE(self, Function f)
    %    
    %
    %
    %In regular mode, teh current value is the converged solution In debug mode,
    %the value can be non-converged
    %
    %Parameters:
    %-----------
    %
    %values:  Optional assignment expressions (e.g. x==3) to overrule the current
    %value
    %
    %
    %> VALUE(self, Coefficient c)
    %> VALUE(self, Function f)
    %------------------------------------------------------------------------
    %
    %> VALUE(self, DM x, {MX} values)
    %> VALUE(self, SX x, {MX} values)
    %> VALUE(self, MX x, {MX} values)
    %------------------------------------------------------------------------
    %
    %
    %Obtain value of expression at the current value
    %
    %In regular mode, teh current value is the converged solution In debug mode,
    %the value can be non-converged
    %
    %Parameters:
    %-----------
    %
    %values:  Optional assignment expressions (e.g. x==3) to overrule the current
    %value
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(404, self, varargin{:});
    end
    function varargout = opti(self,varargin)
    %OPTI 
    %
    %  OptiSpline = OPTI(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(405, self, varargin{:});
    end
    function self = OptiSplineSol(varargin)
    %OPTISPLINESOL 
    %
    %  new_obj = OPTISPLINESOL(OptiSpline opti)
    %
    %
      self@casadi.OptiSol(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(406, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(407, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
