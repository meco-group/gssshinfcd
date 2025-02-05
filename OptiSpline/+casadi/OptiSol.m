classdef  OptiSol < casadi.PrintableCommon
    %OPTISOL A simplified interface for NLP modeling/solving.
    %
    %
    %
    %This class offers a view with solution retrieval facilities The API is
    %guaranteed to be stable.
    %
    %Joris Gillis, Erik Lambrechts
    %
    %C++ includes: optistack.hpp 
    %
  methods
    function varargout = type_name(self,varargin)
    %TYPE_NAME 
    %
    %  char = TYPE_NAME(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1219, self, varargin{:});
    end
    function varargout = disp(self,varargin)
    %DISP 
    %
    %  std::ostream & = DISP(self, bool more)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1220, self, varargin{:});
    end
    function varargout = str(self,varargin)
    %STR 
    %
    %  char = STR(self, bool more)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1221, self, varargin{:});
    end
    function varargout = value(self,varargin)
    %VALUE Obtain value of expression at the current value
    %
    %  double = VALUE(self, DM x, {MX} values)
    %  double = VALUE(self, SX x, {MX} values)
    %  double = VALUE(self, MX x, {MX} values)
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
    %
      [varargout{1:nargout}] = casadiMEX(1222, self, varargin{:});
    end
    function varargout = value_variables(self,varargin)
    %VALUE_VARIABLES get assignment expressions for the optimal solution
    %
    %  {MX} = VALUE_VARIABLES(self)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1223, self, varargin{:});
    end
    function varargout = value_parameters(self,varargin)
    %VALUE_PARAMETERS 
    %
    %  {MX} = VALUE_PARAMETERS(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1224, self, varargin{:});
    end
    function varargout = stats(self,varargin)
    %STATS Get statistics.
    %
    %  struct = STATS(self)
    %
    %
    %nlpsol stats are passed as-is. No stability can be guaranteed about this
    %part of the API
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(1225, self, varargin{:});
    end
    function varargout = opti(self,varargin)
    %OPTI 
    %
    %  Opti = OPTI(self)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1226, self, varargin{:});
    end
    function self = OptiSol(varargin)
    %OPTISOL 
    %
    %  new_obj = OPTISOL()
    %
    %
      self@casadi.PrintableCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(1227, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        casadiMEX(1228, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
