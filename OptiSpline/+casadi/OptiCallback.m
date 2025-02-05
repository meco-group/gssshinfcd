classdef  OptiCallback < SwigRef
    %OPTICALLBACK C++ includes: optistack.hpp 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = casadiMEX(3, self);
    end
    function self = OptiCallback(varargin)
    %OPTICALLBACK 
    %
    %  new_obj = OPTICALLBACK(self)
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        if strcmp(class(self),'director_basic.OptiCallback')
          tmp = casadiMEX(1177, 0, varargin{:});
        else
          tmp = casadiMEX(1177, self, varargin{:});
        end
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function varargout = call(self,varargin)
    %CALL 
    %
    %  CALL(self, int i)
    %
    %
      [varargout{1:nargout}] = casadiMEX(1178, self, varargin{:});
    end
    function delete(self)
      if self.swigPtr
        casadiMEX(1179, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
