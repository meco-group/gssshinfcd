classdef  Parameter < splines.Polynomial
    %PARAMETER 
    %
    %
    %
  methods
    function varargout = type(self,varargin)
    %TYPE 
    %
    %  char = TYPE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(369, self, varargin{:});
    end
    function varargout = to_string(self,varargin)
    %TO_STRING 
    %
    %  char = TO_STRING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(370, self, varargin{:});
    end
    function varargout = name(self,varargin)
    %NAME 
    %
    %  char = NAME(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(371, self, varargin{:});
    end
    function self = Parameter(varargin)
    %PARAMETER 
    %
    %  new_obj = PARAMETER()
    %  new_obj = PARAMETER(char a)
    %
    %
      self@splines.Polynomial(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(372, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(373, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
