classdef  Interval < splines.Domain
    %INTERVAL 
    %
    %
    %
  methods
    function varargout = min(self,varargin)
    %MIN 
    %
    %  AnyScalar = MIN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(197, self, varargin{:});
    end
    function varargout = max(self,varargin)
    %MAX 
    %
    %  AnyScalar = MAX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(198, self, varargin{:});
    end
    function self = Interval(varargin)
    %INTERVAL 
    %
    %  new_obj = INTERVAL()
    %  new_obj = INTERVAL([AnyScalar] data)
    %  new_obj = INTERVAL(Domain dom)
    %  new_obj = INTERVAL(AnyScalar min, AnyScalar max)
    %
    %
      self@splines.Domain(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(199, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(200, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
