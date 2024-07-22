classdef  BSplineBasis < splines.UnivariateBasis
    %BSPLINEBASIS 
    %
    %
    %
  methods
    function varargout = knots(self,varargin)
    %KNOTS 
    %
    %  [AnyScalar] = KNOTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(242, self, varargin{:});
    end
    function varargout = greville(self,varargin)
    %GREVILLE 
    %
    %  [AnyScalar] = GREVILLE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(243, self, varargin{:});
    end
    function self = BSplineBasis(varargin)
    %BSPLINEBASIS 
    %
    %  new_obj = BSPLINEBASIS(AnyVector knots, int degree)
    %  new_obj = BSPLINEBASIS(AnyVector bounds, int degree, int numberOfIntervals)
    %
    %
      self@splines.UnivariateBasis(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(244, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(245, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
    function varargout = from_single(varargin)
    %FROM_SINGLE 
    %
    %  BSplineBasis = FROM_SINGLE(AnyVector knots, int degree)
    %
    %
     [varargout{1:nargout}] = splinesMEX(241, varargin{:});
    end
  end
end
