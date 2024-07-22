classdef  UnivariateBasis < splines.Basis
    %UNIVARIATEBASIS 
    %
    %
    %
  methods
    function varargout = degree(self,varargin)
    %DEGREE 
    %
    %  int = DEGREE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(236, self, varargin{:});
    end
    function self = UnivariateBasis(varargin)
    %UNIVARIATEBASIS 
    %
    %  new_obj = UNIVARIATEBASIS()
    %
    %
      self@splines.Basis(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(237, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(238, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
