classdef  MonomialBasis < splines.UnivariateBasis
    %MONOMIALBASIS 
    %
    %
    %
  methods
    function self = MonomialBasis(varargin)
    %MONOMIALBASIS 
    %
    %  new_obj = MONOMIALBASIS(int degree)
    %
    %
      self@splines.UnivariateBasis(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(239, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(240, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
