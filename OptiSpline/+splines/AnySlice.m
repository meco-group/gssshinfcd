classdef  AnySlice < SwigRef
    %ANYSLICE 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function v = start_(self)
      v = splinesMEX(39, self);
    end
    function v = stop_(self)
      v = splinesMEX(40, self);
    end
    function v = step_(self)
      v = splinesMEX(41, self);
    end
    function v = indices_(self)
      v = splinesMEX(42, self);
    end
    function v = has_indices_(self)
      v = splinesMEX(43, self);
    end
    function varargout = indices(self,varargin)
    %INDICES 
    %
    %  [int] = INDICES(self, int len)
    %
    %
      [varargout{1:nargout}] = splinesMEX(44, self, varargin{:});
    end
    function self = AnySlice(varargin)
    %ANYSLICE 
    %
    %  new_obj = ANYSLICE()
    %  new_obj = ANYSLICE(int i)
    %  new_obj = ANYSLICE([int] indices)
    %  new_obj = ANYSLICE(int start, int stop, int step)
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(45, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(46, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
