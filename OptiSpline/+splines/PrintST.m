classdef  PrintST < SwigRef
    %PRINTST 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function self = PrintST(varargin)
    %PRINTST 
    %
    %  new_obj = PRINTST()
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(27, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(28, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
