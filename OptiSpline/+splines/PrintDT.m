classdef  PrintDT < SwigRef
    %PRINTDT 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function self = PrintDT(varargin)
    %PRINTDT 
    %
    %  new_obj = PRINTDT()
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(25, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(26, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
