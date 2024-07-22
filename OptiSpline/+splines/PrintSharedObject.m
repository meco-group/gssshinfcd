classdef  PrintSharedObject < SwigRef
    %PRINTSHAREDOBJECT 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function self = PrintSharedObject(varargin)
    %PRINTSHAREDOBJECT 
    %
    %  new_obj = PRINTSHAREDOBJECT()
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(21, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(22, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
