classdef  SharedObject < splines.PrintSharedObject
    %SHAREDOBJECT 
    %
    %
    %
  methods
    function varargout = to_string(self,varargin)
    %TO_STRING 
    %
    %  char = TO_STRING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(29, self, varargin{:});
    end

    function [] = disp(self)
      disp(self.to_string())
    end
    function self = subsasgn(self,varargin)
        error('Not supported: subsasgn');
    end
      function self = SharedObject(varargin)
    %SHAREDOBJECT 
    %
    %  new_obj = SHAREDOBJECT()
    %
    %
      self@splines.PrintSharedObject(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(30, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(31, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
