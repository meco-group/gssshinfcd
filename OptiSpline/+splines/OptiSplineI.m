classdef  OptiSplineI < SwigRef
    %OPTISPLINEI 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function varargout = set_initial(self,varargin)
    %SET_INITIAL 
    %
    %  SET_INITIAL(self, MTensor t, DTensor d)
    %  SET_INITIAL(self, Coefficient c, DTensor d)
    %  SET_INITIAL(self, Function f, Function g)
    %
    %
      [varargout{1:nargout}] = splinesMEX(374, self, varargin{:});
    end
    function varargout = set_value(self,varargin)
    %SET_VALUE 
    %
    %  SET_VALUE(self, MTensor t, DTensor d)
    %  SET_VALUE(self, Coefficient c, DTensor d)
    %  SET_VALUE(self, Function f, Function g)
    %
    %
      [varargout{1:nargout}] = splinesMEX(375, self, varargin{:});
    end
    function varargout = update_user_dict(self,varargin)
    %UPDATE_USER_DICT 
    %
    %  UPDATE_USER_DICT(self, MTensor t, struct meta)
    %  UPDATE_USER_DICT(self, Function f, struct meta)
    %  UPDATE_USER_DICT(self, Coefficient c, struct meta)
    %
    %
      [varargout{1:nargout}] = splinesMEX(376, self, varargin{:});
    end
    function varargout = Function(self,varargin)
    %FUNCTION 
    %
    %  Function = FUNCTION(self, TensorBasis b, [int] shape, char attribute, char opti_type)
    %
    %
      [varargout{1:nargout}] = splinesMEX(377, self, varargin{:});
    end
    function self = OptiSplineI(varargin)
    %OPTISPLINEI 
    %
    %  new_obj = OPTISPLINEI()
    %
    %
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(378, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(379, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
