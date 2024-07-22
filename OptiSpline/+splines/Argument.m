classdef  Argument < splines.SharedObject
    %ARGUMENT 
    %
    %
    %
  methods
    function varargout = type(self,varargin)
    %TYPE 
    %
    %  char = TYPE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(33, self, varargin{:});
    end
    function varargout = is_all(self,varargin)
    %IS_ALL 
    %
    %  bool = IS_ALL(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(34, self, varargin{:});
    end
    function self = Argument(varargin)
    %ARGUMENT 
    %
    %  new_obj = ARGUMENT()
    %  new_obj = ARGUMENT(Parameter const & parameter)
    %  new_obj = ARGUMENT(int index)
    %  new_obj = ARGUMENT(char name)
    %
    %
      self@splines.SharedObject(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(36, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(37, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
    function varargout = from_vector(varargin)
    %FROM_VECTOR 
    %
    %  [index] = FROM_VECTOR(std::vector< Parameter,std::allocator< spline::Parameter > > const & ind)
    %  [index] = FROM_VECTOR([char] ind)
    %  [index] = FROM_VECTOR([int] ind)
    %
    %
     [varargout{1:nargout}] = splinesMEX(32, varargin{:});
    end
    function varargout = concrete(varargin)
    %CONCRETE 
    %
    %  int = CONCRETE(self, [char] args)
    %  [int] = CONCRETE([index] args, TensorBasis const & tb)
    %  [int] = CONCRETE([index] args, [char] strings)
    %  int = CONCRETE(index arg, TensorBasis const & tb)
    %
    %
     [varargout{1:nargout}] = splinesMEX(35, varargin{:});
    end
  end
end
