classdef  OptiSplineAdvanced < casadi.OptiAdvanced & splines.OptiSplineAdvancedI
    %OPTISPLINEADVANCED 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function varargout = solve(self,varargin)
    %SOLVE Crunch the numbers; solve the problem.
    %
    %  OptiSplineSol = SOLVE(self)
    %
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(395, self, varargin{:});
    end
    function varargout = value(self,varargin)
    %VALUE Obtain value of expression at the current value
    %
    %  double = VALUE(self, DM x, {MX} values)
    %  double = VALUE(self, SX x, {MX} values)
    %  AnyTensor = VALUE(self, AnyTensor c)
    %    
    %  double = VALUE(self, MX x, {MX} values)
    %  Coefficient = VALUE(self, Coefficient c)
    %    
    %  Function = VALUE(self, Function f)
    %    
    %
    %
    %In regular mode, teh current value is the converged solution In debug mode,
    %the value can be non-converged
    %
    %Parameters:
    %-----------
    %
    %values:  Optional assignment expressions (e.g. x==3) to overrule the current
    %value
    %
    %
    %> VALUE(self, AnyTensor c)
    %> VALUE(self, Coefficient c)
    %> VALUE(self, Function f)
    %------------------------------------------------------------------------
    %
    %> VALUE(self, DM x, {MX} values)
    %> VALUE(self, SX x, {MX} values)
    %> VALUE(self, MX x, {MX} values)
    %------------------------------------------------------------------------
    %
    %
    %Obtain value of expression at the current value
    %
    %In regular mode, teh current value is the converged solution In debug mode,
    %the value can be non-converged
    %
    %Parameters:
    %-----------
    %
    %values:  Optional assignment expressions (e.g. x==3) to overrule the current
    %value
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(396, self, varargin{:});
    end
    function varargout = symvar(self,varargin)
    %SYMVAR Get symbols present in expression.
    %
    %  {MX} = SYMVAR(self)
    %  {MX} = SYMVAR(self, MTensor m)
    %    
    %  {MX} = SYMVAR(self, MX expr)
    %  {MX} = SYMVAR(self, MTensor expr, casadi::VariableType type)
    %    
    %  {MX} = SYMVAR(self, MX expr, casadi::VariableType type)
    %
    %
    %Returned vector is ordered according to the order of variable()/parameter()
    %calls used to create the variables
    %
    %
    %> SYMVAR(self, MTensor m)
    %> SYMVAR(self, MTensor expr, casadi::VariableType type)
    %------------------------------------------------------------------------
    %
    %> SYMVAR(self)
    %> SYMVAR(self, MX expr)
    %> SYMVAR(self, MX expr, casadi::VariableType type)
    %------------------------------------------------------------------------
    %
    %
    %Get symbols present in expression.
    %
    %Returned vector is ordered according to the order of variable()/parameter()
    %calls used to create the variables
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(397, self, varargin{:});
    end
    function varargout = matlab_dump(self,varargin)
    %MATLAB_DUMP 
    %
    %  MATLAB_DUMP(self, Function f, char fname)
    %
    %
      [varargout{1:nargout}] = splinesMEX(398, self, varargin{:});
    end
    function varargout = set_initial(self,varargin)
    %SET_INITIAL Set initial guess for decision variables
    %
    %  SET_INITIAL(self, {MX} assignments)
    %  SET_INITIAL(self, MX x, DM v)
    %  SET_INITIAL(self, MTensor t, DTensor d)
    %    
    %  SET_INITIAL(self, Coefficient c, DTensor d)
    %    
    %  SET_INITIAL(self, Function f, Function g)
    %    
    %
    %
    %::
    %
    %  * opti.set_initial(x, 2)
    %  * opti.set_initial(10*x(1), 2)
    %  * 
    %
    %
    %
    %
    %> SET_INITIAL(self, MTensor t, DTensor d)
    %> SET_INITIAL(self, Coefficient c, DTensor d)
    %> SET_INITIAL(self, Function f, Function g)
    %------------------------------------------------------------------------
    %
    %> SET_INITIAL(self, {MX} assignments)
    %> SET_INITIAL(self, MX x, DM v)
    %------------------------------------------------------------------------
    %
    %
    %Set initial guess for decision variables
    %
    %::
    %
    %  * opti.set_initial(x, 2)
    %  * opti.set_initial(10*x(1), 2)
    %  * 
    %
    %
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(399, self, varargin{:});
    end
    function varargout = set_value(self,varargin)
    %SET_VALUE Set value of parameter.
    %
    %  SET_VALUE(self, {MX} assignments)
    %  SET_VALUE(self, MX x, DM v)
    %  SET_VALUE(self, MTensor t, DTensor d)
    %    
    %  SET_VALUE(self, Coefficient c, DTensor d)
    %    
    %  SET_VALUE(self, Function f, Function g)
    %    
    %
    %
    %Each parameter must be given a value before 'solve' can be called
    %
    %
    %> SET_VALUE(self, MTensor t, DTensor d)
    %> SET_VALUE(self, Coefficient c, DTensor d)
    %> SET_VALUE(self, Function f, Function g)
    %------------------------------------------------------------------------
    %
    %> SET_VALUE(self, {MX} assignments)
    %> SET_VALUE(self, MX x, DM v)
    %------------------------------------------------------------------------
    %
    %
    %Set value of parameter.
    %
    %Each parameter must be given a value before 'solve' can be called
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(400, self, varargin{:});
    end
    function varargout = update_user_dict(self,varargin)
    %UPDATE_USER_DICT 
    %
    %  UPDATE_USER_DICT(self, MTensor t, struct meta)
    %  UPDATE_USER_DICT(self, MX m, struct meta)
    %    add user data Add arbitrary data in the form of a dictionary to symbols or
    %  UPDATE_USER_DICT(self, {MX} m, struct meta)
    %  UPDATE_USER_DICT(self, Function f, struct meta)
    %  UPDATE_USER_DICT(self, Coefficient c, struct meta)
    %
    %> UPDATE_USER_DICT(self, MTensor t, struct meta)
    %> UPDATE_USER_DICT(self, Function f, struct meta)
    %> UPDATE_USER_DICT(self, Coefficient c, struct meta)
    %------------------------------------------------------------------------
    %
    %> UPDATE_USER_DICT(self, {MX} m, struct meta)
    %------------------------------------------------------------------------
    %
    %
    %
    %> UPDATE_USER_DICT(self, MX m, struct meta)
    %------------------------------------------------------------------------
    %
    %
    %add user data Add arbitrary data in the form of a dictionary to symbols or
    %constraints
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(401, self, varargin{:});
    end
    function self = OptiSplineAdvanced(varargin)
    %OPTISPLINEADVANCED 
    %
    %  new_obj = OPTISPLINEADVANCED(OptiSpline rhs)
    %
    %
      self@casadi.OptiAdvanced(SwigRef.Null);
      self@splines.OptiSplineAdvancedI(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(402, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(403, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
