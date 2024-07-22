classdef  OptiSpline < casadi.Opti & splines.OptiSplineI
    %OPTISPLINE 
    %
    %
    %
  methods
    function this = swig_this(self)
      this = splinesMEX(3, self);
    end
    function varargout = copy(self,varargin)
    %COPY Get a copy of the.
    %
    %  OptiSpline = COPY(self)
    %
    %
    %The copy is effectively a deep copy: Updating the state of the copy does not
    %update the original.
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(386, self, varargin{:});
    end
    function varargout = debug(self,varargin)
    %DEBUG Get a copy with advanced functionality.
    %
    %  OptiSplineAdvanced = DEBUG(self)
    %
    %
    %You get access to more methods, but you have no guarantees about API
    %stability
    %
    %The copy is effectively a deep copy: Updating the state of the copy does not
    %update the original.
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(387, self, varargin{:});
    end
    function varargout = advanced(self,varargin)
    %ADVANCED Get a copy with advanced functionality.
    %
    %  OptiSplineAdvanced = ADVANCED(self)
    %
    %
    %You get access to more methods, but you have no guarantees about API
    %stability
    %
    %The copy is effectively a deep copy: Updating the state of the copy does not
    %update the original.
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(388, self, varargin{:});
    end
    function varargout = solve(self,varargin)
    %SOLVE Crunch the numbers; solve the problem.
    %
    %  OptiSplineSol = SOLVE(self)
    %
    %
    %
    %
      [varargout{1:nargout}] = splinesMEX(389, self, varargin{:});
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
      [varargout{1:nargout}] = splinesMEX(390, self, varargin{:});
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
      [varargout{1:nargout}] = splinesMEX(391, self, varargin{:});
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
      [varargout{1:nargout}] = splinesMEX(392, self, varargin{:});
    end

      function out = variable(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', st(1));
        else
          meta = struct;
        end
        out = self.internal_variable(varargin{:});
        self.update_user_dict(out, meta);
      end
      function out = parameter(self, varargin)
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', st(1));
        else
          meta = struct;
        end
        out = self.internal_parameter(varargin{:});
        self.update_user_dict(out, meta);
      end
      function [] = subject_to(self, varargin)
        if length(varargin)==0
          self.internal_subject_to();
          return;
        end
        st = dbstack('-completenames',1);
        if length(st)>0
          meta = struct('stacktrace', st(1));
        else
          meta = struct;
        end
        self.internal_subject_to(varargin{:});
        self.update_user_dict(varargin{1}, meta);
      end
    
    function [] = callback(self, varargin)
      casadi.OptiCallbackHelper.callback_setup(self, varargin{:})
    end
      function self = OptiSpline(varargin)
    %OPTISPLINE 
    %
    %  new_obj = OPTISPLINE()
    %  new_obj = OPTISPLINE(OptiSplineAdvanced rhs)
    %
    %
      self@casadi.Opti(SwigRef.Null);
      self@splines.OptiSplineI(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(393, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(394, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
