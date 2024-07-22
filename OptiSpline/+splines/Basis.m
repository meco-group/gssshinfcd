classdef  Basis < splines.SharedObject
    %BASIS 
    %
    %
    %
  methods
    function varargout = plus(self,varargin)
    %PLUS 
    %
    %  Basis = PLUS(self, Basis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(213, self, varargin{:});
    end
    function varargout = mtimes(self,varargin)
    %MTIMES 
    %
    %  Basis = MTIMES(self, Basis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(214, self, varargin{:});
    end
    function varargout = eval(self,varargin)
    %EVAL 
    %
    %  AnyTensor = EVAL(self, AnyTensor x)
    %
    %
      [varargout{1:nargout}] = splinesMEX(215, self, varargin{:});
    end
    function varargout = list_eval(self,varargin)
    %LIST_EVAL 
    %
    %  AnyTensor = LIST_EVAL(self, AnyTensor x)
    %
    %
      [varargout{1:nargout}] = splinesMEX(216, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  bool = EQ(self, Basis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(217, self, varargin{:});
    end
    function varargout = type(self,varargin)
    %TYPE 
    %
    %  char = TYPE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(218, self, varargin{:});
    end
    function varargout = domain(self,varargin)
    %DOMAIN 
    %
    %  Domain = DOMAIN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(219, self, varargin{:});
    end
    function varargout = evaluation_grid(self,varargin)
    %EVALUATION_GRID 
    %
    %  AnyTensor = EVALUATION_GRID(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(220, self, varargin{:});
    end
    function varargout = insert_knots(self,varargin)
    %INSERT_KNOTS 
    %
    %  [Basis , AnyTensor OUTPUT] = INSERT_KNOTS(self, AnyVector new_knots)
    %
    %
      [varargout{1:nargout}] = splinesMEX(221, self, varargin{:});
    end
    function varargout = midpoint_refinement(self,varargin)
    %MIDPOINT_REFINEMENT 
    %
    %  [Basis , AnyTensor OUTPUT] = MIDPOINT_REFINEMENT(self, int refinement)
    %
    %
      [varargout{1:nargout}] = splinesMEX(222, self, varargin{:});
    end
    function varargout = degree_elevation(self,varargin)
    %DEGREE_ELEVATION 
    %
    %  [Basis , AnyTensor OUTPUT] = DEGREE_ELEVATION(self, int elevation)
    %
    %
      [varargout{1:nargout}] = splinesMEX(223, self, varargin{:});
    end
    function varargout = kick_boundary(self,varargin)
    %KICK_BOUNDARY 
    %
    %  [Basis , AnyTensor OUTPUT] = KICK_BOUNDARY(self, Domain boundary)
    %
    %
      [varargout{1:nargout}] = splinesMEX(224, self, varargin{:});
    end
    function varargout = project_to(self,varargin)
    %PROJECT_TO 
    %
    %  AnyTensor = PROJECT_TO(self, Basis b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(225, self, varargin{:});
    end
    function varargout = transform_to(self,varargin)
    %TRANSFORM_TO 
    %
    %  [Basis , AnyTensor OUTPUT] = TRANSFORM_TO(self, Basis b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(226, self, varargin{:});
    end
    function varargout = basis_functions(self,varargin)
    %BASIS_FUNCTIONS 
    %
    %  Function = BASIS_FUNCTIONS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(227, self, varargin{:});
    end
    function varargout = dimension(self,varargin)
    %DIMENSION 
    %
    %  int = DIMENSION(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(228, self, varargin{:});
    end
    function varargout = n_inputs(self,varargin)
    %N_INPUTS 
    %
    %  int = N_INPUTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(229, self, varargin{:});
    end
    function varargout = derivative(self,varargin)
    %DERIVATIVE 
    %
    %  [Basis , AnyTensor OUTPUT] = DERIVATIVE(self)
    %  [Basis , AnyTensor OUTPUT] = DERIVATIVE(self, int order)
    %
    %
      [varargout{1:nargout}] = splinesMEX(230, self, varargin{:});
    end
    function varargout = antiderivative(self,varargin)
    %ANTIDERIVATIVE 
    %
    %  [Basis , AnyTensor OUTPUT] = ANTIDERIVATIVE(self)
    %  [Basis , AnyTensor OUTPUT] = ANTIDERIVATIVE(self, int order)
    %
    %
      [varargout{1:nargout}] = splinesMEX(231, self, varargin{:});
    end
    function varargout = integral(self,varargin)
    %INTEGRAL 
    %
    %  AnyTensor = INTEGRAL(self, Domain domain)
    %
    %
      [varargout{1:nargout}] = splinesMEX(232, self, varargin{:});
    end
    function varargout = const_coeff_tensor(self,varargin)
    %CONST_COEFF_TENSOR 
    %
    %  AnyTensor = CONST_COEFF_TENSOR(self, AnyTensor t)
    %
    %
      [varargout{1:nargout}] = splinesMEX(233, self, varargin{:});
    end

      function [varargout] = horzcat(self,varargin)
        error('Cannot concatenate Basis objects; perhaps you meant to group them with {}?');
      end
      function [varargout] = vertcat(self,varargin)
        error('Cannot concatenate Basis objects; perhaps you meant to group them with {}?');
      end
       function self = Basis(varargin)
    %BASIS 
    %
    %  new_obj = BASIS()
    %
    %
      self@splines.SharedObject(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(234, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(235, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
