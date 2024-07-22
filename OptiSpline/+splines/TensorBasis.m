classdef  TensorBasis < splines.SharedObject
    %TENSORBASIS 
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
      [varargout{1:nargout}] = splinesMEX(246, self, varargin{:});
    end
    function varargout = n_basis(self,varargin)
    %N_BASIS 
    %
    %  int = N_BASIS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(247, self, varargin{:});
    end
    function varargout = n_inputs(self,varargin)
    %N_INPUTS 
    %
    %  int = N_INPUTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(248, self, varargin{:});
    end
    function varargout = n_inputs_list(self,varargin)
    %N_INPUTS_LIST 
    %
    %  [int] = N_INPUTS_LIST(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(249, self, varargin{:});
    end
    function varargout = dimension(self,varargin)
    %DIMENSION 
    %
    %  [int] = DIMENSION(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(250, self, varargin{:});
    end
    function varargout = arguments(self,varargin)
    %ARGUMENTS 
    %
    %  [char] = ARGUMENTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(251, self, varargin{:});
    end
    function varargout = arguments_temp(self,varargin)
    %ARGUMENTS_TEMP 
    %
    %  [index] = ARGUMENTS_TEMP(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(252, self, varargin{:});
    end
    function varargout = argument(self,varargin)
    %ARGUMENT 
    %
    %  char = ARGUMENT(self, int index)
    %
    %
      [varargout{1:nargout}] = splinesMEX(253, self, varargin{:});
    end
    function varargout = index_argument(self,varargin)
    %INDEX_ARGUMENT 
    %
    %  int = INDEX_ARGUMENT(self, index a)
    %
    %
      [varargout{1:nargout}] = splinesMEX(254, self, varargin{:});
    end
    function varargout = indexArgument(self,varargin)
    %INDEXARGUMENT 
    %
    %  int = INDEXARGUMENT(self, char a)
    %
    %
      [varargout{1:nargout}] = splinesMEX(255, self, varargin{:});
    end
    function varargout = hasArguments(self,varargin)
    %HASARGUMENTS 
    %
    %  bool = HASARGUMENTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(256, self, varargin{:});
    end
    function varargout = valid_argument(self,varargin)
    %VALID_ARGUMENT 
    %
    %  bool = VALID_ARGUMENT(self, index a)
    %
    %
      [varargout{1:nargout}] = splinesMEX(257, self, varargin{:});
    end
    function varargout = valid_argument_list(self,varargin)
    %VALID_ARGUMENT_LIST 
    %
    %  bool = VALID_ARGUMENT_LIST(self, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(258, self, varargin{:});
    end
    function varargout = domain(self,varargin)
    %DOMAIN 
    %
    %  TensorDomain = DOMAIN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(259, self, varargin{:});
    end
    function varargout = basis(self,varargin)
    %BASIS 
    %
    %  Basis = BASIS(self)
    %  Basis = BASIS(self, index index)
    %  Basis = BASIS(self, char a)
    %
    %
      [varargout{1:nargout}] = splinesMEX(260, self, varargin{:});
    end
    function varargout = bases(self,varargin)
    %BASES 
    %
    %  [Basis] = BASES(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(261, self, varargin{:});
    end
    function varargout = add_basis(self,varargin)
    %ADD_BASIS 
    %
    %  TensorBasis = ADD_BASIS(self, TensorBasis basis)
    %  TensorBasis = ADD_BASIS(self, Basis basis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(262, self, varargin{:});
    end
    function varargout = substitute_bases(self,varargin)
    %SUBSTITUTE_BASES 
    %
    %  TensorBasis = SUBSTITUTE_BASES(self, [index] arg_ind, [Basis] bases)
    %
    %
      [varargout{1:nargout}] = splinesMEX(263, self, varargin{:});
    end
    function varargout = plus(self,varargin)
    %PLUS 
    %
    %  TensorBasis = PLUS(self, TensorBasis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(264, self, varargin{:});
    end
    function varargout = mtimes(self,varargin)
    %MTIMES 
    %
    %  TensorBasis = MTIMES(self, TensorBasis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(265, self, varargin{:});
    end
    function varargout = eval(self,varargin)
    %EVAL 
    %
    %  AnyTensor = EVAL(self, [AnyScalar] x, [index] arg_ind, bool reorder_output)
    %
    %
      [varargout{1:nargout}] = splinesMEX(266, self, varargin);
    end
    function varargout = grid_eval(self,varargin)
    %GRID_EVAL 
    %
    %  AnyTensor = GRID_EVAL(self, [AnyTensor] x, [index] arg_ind, bool reorder_output)
    %
    %
      [varargout{1:nargout}] = splinesMEX(267, self, varargin{:});
    end
    function varargout = evaluation_grid(self,varargin)
    %EVALUATION_GRID 
    %
    %  [AnyTensor] = EVALUATION_GRID(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(268, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  bool = EQ(self, TensorBasis rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(269, self, varargin{:});
    end
    function varargout = totalNumberBasisFunctions(self,varargin)
    %TOTALNUMBERBASISFUNCTIONS 
    %
    %  int = TOTALNUMBERBASISFUNCTIONS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(270, self, varargin{:});
    end
    function varargout = const_coeff_tensor(self,varargin)
    %CONST_COEFF_TENSOR 
    %
    %  AnyTensor = CONST_COEFF_TENSOR(self, AnyTensor t)
    %
    %
      [varargout{1:nargout}] = splinesMEX(271, self, varargin{:});
    end
    function varargout = insert_knots(self,varargin)
    %INSERT_KNOTS 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = INSERT_KNOTS(self, [AnyVector] new_knots, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(272, self, varargin{:});
    end
    function varargout = midpoint_refinement(self,varargin)
    %MIDPOINT_REFINEMENT 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = MIDPOINT_REFINEMENT(self, [int] refinement, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(273, self, varargin{:});
    end
    function varargout = degree_elevation(self,varargin)
    %DEGREE_ELEVATION 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = DEGREE_ELEVATION(self, [int] elevation, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(274, self, varargin{:});
    end
    function varargout = kick_boundary(self,varargin)
    %KICK_BOUNDARY 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = KICK_BOUNDARY(self, TensorDomain boundary, [char] args)
    %  [TensorBasis , [AnyTensor] OUTPUT] = KICK_BOUNDARY(self, TensorDomain boundary, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(275, self, varargin{:});
    end
    function varargout = derivative(self,varargin)
    %DERIVATIVE 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = DERIVATIVE(self, [index] arg)
    %  [TensorBasis , [AnyTensor] OUTPUT] = DERIVATIVE(self, [int] order, [index] arg)
    %
    %
      [varargout{1:nargout}] = splinesMEX(276, self, varargin{:});
    end
    function varargout = antiderivative(self,varargin)
    %ANTIDERIVATIVE 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = ANTIDERIVATIVE(self, [index] arg)
    %  [TensorBasis , [AnyTensor] OUTPUT] = ANTIDERIVATIVE(self, [int] order, [index] arg)
    %
    %
      [varargout{1:nargout}] = splinesMEX(277, self, varargin{:});
    end
    function varargout = integral(self,varargin)
    %INTEGRAL 
    %
    %  [AnyTensor] = INTEGRAL(self, TensorDomain domain)
    %
    %
      [varargout{1:nargout}] = splinesMEX(278, self, varargin{:});
    end
    function varargout = partial_integral(self,varargin)
    %PARTIAL_INTEGRAL 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = PARTIAL_INTEGRAL(self, TensorDomain domain, [char] args)
    %  [TensorBasis , [AnyTensor] OUTPUT] = PARTIAL_INTEGRAL(self, TensorDomain domain, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(279, self, varargin{:});
    end
    function varargout = project_to(self,varargin)
    %PROJECT_TO 
    %
    %  AnyTensor = PROJECT_TO(self, TensorBasis b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(280, self, varargin{:});
    end
    function varargout = transform_to(self,varargin)
    %TRANSFORM_TO 
    %
    %  [TensorBasis , [AnyTensor] OUTPUT] = TRANSFORM_TO(self, TensorBasis b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(281, self, varargin{:});
    end
    function varargout = basis_functions(self,varargin)
    %BASIS_FUNCTIONS 
    %
    %  Function = BASIS_FUNCTIONS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(282, self, varargin{:});
    end
    function varargout = get_permutation(self,varargin)
    %GET_PERMUTATION 
    %
    %  [int] = GET_PERMUTATION(self, TensorBasis grid)
    %
    %
      [varargout{1:nargout}] = splinesMEX(283, self, varargin{:});
    end
    function self = TensorBasis(varargin)
    %TENSORBASIS 
    %
    %  new_obj = TENSORBASIS()
    %  new_obj = TENSORBASIS(Basis bases_)
    %  new_obj = TENSORBASIS([Basis] bases_)
    %  new_obj = TENSORBASIS([TensorBasis] allBasis)
    %  new_obj = TENSORBASIS(TensorBasis tensor, [char] args)
    %  new_obj = TENSORBASIS(Basis bases_, char args)
    %  new_obj = TENSORBASIS([Basis] bases_, [char] args)
    %
    %
      self@splines.SharedObject(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(284, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(285, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
  end
end
