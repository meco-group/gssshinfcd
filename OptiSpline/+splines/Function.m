classdef (InferiorClasses = {?casadi.DM,?casadi.SX,?casadi.MX,?DTensor,?STensor,?MTensor}) Function < splines.SharedObject
    %FUNCTION 
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
      [varargout{1:nargout}] = splinesMEX(305, self, varargin{:});
    end
    function varargout = ge(self,varargin)
    %GE 
    %
    %  MX = GE(self, MX x)
    %  MX = GE(self, Function f)
    %
    %
      [varargout{1:nargout}] = splinesMEX(306, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  MX = EQ(self, MX x)
    %  MX = EQ(self, Function f)
    %
    %
      [varargout{1:nargout}] = splinesMEX(307, self, varargin{:});
    end
    function varargout = eval_internal(self,varargin)
    %EVAL_INTERNAL 
    %
    %  AnyTensor = EVAL_INTERNAL(self, AnyTensor x, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(308, self, varargin{:});
    end
    function varargout = list_eval(self,varargin)
    %LIST_EVAL 
    %
    %  AnyTensor = LIST_EVAL(self, AnyTensor x, [index] args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(309, self, varargin{:});
    end
    function varargout = grid_eval(self,varargin)
    %GRID_EVAL 
    %
    %  AnyTensor = GRID_EVAL(self, [AnyTensor] x, [index] args, bool squeeze_return)
    %
    %
      [varargout{1:nargout}] = splinesMEX(310, self, varargin{:});
    end
    function varargout = partial_eval(self,varargin)
    %PARTIAL_EVAL 
    %
    %  Function = PARTIAL_EVAL(self, AnyTensor x, index args)
    %
    %
      [varargout{1:nargout}] = splinesMEX(311, self, varargin{:});
    end
    function varargout = uminus(self,varargin)
    %UMINUS 
    %
    %  Function = UMINUS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(312, self, varargin{:});
    end
    function varargout = pow(self,varargin)
    %POW 
    %
    %  Function = POW(self, int power)
    %
    %
      [varargout{1:nargout}] = splinesMEX(313, self, varargin{:});
    end
    function varargout = mpow(self,varargin)
    %MPOW 
    %
    %  Function = MPOW(self, int power)
    %
    %
      [varargout{1:nargout}] = splinesMEX(314, self, varargin{:});
    end
    function varargout = slice(self,varargin)
    %SLICE 
    %
    %  Function = SLICE(self, AnySlice i)
    %  Function = SLICE(self, AnySlice i, AnySlice j)
    %
    %
      [varargout{1:nargout}] = splinesMEX(315, self, varargin{:});
    end
    function varargout = transpose(self,varargin)
    %TRANSPOSE 
    %
    %  Function = TRANSPOSE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(316, self, varargin{:});
    end
    function varargout = trace(self,varargin)
    %TRACE 
    %
    %  Function = TRACE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(317, self, varargin{:});
    end
    function varargout = sum(self,varargin)
    %SUM 
    %
    %  Function = SUM(self)
    %  Function = SUM(self, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(318, self, varargin{:});
    end
    function varargout = vertcat(varargin)
    %VERTCAT 
    %
    %  Function = VERTCAT([Function] f)
    %
    %
     [varargout{1:nargout}] = splinesMEX(319, varargin);
    end
    function varargout = horzcat(varargin)
    %HORZCAT 
    %
    %  Function = HORZCAT([Function] f)
    %
    %
     [varargout{1:nargout}] = splinesMEX(320, varargin);
    end
    function varargout = blkdiag(varargin)
    %BLKDIAG 
    %
    %  Function = BLKDIAG([Function] f)
    %
    %
     [varargout{1:nargout}] = splinesMEX(323, varargin);
    end
    function varargout = coeff(self,varargin)
    %COEFF 
    %
    %  Coefficient = COEFF(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(324, self, varargin{:});
    end
    function varargout = coeff_tensor(self,varargin)
    %COEFF_TENSOR 
    %
    %  AnyTensor = COEFF_TENSOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(325, self, varargin{:});
    end
    function varargout = data(self,varargin)
    %DATA 
    %
    %  AnyTensor = DATA(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(326, self, varargin{:});
    end
    function varargout = basis(self,varargin)
    %BASIS 
    %
    %  Basis = BASIS(self)
    %  Basis = BASIS(self, index i)
    %
    %
      [varargout{1:nargout}] = splinesMEX(327, self, varargin{:});
    end
    function varargout = tensor_basis(self,varargin)
    %TENSOR_BASIS 
    %
    %  TensorBasis = TENSOR_BASIS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(328, self, varargin{:});
    end
    function varargout = domain(self,varargin)
    %DOMAIN 
    %
    %  TensorDomain = DOMAIN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(329, self, varargin{:});
    end
    function varargout = is_scalar(self,varargin)
    %IS_SCALAR 
    %
    %  bool = IS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(330, self, varargin{:});
    end
    function varargout = is_vector(self,varargin)
    %IS_VECTOR 
    %
    %  bool = IS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(331, self, varargin{:});
    end
    function varargout = is_column(self,varargin)
    %IS_COLUMN 
    %
    %  bool = IS_COLUMN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(332, self, varargin{:});
    end
    function varargout = shape(self,varargin)
    %SHAPE 
    %
    %  [int] = SHAPE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(333, self, varargin{:});
    end
    function varargout = reshape(self,varargin)
    %RESHAPE 
    %
    %  Function = RESHAPE(self, [int] shape)
    %
    %
      [varargout{1:nargout}] = splinesMEX(334, self, varargin{:});
    end
    function varargout = transform_to(self,varargin)
    %TRANSFORM_TO 
    %
    %  Function = TRANSFORM_TO(self, Basis basis)
    %  Function = TRANSFORM_TO(self, TensorBasis basis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(335, self, varargin{:});
    end
    function varargout = project_to(self,varargin)
    %PROJECT_TO 
    %
    %  Function = PROJECT_TO(self, Basis basis)
    %  Function = PROJECT_TO(self, TensorBasis basis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(336, self, varargin{:});
    end
    function varargout = n_inputs(self,varargin)
    %N_INPUTS 
    %
    %  int = N_INPUTS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(337, self, varargin{:});
    end
    function varargout = insert_knots(self,varargin)
    %INSERT_KNOTS 
    %
    %  Function = INSERT_KNOTS(self, AnyVector new_knots, index arg)
    %  Function = INSERT_KNOTS(self, [AnyVector] new_knots, [index] arg)
    %
    %
      [varargout{1:nargout}] = splinesMEX(338, self, varargin{:});
    end
    function varargout = midpoint_refinement(self,varargin)
    %MIDPOINT_REFINEMENT 
    %
    %  Function = MIDPOINT_REFINEMENT(self, int refinement, index arg)
    %  Function = MIDPOINT_REFINEMENT(self, [int] refinement, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(339, self, varargin{:});
    end
    function varargout = degree_elevation(self,varargin)
    %DEGREE_ELEVATION 
    %
    %  Function = DEGREE_ELEVATION(self, int order, index arg)
    %  Function = DEGREE_ELEVATION(self, [int] orders, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(340, self, varargin{:});
    end
    function varargout = kick_boundary(self,varargin)
    %KICK_BOUNDARY 
    %
    %  Function = KICK_BOUNDARY(self, TensorDomain boundary)
    %  Function = KICK_BOUNDARY(self, TensorDomain boundary, index arg_ind)
    %  Function = KICK_BOUNDARY(self, TensorDomain boundary, [char] args)
    %  Function = KICK_BOUNDARY(self, TensorDomain boundary, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(341, self, varargin{:});
    end
    function varargout = derivative(self,varargin)
    %DERIVATIVE 
    %
    %  Function = DERIVATIVE(self, int order, index arg)
    %  Function = DERIVATIVE(self, [int] orders, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(342, self, varargin{:});
    end
    function varargout = antiderivative(self,varargin)
    %ANTIDERIVATIVE 
    %
    %  Function = ANTIDERIVATIVE(self, int order, index arg)
    %  Function = ANTIDERIVATIVE(self, [int] orders, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(343, self, varargin{:});
    end
    function varargout = jacobian(self,varargin)
    %JACOBIAN 
    %
    %  [Function] = JACOBIAN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(344, self, varargin{:});
    end
    function varargout = hessian(self,varargin)
    %HESSIAN 
    %
    %  [[Function]] = HESSIAN(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(345, self, varargin{:});
    end
    function varargout = jacobian_matrix(self,varargin)
    %JACOBIAN_MATRIX 
    %
    %  Function = JACOBIAN_MATRIX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(346, self, varargin{:});
    end
    function varargout = hessian_matrix(self,varargin)
    %HESSIAN_MATRIX 
    %
    %  Function = HESSIAN_MATRIX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(347, self, varargin{:});
    end
    function varargout = integral(self,varargin)
    %INTEGRAL 
    %
    %  AnyTensor = INTEGRAL(self)
    %  AnyTensor = INTEGRAL(self, TensorDomain domain)
    %
    %
      [varargout{1:nargout}] = splinesMEX(348, self, varargin{:});
    end
    function varargout = partial_integral(self,varargin)
    %PARTIAL_INTEGRAL 
    %
    %  Function = PARTIAL_INTEGRAL(self, TensorDomain domain, [char] args)
    %  Function = PARTIAL_INTEGRAL(self, TensorDomain domain, [index] arg_ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(349, self, varargin{:});
    end
    function varargout = to_casadi(self,varargin)
    %TO_CASADI 
    %
    %  Function = TO_CASADI(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(351, self, varargin{:});
    end
    function varargout = fast_eval(self,varargin)
    %FAST_EVAL 
    %
    %  AnyTensor = FAST_EVAL(self, AnyTensor xy)
    %
    %
      [varargout{1:nargout}] = splinesMEX(352, self, varargin{:});
    end
    function varargout = fast_jac(self,varargin)
    %FAST_JAC 
    %
    %  DM = FAST_JAC(self, AnyTensor xy)
    %
    %
      [varargout{1:nargout}] = splinesMEX(353, self, varargin{:});
    end

   function varargout = subsref(self,s)
      if numel(s)==1 && strcmp(s(1).type,'()')
        [varargout{1:nargout}] = paren(self, s(1).subs);
      elseif numel(s)>1 && strcmp(s(1).type,'()')
        r = paren(self, s(1).subs);
        [varargout{1:nargout}] = subsref(r,s(2:end));
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
   end
    function ind = end(self,i,n)
      if n==1
        ind = size(self,1)*size(self,2);
      else
        ind = size(self, i);
      end
    end
    
    function self = subsasgn(self,varargin)
        error('Not supported: subsasgn');
    end
    function varargout = size(self, varargin)
      siz = shape(self);
      if nargin > 1
        if length(siz) < varargin{1}
          siz = 1;
        else
          siz = siz(varargin{1});
        end
      end
      varargout = cell(nargout);
      switch nargout
        case 0, disp(siz);
        case 1, varargout{1} = siz;
        case 2, varargout{1} = siz(1); varargout{2} = siz(2);
      end
    end
    function self = mrdivide(self,other)
      if ~isnumeric(other)
        [b,other] = isconst(other);
        assert(b,'Cannot divide by a spline');
      end
      self = mtimes(self,inv(other));
    end
    function self = mldivide(self,other)
      if ~isnumeric(self)
        [b,self] = isconst(self);
        assert(b,'Cannot divide by a spline');
      end
      self = mtimes(inv(self),other);
    end
  
   function out = index_helper(self, varargin)
      args = [];
      for i=1:numel(varargin)
        if strcmp(varargin{i},':')
          args = [args -1];
        else
          args = [args varargin{i}-1];
        end
      end
      out = index(self,args);
   end

  function out = paren(self,args)
    out = self.slice(args{:});
  end

  function out = eval(self,varargin)
      assert(numel(varargin)>0);
      if (numel(varargin)==1)
        out=self.eval_internal(varargin{1});
      else
        arg2 = varargin{2};
        arguments = false;
        if ischar(arg2)
          arguments = true;
        end
        if iscell(arg2) && numel(arg2)>0 && ischar(arg2{1})
          arguments = true;
        end

        if (arguments)
            out = self.eval_internal(varargin{:});
        else
            out = self.eval_internal([varargin{:}]);
        end
      end

  end

      function varargout = power(varargin)
    %POWER 
    %
    %  Function = POWER(Function lhs, int rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(354, varargin{:});
    end
    function varargout = mpower(varargin)
    %MPOWER 
    %
    %  Function = MPOWER(Function lhs, int rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(355, varargin{:});
    end
    function varargout = plus(varargin)
    %PLUS 
    %
    %  Function = PLUS(self, AnyTensor t)
    %  Function = PLUS(self, Function f)
    %  Function = PLUS(AnyTensor lhs, Function rhs)
    %  Function = PLUS(Function lhs, AnyTensor rhs)
    %  Function = PLUS(Function lhs, Function rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(356, varargin{:});
    end
    function varargout = minus(varargin)
    %MINUS 
    %
    %  Function = MINUS(self, AnyTensor t)
    %  Function = MINUS(self, Function f)
    %  Function = MINUS(AnyTensor lhs, Function rhs)
    %  Function = MINUS(Function lhs, AnyTensor rhs)
    %  Function = MINUS(Function lhs, Function rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(357, varargin{:});
    end
    function varargout = times(varargin)
    %TIMES 
    %
    %  Function = TIMES(self, AnyTensor t)
    %  Function = TIMES(self, Function f)
    %  Function = TIMES(AnyTensor lhs, Function rhs)
    %  Function = TIMES(Function lhs, AnyTensor rhs)
    %  Function = TIMES(Function lhs, Function rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(358, varargin{:});
    end
    function varargout = mtimes(varargin)
    %MTIMES 
    %
    %  Function = MTIMES(self, AnyTensor f)
    %  Function = MTIMES(self, Function f)
    %  Function = MTIMES(AnyTensor lhs, Function rhs)
    %  Function = MTIMES(Function lhs, AnyTensor rhs)
    %  Function = MTIMES(Function lhs, Function rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(359, varargin{:});
    end
    function varargout = rmtimes(varargin)
    %RMTIMES 
    %
    %  Function = RMTIMES(self, AnyTensor f)
    %  Function = RMTIMES(AnyTensor lhs, Function rhs)
    %  Function = RMTIMES(Function lhs, AnyTensor rhs)
    %  Function = RMTIMES(Function lhs, Function rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(360, varargin{:});
    end
    function varargout = le(varargin)
    %LE 
    %
    %  MX = LE(self, MX x)
    %  MX = LE(self, Function f)
    %  MX = LE(MX lhs, Function rhs)
    %  MX = LE(Function lhs, MX rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(361, varargin{:});
    end
    function varargout = ctranspose(self,varargin)
    %CTRANSPOSE 
    %
    %  Function = CTRANSPOSE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(362, self, varargin{:});
    end
    function self = Function(varargin)
    %FUNCTION 
    %
    %  new_obj = FUNCTION()
    %  new_obj = FUNCTION(AnyTensor tensor)
    %  new_obj = FUNCTION(AnyScalar value, [int] shape)
    %  new_obj = FUNCTION(TensorBasis basis, Coefficient coeff)
    %
    %
      self@splines.SharedObject(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(363, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(364, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
    function varargout = blockcat(varargin)
    %BLOCKCAT 
    %
    %  Function = BLOCKCAT([[Function]] f)
    %
    %
     [varargout{1:nargout}] = splinesMEX(321, varargin{:});
    end
    function varargout = cat(varargin)
    %CAT 
    %
    %  Function = CAT(index index, [Function] f)
    %
    %
     [varargout{1:nargout}] = splinesMEX(322, varargin{:});
    end
    function varargout = linear(varargin)
    %LINEAR 
    %
    %  Function = LINEAR(AnyVector x, AnyVector y)
    %
    %
     [varargout{1:nargout}] = splinesMEX(350, varargin{:});
    end
  end
end
