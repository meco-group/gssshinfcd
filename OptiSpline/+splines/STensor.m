classdef (InferiorClasses = {?casadi.DM,?casadi.SX,?casadi.MX}) STensor < splines.PrintST
    %STENSOR 
    %
    %
    %
  methods
    function self = STensor(varargin)
    %STENSOR 
    %
    %  new_obj = STENSOR(double a)
    %  new_obj = STENSOR(SX data)
    %  new_obj = STENSOR(STensor t)
    %  new_obj = STENSOR(SX data, [int] dims)
    %
    %
      self@splines.PrintST(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(98, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(99, self);
        self.swigPtr=[];
      end
    end
    function varargout = is_scalar(self,varargin)
    %IS_SCALAR 
    %
    %  bool = IS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(100, self, varargin{:});
    end
    function varargout = as_scalar(self,varargin)
    %AS_SCALAR 
    %
    %  SX = AS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(101, self, varargin{:});
    end
    function varargout = is_vector(self,varargin)
    %IS_VECTOR 
    %
    %  bool = IS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(102, self, varargin{:});
    end
    function varargout = as_vector(self,varargin)
    %AS_VECTOR 
    %
    %  STensor = AS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(103, self, varargin{:});
    end
    function varargout = data(self,varargin)
    %DATA 
    %
    %  SX = DATA(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(104, self, varargin{:});
    end
    function varargout = matrix(self,varargin)
    %MATRIX 
    %
    %  SX = MATRIX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(105, self, varargin{:});
    end
    function varargout = squeeze(self,varargin)
    %SQUEEZE 
    %
    %  STensor = SQUEEZE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(106, self, varargin{:});
    end
    function varargout = squeeze_tailing(self,varargin)
    %SQUEEZE_TAILING 
    %
    %  STensor = SQUEEZE_TAILING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(107, self, varargin{:});
    end
    function varargout = shape(self,varargin)
    %SHAPE 
    %
    %  STensor = SHAPE(self, [int] dims)
    %
    %
      [varargout{1:nargout}] = splinesMEX(108, self, varargin{:});
    end
    function varargout = numel(self,varargin)
    %NUMEL 
    %
    %  int = NUMEL(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(109, self, varargin{:});
    end
    function varargout = diff(self,varargin)
    %DIFF 
    %
    %  STensor = DIFF(self, int n, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(112, self, varargin{:});
    end
    function varargout = sum(self,varargin)
    %SUM 
    %
    %  STensor = SUM(self, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(113, self, varargin{:});
    end
    function varargout = get_slice(self,varargin)
    %GET_SLICE 
    %
    %  STensor = GET_SLICE(self, AnySlice i)
    %  STensor = GET_SLICE(self, AnySlice i, AnySlice j)
    %
    %
      [varargout{1:nargout}] = splinesMEX(114, self, varargin{:});
    end
    function varargout = n_dims(self,varargin)
    %N_DIMS 
    %
    %  int = N_DIMS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(117, self, varargin{:});
    end
    function varargout = dims(self,varargin)
    %DIMS 
    %
    %  [int] = DIMS(self)
    %  int = DIMS(self, int i)
    %
    %
      [varargout{1:nargout}] = splinesMEX(118, self, varargin{:});
    end
    function varargout = empty(self,varargin)
    %EMPTY 
    %
    %  bool = EMPTY(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(119, self, varargin{:});
    end
    function varargout = solve(self,varargin)
    %SOLVE 
    %
    %  STensor = SOLVE(self, STensor B)
    %
    %
      [varargout{1:nargout}] = splinesMEX(121, self, varargin{:});
    end
    function varargout = uminus(self,varargin)
    %UMINUS 
    %
    %  STensor = UMINUS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(122, self, varargin{:});
    end
    function varargout = le(self,varargin)
    %LE 
    %
    %  STensor = LE(self, STensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(123, self, varargin{:});
    end
    function varargout = ge(self,varargin)
    %GE 
    %
    %  STensor = GE(self, STensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(124, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  STensor = EQ(self, STensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(125, self, varargin{:});
    end
    function varargout = paren(self,varargin)
    %PAREN 
    %
    %  STensor = PAREN(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(126, self, varargin{:});
    end
    function varargout = transform(self,varargin)
    %TRANSFORM 
    %
    %  STensor = TRANSFORM(self, STensor tr, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(127, self, varargin{:});
    end
    function varargout = index(self,varargin)
    %INDEX 
    %
    %  STensor = INDEX(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(128, self, varargin{:});
    end
    function varargout = index_anyslice(self,varargin)
    %INDEX_ANYSLICE 
    %
    %  STensor = INDEX_ANYSLICE(self, [AnySlice] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(129, self, varargin{:});
    end
    function varargout = reorder_dims(self,varargin)
    %REORDER_DIMS 
    %
    %  STensor = REORDER_DIMS(self, [index] order)
    %
    %
      [varargout{1:nargout}] = splinesMEX(130, self, varargin{:});
    end
    function varargout = einstein(self,varargin)
    %EINSTEIN 
    %
    %  STensor = EINSTEIN(self, [int] a_e, [int] c_e)
    %  STensor = EINSTEIN(self, STensor B, [int] a, [int] b, [int] c)
    %
    %
      [varargout{1:nargout}] = splinesMEX(131, self, varargin{:});
    end
    function varargout = outer_product(self,varargin)
    %OUTER_PRODUCT 
    %
    %  STensor = OUTER_PRODUCT(self, STensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(132, self, varargin{:});
    end
    function varargout = trailing_mtimes(self,varargin)
    %TRAILING_MTIMES 
    %
    %  STensor = TRAILING_MTIMES(self, STensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(133, self, varargin{:});
    end
    function varargout = trailing_rmtimes(self,varargin)
    %TRAILING_RMTIMES 
    %
    %  STensor = TRAILING_RMTIMES(self, STensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(134, self, varargin{:});
    end
    function varargout = inner(self,varargin)
    %INNER 
    %
    %  STensor = INNER(self, STensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(135, self, varargin{:});
    end
    function varargout = partial_product(self,varargin)
    %PARTIAL_PRODUCT 
    %
    %  STensor = PARTIAL_PRODUCT(self, STensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(136, self, varargin{:});
    end
    function varargout = to_string(self,varargin)
    %TO_STRING 
    %
    %  char = TO_STRING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(137, self, varargin{:});
    end

    function [] = disp(self)
      disp(self.to_string())
    end
  
   function varargout = subsref(self,s)
      if numel(s)==1 & s.type=='()'
        [varargout{1:nargout}]= index_helper(self, s.subs{:});
      else
        [varargout{1:nargout}] = builtin('subsref',self,s);
      end
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
    function varargout = plus(varargin)
    %PLUS 
    %
    %  STensor = PLUS(self, STensor rhs)
    %  AnyTensor = PLUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(138, varargin{:});
    end
    function varargout = minus(varargin)
    %MINUS 
    %
    %  STensor = MINUS(self, STensor rhs)
    %  AnyTensor = MINUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(139, varargin{:});
    end
    function varargout = times(varargin)
    %TIMES 
    %
    %  AnyTensor = TIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(140, varargin{:});
    end
    function varargout = mtimes(varargin)
    %MTIMES 
    %
    %  STensor = MTIMES(self, STensor rhs)
    %  AnyTensor = MTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(141, varargin{:});
    end
    function varargout = rmtimes(varargin)
    %RMTIMES 
    %
    %  AnyTensor = RMTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(142, varargin{:});
    end
  end
  methods(Static)
    function varargout = concat(varargin)
    %CONCAT 
    %
    %  STensor = CONCAT([STensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(95, varargin{:});
    end
    function varargout = repeat(varargin)
    %REPEAT 
    %
    %  STensor = REPEAT(STensor e, [int] factors)
    %
    %
     [varargout{1:nargout}] = splinesMEX(96, varargin{:});
    end
    function varargout = pack(varargin)
    %PACK 
    %
    %  STensor = PACK([STensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(97, varargin{:});
    end
    function varargout = normalize_dim(varargin)
    %NORMALIZE_DIM 
    %
    %  [int,int] = NORMALIZE_DIM([int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(110, varargin{:});
    end
    function varargout = binary_dims(varargin)
    %BINARY_DIMS 
    %
    %  [int] = BINARY_DIMS(STensor a, STensor b)
    %
    %
     [varargout{1:nargout}] = splinesMEX(111, varargin{:});
    end
    function varargout = get(varargin)
    %GET 
    %
    %  SX = GET(SX data, [int] dims, [int] ind)
    %
    %
     [varargout{1:nargout}] = splinesMEX(115, varargin{:});
    end
    function varargout = set(varargin)
    %SET 
    %
    %  SET(SX data, [int] dims, [int] ind, SX rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(116, varargin{:});
    end
    function varargout = sym(varargin)
    %SYM 
    %
    %  STensor = SYM(char name, [int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(120, varargin{:});
    end
  end
end
