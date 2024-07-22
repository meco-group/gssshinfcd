classdef (InferiorClasses = {?casadi.DM,?casadi.SX,?casadi.MX}) DTensor < splines.PrintDT
    %DTENSOR 
    %
    %
    %
  methods
    function self = DTensor(varargin)
    %DTENSOR 
    %
    %  new_obj = DTENSOR(double a)
    %  new_obj = DTENSOR(DM data)
    %  new_obj = DTENSOR(DTensor t)
    %  new_obj = DTENSOR(DM data, [int] dims)
    %
    %
      self@splines.PrintDT(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(50, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(51, self);
        self.swigPtr=[];
      end
    end
    function varargout = is_scalar(self,varargin)
    %IS_SCALAR 
    %
    %  bool = IS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(52, self, varargin{:});
    end
    function varargout = as_scalar(self,varargin)
    %AS_SCALAR 
    %
    %  DM = AS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(53, self, varargin{:});
    end
    function varargout = is_vector(self,varargin)
    %IS_VECTOR 
    %
    %  bool = IS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(54, self, varargin{:});
    end
    function varargout = as_vector(self,varargin)
    %AS_VECTOR 
    %
    %  DTensor = AS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(55, self, varargin{:});
    end
    function varargout = data(self,varargin)
    %DATA 
    %
    %  DM = DATA(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(56, self, varargin{:});
    end
    function varargout = matrix(self,varargin)
    %MATRIX 
    %
    %  DM = MATRIX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(57, self, varargin{:});
    end
    function varargout = squeeze(self,varargin)
    %SQUEEZE 
    %
    %  DTensor = SQUEEZE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(58, self, varargin{:});
    end
    function varargout = squeeze_tailing(self,varargin)
    %SQUEEZE_TAILING 
    %
    %  DTensor = SQUEEZE_TAILING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(59, self, varargin{:});
    end
    function varargout = shape(self,varargin)
    %SHAPE 
    %
    %  DTensor = SHAPE(self, [int] dims)
    %
    %
      [varargout{1:nargout}] = splinesMEX(60, self, varargin{:});
    end
    function varargout = numel(self,varargin)
    %NUMEL 
    %
    %  int = NUMEL(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(61, self, varargin{:});
    end
    function varargout = diff(self,varargin)
    %DIFF 
    %
    %  DTensor = DIFF(self, int n, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(64, self, varargin{:});
    end
    function varargout = sum(self,varargin)
    %SUM 
    %
    %  DTensor = SUM(self, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(65, self, varargin{:});
    end
    function varargout = get_slice(self,varargin)
    %GET_SLICE 
    %
    %  DTensor = GET_SLICE(self, AnySlice i)
    %  DTensor = GET_SLICE(self, AnySlice i, AnySlice j)
    %
    %
      [varargout{1:nargout}] = splinesMEX(66, self, varargin{:});
    end
    function varargout = n_dims(self,varargin)
    %N_DIMS 
    %
    %  int = N_DIMS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(69, self, varargin{:});
    end
    function varargout = dims(self,varargin)
    %DIMS 
    %
    %  [int] = DIMS(self)
    %  int = DIMS(self, int i)
    %
    %
      [varargout{1:nargout}] = splinesMEX(70, self, varargin{:});
    end
    function varargout = empty(self,varargin)
    %EMPTY 
    %
    %  bool = EMPTY(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(71, self, varargin{:});
    end
    function varargout = solve(self,varargin)
    %SOLVE 
    %
    %  DTensor = SOLVE(self, DTensor B)
    %
    %
      [varargout{1:nargout}] = splinesMEX(73, self, varargin{:});
    end
    function varargout = uminus(self,varargin)
    %UMINUS 
    %
    %  DTensor = UMINUS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(74, self, varargin{:});
    end
    function varargout = le(self,varargin)
    %LE 
    %
    %  DTensor = LE(self, DTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(75, self, varargin{:});
    end
    function varargout = ge(self,varargin)
    %GE 
    %
    %  DTensor = GE(self, DTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(76, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  DTensor = EQ(self, DTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(77, self, varargin{:});
    end
    function varargout = paren(self,varargin)
    %PAREN 
    %
    %  DTensor = PAREN(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(78, self, varargin{:});
    end
    function varargout = transform(self,varargin)
    %TRANSFORM 
    %
    %  DTensor = TRANSFORM(self, DTensor tr, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(79, self, varargin{:});
    end
    function varargout = index(self,varargin)
    %INDEX 
    %
    %  DTensor = INDEX(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(80, self, varargin{:});
    end
    function varargout = index_anyslice(self,varargin)
    %INDEX_ANYSLICE 
    %
    %  DTensor = INDEX_ANYSLICE(self, [AnySlice] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(81, self, varargin{:});
    end
    function varargout = reorder_dims(self,varargin)
    %REORDER_DIMS 
    %
    %  DTensor = REORDER_DIMS(self, [index] order)
    %
    %
      [varargout{1:nargout}] = splinesMEX(82, self, varargin{:});
    end
    function varargout = einstein(self,varargin)
    %EINSTEIN 
    %
    %  DTensor = EINSTEIN(self, [int] a_e, [int] c_e)
    %  DTensor = EINSTEIN(self, DTensor B, [int] a, [int] b, [int] c)
    %
    %
      [varargout{1:nargout}] = splinesMEX(83, self, varargin{:});
    end
    function varargout = outer_product(self,varargin)
    %OUTER_PRODUCT 
    %
    %  DTensor = OUTER_PRODUCT(self, DTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(84, self, varargin{:});
    end
    function varargout = trailing_mtimes(self,varargin)
    %TRAILING_MTIMES 
    %
    %  DTensor = TRAILING_MTIMES(self, DTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(85, self, varargin{:});
    end
    function varargout = trailing_rmtimes(self,varargin)
    %TRAILING_RMTIMES 
    %
    %  DTensor = TRAILING_RMTIMES(self, DTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(86, self, varargin{:});
    end
    function varargout = inner(self,varargin)
    %INNER 
    %
    %  DTensor = INNER(self, DTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(87, self, varargin{:});
    end
    function varargout = partial_product(self,varargin)
    %PARTIAL_PRODUCT 
    %
    %  DTensor = PARTIAL_PRODUCT(self, DTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(88, self, varargin{:});
    end
    function varargout = to_string(self,varargin)
    %TO_STRING 
    %
    %  char = TO_STRING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(89, self, varargin{:});
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
    %  DTensor = PLUS(self, DTensor rhs)
    %  AnyTensor = PLUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(90, varargin{:});
    end
    function varargout = minus(varargin)
    %MINUS 
    %
    %  DTensor = MINUS(self, DTensor rhs)
    %  AnyTensor = MINUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(91, varargin{:});
    end
    function varargout = times(varargin)
    %TIMES 
    %
    %  AnyTensor = TIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(92, varargin{:});
    end
    function varargout = mtimes(varargin)
    %MTIMES 
    %
    %  DTensor = MTIMES(self, DTensor rhs)
    %  AnyTensor = MTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(93, varargin{:});
    end
    function varargout = rmtimes(varargin)
    %RMTIMES 
    %
    %  AnyTensor = RMTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(94, varargin{:});
    end
  end
  methods(Static)
    function varargout = concat(varargin)
    %CONCAT 
    %
    %  DTensor = CONCAT([DTensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(47, varargin{:});
    end
    function varargout = repeat(varargin)
    %REPEAT 
    %
    %  DTensor = REPEAT(DTensor e, [int] factors)
    %
    %
     [varargout{1:nargout}] = splinesMEX(48, varargin{:});
    end
    function varargout = pack(varargin)
    %PACK 
    %
    %  DTensor = PACK([DTensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(49, varargin{:});
    end
    function varargout = normalize_dim(varargin)
    %NORMALIZE_DIM 
    %
    %  [int,int] = NORMALIZE_DIM([int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(62, varargin{:});
    end
    function varargout = binary_dims(varargin)
    %BINARY_DIMS 
    %
    %  [int] = BINARY_DIMS(DTensor a, DTensor b)
    %
    %
     [varargout{1:nargout}] = splinesMEX(63, varargin{:});
    end
    function varargout = get(varargin)
    %GET 
    %
    %  DM = GET(DM data, [int] dims, [int] ind)
    %
    %
     [varargout{1:nargout}] = splinesMEX(67, varargin{:});
    end
    function varargout = set(varargin)
    %SET 
    %
    %  SET(DM data, [int] dims, [int] ind, DM rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(68, varargin{:});
    end
    function varargout = sym(varargin)
    %SYM 
    %
    %  DTensor = SYM(char name, [int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(72, varargin{:});
    end
  end
end
