classdef (InferiorClasses = {?casadi.DM,?casadi.SX,?casadi.MX}) MTensor < splines.PrintMT
    %MTENSOR 
    %
    %
    %
  methods
    function self = MTensor(varargin)
    %MTENSOR 
    %
    %  new_obj = MTENSOR(double a)
    %  new_obj = MTENSOR(MTensor t)
    %  new_obj = MTENSOR(MX data)
    %  new_obj = MTENSOR(MX data, [int] dims)
    %
    %
      self@splines.PrintMT(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = splinesMEX(146, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        splinesMEX(147, self);
        self.swigPtr=[];
      end
    end
    function varargout = is_scalar(self,varargin)
    %IS_SCALAR 
    %
    %  bool = IS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(148, self, varargin{:});
    end
    function varargout = as_scalar(self,varargin)
    %AS_SCALAR 
    %
    %  MX = AS_SCALAR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(149, self, varargin{:});
    end
    function varargout = is_vector(self,varargin)
    %IS_VECTOR 
    %
    %  bool = IS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(150, self, varargin{:});
    end
    function varargout = as_vector(self,varargin)
    %AS_VECTOR 
    %
    %  MTensor = AS_VECTOR(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(151, self, varargin{:});
    end
    function varargout = data(self,varargin)
    %DATA 
    %
    %  MX = DATA(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(152, self, varargin{:});
    end
    function varargout = matrix(self,varargin)
    %MATRIX 
    %
    %  MX = MATRIX(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(153, self, varargin{:});
    end
    function varargout = squeeze(self,varargin)
    %SQUEEZE 
    %
    %  MTensor = SQUEEZE(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(154, self, varargin{:});
    end
    function varargout = squeeze_tailing(self,varargin)
    %SQUEEZE_TAILING 
    %
    %  MTensor = SQUEEZE_TAILING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(155, self, varargin{:});
    end
    function varargout = shape(self,varargin)
    %SHAPE 
    %
    %  MTensor = SHAPE(self, [int] dims)
    %
    %
      [varargout{1:nargout}] = splinesMEX(156, self, varargin{:});
    end
    function varargout = numel(self,varargin)
    %NUMEL 
    %
    %  int = NUMEL(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(157, self, varargin{:});
    end
    function varargout = diff(self,varargin)
    %DIFF 
    %
    %  MTensor = DIFF(self, int n, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(160, self, varargin{:});
    end
    function varargout = sum(self,varargin)
    %SUM 
    %
    %  MTensor = SUM(self, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(161, self, varargin{:});
    end
    function varargout = get_slice(self,varargin)
    %GET_SLICE 
    %
    %  MTensor = GET_SLICE(self, AnySlice i)
    %  MTensor = GET_SLICE(self, AnySlice i, AnySlice j)
    %
    %
      [varargout{1:nargout}] = splinesMEX(162, self, varargin{:});
    end
    function varargout = n_dims(self,varargin)
    %N_DIMS 
    %
    %  int = N_DIMS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(165, self, varargin{:});
    end
    function varargout = dims(self,varargin)
    %DIMS 
    %
    %  [int] = DIMS(self)
    %  int = DIMS(self, int i)
    %
    %
      [varargout{1:nargout}] = splinesMEX(166, self, varargin{:});
    end
    function varargout = empty(self,varargin)
    %EMPTY 
    %
    %  bool = EMPTY(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(167, self, varargin{:});
    end
    function varargout = solve(self,varargin)
    %SOLVE 
    %
    %  MTensor = SOLVE(self, MTensor B)
    %
    %
      [varargout{1:nargout}] = splinesMEX(169, self, varargin{:});
    end
    function varargout = uminus(self,varargin)
    %UMINUS 
    %
    %  MTensor = UMINUS(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(170, self, varargin{:});
    end
    function varargout = le(self,varargin)
    %LE 
    %
    %  MTensor = LE(self, MTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(171, self, varargin{:});
    end
    function varargout = ge(self,varargin)
    %GE 
    %
    %  MTensor = GE(self, MTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(172, self, varargin{:});
    end
    function varargout = eq(self,varargin)
    %EQ 
    %
    %  MTensor = EQ(self, MTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(173, self, varargin{:});
    end
    function varargout = paren(self,varargin)
    %PAREN 
    %
    %  MTensor = PAREN(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(174, self, varargin{:});
    end
    function varargout = transform(self,varargin)
    %TRANSFORM 
    %
    %  MTensor = TRANSFORM(self, MTensor tr, index axis)
    %
    %
      [varargout{1:nargout}] = splinesMEX(175, self, varargin{:});
    end
    function varargout = index(self,varargin)
    %INDEX 
    %
    %  MTensor = INDEX(self, [int] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(176, self, varargin{:});
    end
    function varargout = index_anyslice(self,varargin)
    %INDEX_ANYSLICE 
    %
    %  MTensor = INDEX_ANYSLICE(self, [AnySlice] ind)
    %
    %
      [varargout{1:nargout}] = splinesMEX(177, self, varargin{:});
    end
    function varargout = reorder_dims(self,varargin)
    %REORDER_DIMS 
    %
    %  MTensor = REORDER_DIMS(self, [index] order)
    %
    %
      [varargout{1:nargout}] = splinesMEX(178, self, varargin{:});
    end
    function varargout = einstein(self,varargin)
    %EINSTEIN 
    %
    %  MTensor = EINSTEIN(self, [int] a_e, [int] c_e)
    %  MTensor = EINSTEIN(self, MTensor B, [int] a, [int] b, [int] c)
    %
    %
      [varargout{1:nargout}] = splinesMEX(179, self, varargin{:});
    end
    function varargout = outer_product(self,varargin)
    %OUTER_PRODUCT 
    %
    %  MTensor = OUTER_PRODUCT(self, MTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(180, self, varargin{:});
    end
    function varargout = trailing_mtimes(self,varargin)
    %TRAILING_MTIMES 
    %
    %  MTensor = TRAILING_MTIMES(self, MTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(181, self, varargin{:});
    end
    function varargout = trailing_rmtimes(self,varargin)
    %TRAILING_RMTIMES 
    %
    %  MTensor = TRAILING_RMTIMES(self, MTensor rhs)
    %
    %
      [varargout{1:nargout}] = splinesMEX(182, self, varargin{:});
    end
    function varargout = inner(self,varargin)
    %INNER 
    %
    %  MTensor = INNER(self, MTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(183, self, varargin{:});
    end
    function varargout = partial_product(self,varargin)
    %PARTIAL_PRODUCT 
    %
    %  MTensor = PARTIAL_PRODUCT(self, MTensor b)
    %
    %
      [varargout{1:nargout}] = splinesMEX(184, self, varargin{:});
    end
    function varargout = to_string(self,varargin)
    %TO_STRING 
    %
    %  char = TO_STRING(self)
    %
    %
      [varargout{1:nargout}] = splinesMEX(185, self, varargin{:});
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
    %  MTensor = PLUS(self, MTensor rhs)
    %  AnyTensor = PLUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(186, varargin{:});
    end
    function varargout = minus(varargin)
    %MINUS 
    %
    %  MTensor = MINUS(self, MTensor rhs)
    %  AnyTensor = MINUS(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(187, varargin{:});
    end
    function varargout = times(varargin)
    %TIMES 
    %
    %  AnyTensor = TIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(188, varargin{:});
    end
    function varargout = mtimes(varargin)
    %MTIMES 
    %
    %  MTensor = MTIMES(self, MTensor rhs)
    %  AnyTensor = MTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(189, varargin{:});
    end
    function varargout = rmtimes(varargin)
    %RMTIMES 
    %
    %  AnyTensor = RMTIMES(AnyTensor lhs, AnyTensor rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(190, varargin{:});
    end
  end
  methods(Static)
    function varargout = concat(varargin)
    %CONCAT 
    %
    %  MTensor = CONCAT([MTensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(143, varargin{:});
    end
    function varargout = repeat(varargin)
    %REPEAT 
    %
    %  MTensor = REPEAT(MTensor e, [int] factors)
    %
    %
     [varargout{1:nargout}] = splinesMEX(144, varargin{:});
    end
    function varargout = pack(varargin)
    %PACK 
    %
    %  MTensor = PACK([MTensor] v, index axis)
    %
    %
     [varargout{1:nargout}] = splinesMEX(145, varargin{:});
    end
    function varargout = normalize_dim(varargin)
    %NORMALIZE_DIM 
    %
    %  [int,int] = NORMALIZE_DIM([int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(158, varargin{:});
    end
    function varargout = binary_dims(varargin)
    %BINARY_DIMS 
    %
    %  [int] = BINARY_DIMS(MTensor a, MTensor b)
    %
    %
     [varargout{1:nargout}] = splinesMEX(159, varargin{:});
    end
    function varargout = get(varargin)
    %GET 
    %
    %  MX = GET(MX data, [int] dims, [int] ind)
    %
    %
     [varargout{1:nargout}] = splinesMEX(163, varargin{:});
    end
    function varargout = set(varargin)
    %SET 
    %
    %  SET(MX data, [int] dims, [int] ind, MX rhs)
    %
    %
     [varargout{1:nargout}] = splinesMEX(164, varargin{:});
    end
    function varargout = sym(varargin)
    %SYM 
    %
    %  MTensor = SYM(char name, [int] dims)
    %
    %
     [varargout{1:nargout}] = splinesMEX(168, varargin{:});
    end
  end
end
