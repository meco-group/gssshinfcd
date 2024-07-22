function varargout = spline_mtimes(varargin)
    %SPLINE_MTIMES 
    %
    %  AnyTensor = SPLINE_MTIMES(AnyTensor lhs, AnyTensor rhs)
    %  Function = SPLINE_MTIMES(AnyTensor lhs, Function rhs)
    %  Function = SPLINE_MTIMES(Function lhs, AnyTensor rhs)
    %  Function = SPLINE_MTIMES(Function lhs, Function rhs)
    %
    %
  [varargout{1:nargout}] = splinesMEX(412, varargin{:});
end
