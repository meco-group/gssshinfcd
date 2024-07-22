function varargout = spline_minus(varargin)
    %SPLINE_MINUS 
    %
    %  AnyTensor = SPLINE_MINUS(AnyTensor lhs, AnyTensor rhs)
    %  Function = SPLINE_MINUS(AnyTensor lhs, Function rhs)
    %  Function = SPLINE_MINUS(Function lhs, AnyTensor rhs)
    %  Function = SPLINE_MINUS(Function lhs, Function rhs)
    %
    %
  [varargout{1:nargout}] = splinesMEX(410, varargin{:});
end
