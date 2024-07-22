function varargout = spline_plus(varargin)
    %SPLINE_PLUS 
    %
    %  AnyTensor = SPLINE_PLUS(AnyTensor lhs, AnyTensor rhs)
    %  Function = SPLINE_PLUS(AnyTensor lhs, Function rhs)
    %  Function = SPLINE_PLUS(Function lhs, AnyTensor rhs)
    %  Function = SPLINE_PLUS(Function lhs, Function rhs)
    %
    %
  [varargout{1:nargout}] = splinesMEX(409, varargin{:});
end
