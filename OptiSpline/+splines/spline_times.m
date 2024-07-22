function varargout = spline_times(varargin)
    %SPLINE_TIMES 
    %
    %  AnyTensor = SPLINE_TIMES(AnyTensor lhs, AnyTensor rhs)
    %  Function = SPLINE_TIMES(AnyTensor lhs, Function rhs)
    %  Function = SPLINE_TIMES(Function lhs, AnyTensor rhs)
    %  Function = SPLINE_TIMES(Function lhs, Function rhs)
    %
    %
  [varargout{1:nargout}] = splinesMEX(411, varargin{:});
end
