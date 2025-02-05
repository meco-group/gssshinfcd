function varargout = reverse(varargin)
    %REVERSE 
    %
    %  {IM}} = REVERSE({IM} ex, {IM} arg, {IM}} v, struct opts)
    %  {{DM}} = REVERSE({DM} ex, {DM} arg, {{DM}} v, struct opts)
    %  {SX}} = REVERSE({SX} ex, {SX} arg, {SX}} v, struct opts)
    %  {{MX}} = REVERSE({MX} ex, {MX} arg, {{MX}} v, struct opts)
    %
    %
  [varargout{1:nargout}] = casadiMEX(918, varargin{:});
end
