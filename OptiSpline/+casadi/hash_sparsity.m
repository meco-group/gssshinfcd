function varargout = hash_sparsity(varargin)
    %HASH_SPARSITY 
    %
    %  std::size_t = HASH_SPARSITY(int nrow, int ncol, int const * colind, int const * row)
    %  std::size_t = HASH_SPARSITY(int nrow, int ncol, [int] colind, [int] row)
    %    Hash a sparsity pattern.
    %
    %> HASH_SPARSITY(int nrow, int ncol, int const * colind, int const * row)
    %------------------------------------------------------------------------
    %
    %
    %
    %> HASH_SPARSITY(int nrow, int ncol, [int] colind, [int] row)
    %------------------------------------------------------------------------
    %
    %
    %Hash a sparsity pattern.
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(166, varargin{:});
end
