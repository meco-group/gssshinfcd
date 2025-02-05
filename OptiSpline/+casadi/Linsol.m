classdef  Linsol < casadi.SharedObject & casadi.PrintableCommon
    %LINSOL Linear solver Create a solver for linear systems of equations Solves the
    %
    %
    %linear system A*X = B or A^T*X = B for X with A square and non- singular.
    %
    %If A is structurally singular, an error will be thrown during init. If A is
    %numerically singular, the prepare step will fail.
    %
    %General information
    %===================
    %
    %
    %
    %List of plugins
    %===============
    %
    %
    %
    %- csparsecholesky
    %
    %- csparse
    %
    %- ma27
    %
    %- lapacklu
    %
    %- lapackqr
    %
    %- ldl
    %
    %- qr
    %
    %- symbolicqr
    %
    %Note: some of the plugins in this list might not be available on your
    %system. Also, there might be extra plugins available to you that are not
    %listed here. You can obtain their documentation with
    %Linsol.doc("myextraplugin")
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %csparsecholesky
    %---------------
    %
    %
    %
    %Linsol with CSparseCholesky Interface
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %csparse
    %-------
    %
    %
    %
    %Linsol with CSparse Interface
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %ma27
    %----
    %
    %
    %
    %Interface to the sparse direct linear solver MA27 Works for symmetric
    %indefinite systems Partly adopted from qpOASES 3.2 Joel Andersson
    %
    %--------------------------------------------------------------------------------
    %
    %lapacklu
    %--------
    %
    %
    %
    %This class solves the linear system A.x=b by making an LU factorization of
    %A: A = L.U, with L lower and U upper triangular
    %
    %>List of available options
    %
    %+-----------------------------+---------+----------------------------------+
    %|             Id              |  Type   |           Description            |
    %+=============================+=========+==================================+
    %| allow_equilibration_failure | OT_BOOL | Non-fatal error when             |
    %|                             |         | equilibration fails              |
    %+-----------------------------+---------+----------------------------------+
    %| equilibration               | OT_BOOL | Equilibrate the matrix           |
    %+-----------------------------+---------+----------------------------------+
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %lapackqr
    %--------
    %
    %
    %
    %This class solves the linear system A.x=b by making an QR factorization of
    %A: A = Q.R, with Q orthogonal and R upper triangular
    %
    %>List of available options
    %
    %+----------+--------+------------------------------------------------------+
    %|    Id    |  Type  |                     Description                      |
    %+==========+========+======================================================+
    %| max_nrhs | OT_INT | Maximum number of right-hand-sides that get          |
    %|          |        | processed in a single pass [default:10].             |
    %+----------+--------+------------------------------------------------------+
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %ldl
    %---
    %
    %
    %
    %Linear solver using sparse direct LDL factorization
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %qr --
    %
    %
    %
    %Linear solver using sparse direct QR factorization
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %
    %
    %--------------------------------------------------------------------------------
    %
    %symbolicqr
    %----------
    %
    %
    %
    %Linear solver for sparse least-squares problems Inspired
    %fromhttps://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/lsqr.py#L96
    %
    %Linsol based on QR factorization with sparsity pattern based reordering
    %without partial pivoting
    %
    %>List of available options
    %
    %+-------+---------+----------------------------------------------------+
    %|  Id   |  Type   |                    Description                     |
    %+=======+=========+====================================================+
    %| fopts | OT_DICT | Options to be passed to generated function objects |
    %+-------+---------+----------------------------------------------------+
    %
    %--------------------------------------------------------------------------------
    %
    %
    %
    %Joel Andersson
    %
    %C++ includes: linsol.hpp 
    %
  methods
    function this = swig_this(self)
      this = casadiMEX(3, self);
    end
    function varargout = plugin_name(self,varargin)
    %PLUGIN_NAME Query plugin name.
    %
    %  char = PLUGIN_NAME(self)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(881, self, varargin{:});
    end
    function varargout = sparsity(self,varargin)
    %SPARSITY Get linear system sparsity.
    %
    %  Sparsity = SPARSITY(self)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(882, self, varargin{:});
    end
    function varargout = sfact(self,varargin)
    %SFACT Symbolic factorization of the linear system, e.g. selecting pivots.
    %
    %  SFACT(self, DM A)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(883, self, varargin{:});
    end
    function varargout = nfact(self,varargin)
    %NFACT Numeric factorization of the linear system.
    %
    %  NFACT(self, DM A)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(884, self, varargin{:});
    end
    function varargout = solve(self,varargin)
    %SOLVE Solve linear system of equations
    %
    %  DM = SOLVE(self, DM A, DM B, bool tr)
    %  MX = SOLVE(self, MX A, MX B, bool tr)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(885, self, varargin{:});
    end
    function varargout = neig(self,varargin)
    %NEIG Number of negative eigenvalues Not available for all solvers.
    %
    %  int = NEIG(self, DM A)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(886, self, varargin{:});
    end
    function varargout = rank(self,varargin)
    %RANK Matrix rank Not available for all solvers.
    %
    %  int = RANK(self, DM A)
    %
    %
    %
    %
      [varargout{1:nargout}] = casadiMEX(887, self, varargin{:});
    end
    function self = Linsol(varargin)
    %LINSOL 
    %
    %  new_obj = LINSOL()
    %    Default constructor.
    %  new_obj = LINSOL(char name, char solver, Sparsity sp, struct opts)
    %    Constructor.
    %
    %> LINSOL()
    %------------------------------------------------------------------------
    %
    %
    %Default constructor.
    %
    %
    %> LINSOL(char name, char solver, Sparsity sp, struct opts)
    %------------------------------------------------------------------------
    %
    %
    %Constructor.
    %
    %
    %
      self@casadi.SharedObject(SwigRef.Null);
      self@casadi.PrintableCommon(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = casadiMEX(888, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        casadiMEX(889, self);
        self.swigPtr=[];
      end
    end
  end
  methods(Static)
    function varargout = type_name(varargin)
    %TYPE_NAME 
    %
    %  char = TYPE_NAME()
    %
    %
     [varargout{1:nargout}] = casadiMEX(876, varargin{:});
    end
    function varargout = test_cast(varargin)
    %TEST_CAST 
    %
    %  bool = TEST_CAST(casadi::SharedObjectInternal const * ptr)
    %
    %
     [varargout{1:nargout}] = casadiMEX(877, varargin{:});
    end
    function varargout = has_plugin(varargin)
    %HAS_PLUGIN 
    %
    %  bool = HAS_PLUGIN(char name)
    %
    %
     [varargout{1:nargout}] = casadiMEX(878, varargin{:});
    end
    function varargout = load_plugin(varargin)
    %LOAD_PLUGIN 
    %
    %  LOAD_PLUGIN(char name)
    %
    %
     [varargout{1:nargout}] = casadiMEX(879, varargin{:});
    end
    function varargout = doc(varargin)
    %DOC 
    %
    %  char = DOC(char name)
    %
    %
     [varargout{1:nargout}] = casadiMEX(880, varargin{:});
    end
  end
end
