# gssshinfcd - Gain-scheduled state-space $\mathcal{H}_\infty$ control design for linear parameter-varying systems

### Gain-scheduled version of [`sshinfcd`](https://github.com/meco-group/sshinfcd). 

Experimental implementation for internal usage. The methodology is based on

* P. Apkarian and R. J. Adams. ['Advanced gain-scheduling techniques for uncertain systems.'](https://doi.org/10.1109/87.654874) IEEE Transactions on Control Systems Technology, vol. 6, no. 1, 1998. 
* W. Van Loock et al. ['Approximate parametric cone programming with applications in control.'](https://doi.org/10.1109/ECC.2016.7810283) 2016 European Control Conference (ECC), pp. 178-183, 2016.
* G. Hilhorst et al. ['Control of Linear Parameter-Varying Systems using B-Splines.'](https://doi.org/10.1109/CDC.2016.7798757) 2016 IEEE 55th Conference on Decision and Control (CDC), pp. 3246-3251, 2016.

and relies on the [OptiSpline](https://github.com/meco-group/optispline) toolbox for efficiently handling tensor-product B-splines. Binaries for Windows are included; they can be replaced with binaries for Linux or macOS if needed. `gssshinfcd` neatly interfaces with [LCToolbox](https://github.com/meco-group/lc_toolbox). 
