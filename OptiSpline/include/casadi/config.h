/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
#ifndef CASADI_CONFIG_H // NOLINT(build/header_guard)
#define CASADI_CONFIG_H // NOLINT(build/header_guard)

#define CASADI_MAJOR_VERSION 3
#define CASADI_MINOR_VERSION 3
#define CASADI_PATCH_VERSION 0
#define CASADI_IS_RELEASE 1
#define CASADI_VERSION_STRING "3.3.0"
#define CASADI_GIT_REVISION "fd942483253551563a308809698fcfcdd916678e"  // NOLINT(whitespace/line_length)
#define CASADI_GIT_DESCRIBE "3.2.0+656.fd9424832"  // NOLINT(whitespace/line_length)
#define CASADI_FEATURE_LIST ""  // NOLINT(whitespace/line_length)
#define CASADI_BUILD_TYPE "Release"  // NOLINT(whitespace/line_length)
#define CASADI_COMPILER_ID "GNU"  // NOLINT(whitespace/line_length)
#define CASADI_COMPILER "/usr/bin/x86_64-w64-mingw32-g++"  // NOLINT(whitespace/line_length)
#define CASADI_COMPILER_FLAGS " -fno-ipa-cp-clone -std=c++11 -fno-ipa-cp-clone   "  // NOLINT(whitespace/line_length)
#define CASADI_MODULES "casadi;casadi_linsol_lapacklu;casadi_linsol_lapackqr;casadi_sundials_common;casadi_integrator_cvodes;casadi_integrator_idas;casadi_rootfinder_kinsol;casadi_nlpsol_ipopt;casadi_nlpsol_bonmin;casadi_conic_qpoases;casadi_conic_clp;casadi_linsol_csparse;casadi_linsol_csparsecholesky;casadi_importer_clang;casadi_linsol_ma27;casadi_xmlfile_tinyxml;casadi_nlpsol_blocksqp;casadi_conic_hpmpc;casadi_dple_slicot;casadi_expm_slicot;casadi_nlpsol_ampl;casadi_conic_nlpsol;casadi_importer_shell;casadi_integrator_rk;casadi_integrator_collocation;casadi_interpolant_linear;casadi_interpolant_bspline;casadi_linsol_symbolicqr;casadi_linsol_qr;casadi_linsol_ldl;casadi_linsol_lsqr;casadi_nlpsol_sqpmethod;casadi_nlpsol_scpgen;casadi_rootfinder_newton;casadi_rootfinder_nlpsol"  // NOLINT(whitespace/line_length)
#define CASADI_PLUGINS "Linsol::lapacklu;Linsol::lapackqr;Integrator::cvodes;Integrator::idas;Rootfinder::kinsol;Nlpsol::ipopt;Nlpsol::bonmin;Conic::qpoases;Conic::clp;Linsol::csparse;Linsol::csparsecholesky;Importer::clang;Linsol::ma27;XmlFile::tinyxml;Nlpsol::blocksqp;Conic::hpmpc;Dple::slicot;Expm::slicot;Nlpsol::ampl;Conic::nlpsol;Importer::shell;Integrator::rk;Integrator::collocation;Interpolant::linear;Interpolant::bspline;Linsol::symbolicqr;Linsol::qr;Linsol::ldl;Linsol::lsqr;Nlpsol::sqpmethod;Nlpsol::scpgen;Rootfinder::newton;Rootfinder::nlpsol"  // NOLINT(whitespace/line_length)
#define CASADI_INSTALL_PREFIX "/home/travis/build/casadi/binaries/casadi/matlab_install"  // NOLINT(whitespace/line_length)

#endif  // CASADI_CONFIG_H // NOLINT(build/header_guard)
