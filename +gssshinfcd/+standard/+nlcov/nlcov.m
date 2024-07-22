function K = nlcov(prob,opts)
% NLCOV Solves the LMI formulation of the standard problem in state-space form
%
% [K,gamma] = NLCOV(prob, opts) returns an LFT representation of the
% parameter-varying controller that solves the SDP formulation of the 
% H-infinity control problem. 

% This file is part of gssshinfcd.
% Copyright (c) 2024, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% gssshinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% gssshinfcd is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
% License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with gssshinfcd. If not, see <https://www.gnu.org/licenses/>.

    assert(prob.Ts==0,'Discrete-time problems are not supported yet.');
    opti = OptiSplineYalmip();
    
    % 1. Decide on the parametrization of the Lyapunov matrices and the objective function
    
        % helper function: remove a basis from a tensor basis
        function ntb = remove(otb,arg)
            args = cellstr(arguments(otb));
            bass = bases(otb); 
            bass(strcmp(arg,args)) = [];
            args(strcmp(arg,args)) = []; 
            if isempty(bass) || isempty(args)
                ntb = splines.TensorBasis({},{});
            else
                ntb = splines.TensorBasis(bass,args); 
            end
        end

        % helper function: create a linear 1D tensor basis
        function ntb = createlinear(arg,range) 
            ntb = splines.TensorBasis({splines.BSplineBasis(range,1,2)},{arg}); 
        end
        
        % Lyapunov matrices
        solvedual = false; 
        if isempty(opts.synthesis.lyapunovbasis) 
            X = opti.variable(prob.n(),prob.n(),'symmetric');
            Y = opti.variable(prob.n(),prob.n(),'symmetric');
            Xs = opti.variable(prob.n(),prob.n(),'symmetric');
            Ys = opti.variable(prob.n(),prob.n(),'symmetric');
            Xsdot = zeros(prob.n(),prob.n());
            Ysdot = zeros(prob.n(),prob.n()); 
        else
            tb = opts.synthesis.lyapunovbasis; 
            
            % edit the tensor basis
            tbc = splines.TensorBasis({},{}); 
            for k=1:length(prob.P.parameters)
                
                if any(isinf(prob.P.parameters(k).rate)) % unbounded => replace parametrization by a constant Lyapunov matrix 
                    tb = remove(tb,prob.P.parameters(k).name); 
                end
                
                if ~any(prob.P.parameters(k).rate) % constant => "constant" Lyapunov matrix can also be a function of this parameter
                    b = bases(tb);
                    b = b(strcmp(prob.P.parameters(k).name,cellstr(arguments(tb)))); 
                    tbc = splines.TensorBasis([bases(tbc),b],[arguments(tbc), {prob.P.parameters(k).name}]); 
                end
                
            end
            
            % generate the Lyapunov matrices
            X = opti.Function(tbc,[prob.n(),prob.n()],'symmetric');
            Y = opti.Function(tbc,[prob.n(),prob.n()],'symmetric');
            Xs = opti.Function(tb,[prob.n(),prob.n()],'symmetric');
            Ys = opti.Function(tb,[prob.n(),prob.n()],'symmetric');
            
            % generate the derivatives
            Xsdot = zeros(prob.n(),prob.n());
            Ysdot = zeros(prob.n(),prob.n());
            for k=1:length(prob.P.parameters)
                if any(prob.P.parameters(k).rate) && ~any(isinf(prob.P.parameters(k).rate)) && ismember(prob.P.parameters(k).name,cellstr(arguments(tb)))
                    solvedual = true; 
                    Xsdot = Xsdot + Xs.derivative(1,{prob.P.parameters(k).name}) * opti.Function(createlinear(strcat('der_',prob.P.parameters(k).name),prob.P.parameters(k).rate));
                    Ysdot = Ysdot + Ys.derivative(1,{prob.P.parameters(k).name}) * opti.Function(createlinear(strcat('der_',prob.P.parameters(k).name),prob.P.parameters(k).rate));
                end
            end
        end     
        
        % controller variables -> same parameterization as plant
        if isempty(opts.synthesis.controllerbasis)
            P = [prob.P.A, prob.P.B ; prob.P.C, prob.P.D];
            TB = P.tensor_basis;
        else
            TB = opts.synthesis.controllerbasis; 
        end
        Ak = opti.Function(TB,[prob.n(),prob.n()],'full');
        Bk = opti.Function(TB,[prob.n(),prob.ny()],'full'); 
        Ck = opti.Function(TB,[prob.nu(),prob.n()],'full'); 
        Dk = opti.Function(TB,[prob.nu(),prob.ny()],'full');
        
        % objective
        isobj = ([prob.specs.weight]~=0);
        if isempty(opts.synthesis.gammabasis)
            sqgamma = opti.variable(sum(isobj),1);
            obj = vertcat(prob.specs(isobj).weight)'*sqgamma;
        else
            sqgamma = opti.Function(opts.synthesis.gammabasis,[sum(isobj),1]);
            obj = vertcat(prob.specs(isobj).weight)'*sqgamma;
            obj = obj.integral(); % L1 norm by default
        end
   
        decvars = struct('X',X, ...
                         'Y',Y, ...
                         'Xs',Xs, ...
                         'Ys',Ys, ...
                         'Xsdot',Xsdot, ...
                         'Ysdot',Ysdot, ...
                         'Ak',Ak, ...
                         'Bk',Bk, ...
                         'Ck',Ck, ...
                         'Dk',Dk, ...
                         'sqgamma',sqgamma, ...
                         'obj', obj); 
        
    % 2. Build the LMIs

        % helper functions
        function h = He(x)
            h = x+x'; 
        end
    
        function e = safe_degree_elevation(o,d)
            if isa(o,'splines.Function')
                e = o.degree_elevation(d);
            else
                e = o;
            end
        end
    
        function e = safe_midpoint_refinement(o,n)
            if isa(o,'splines.Function')
                e = o.degree_elevation(n);
            else
                e = o;
            end
        end
    
            
        % 2.1 The "primal" version (X parameter-varying, Y "constant")

        % build the set of LMIs
        LMIs = cell(length(prob.specs),1); 
        
            % H-infinity performance
            j = 0; 
            for k=1:length(prob.specs)
                Bw = prob.Bw(); Bw = Bw(:,prob.specs(k).in); nwk = length(prob.specs(k).in);
                Cz = prob.Cz(); Cz = Cz(prob.specs(k).out,:); nzk = length(prob.specs(k).out); 
                Dzw = prob.Dzw(); Dzw = Dzw(prob.specs(k).out,prob.specs(k).in);
                Dzu = prob.Dzu(); Dzu = Dzu(prob.specs(k).out,:);
                Dyw = prob.Dyw(); Dyw = Dyw(:,prob.specs(k).in); 
                LMIs{k} = He([Xs*prob.A()+Bk*prob.Cy()                  zeros(prob.n())               zeros(prob.n(),nwk)         zeros(prob.n(),nzk);
                              Ak'+prob.A()+prob.Bu()*Dk*prob.Cy()       prob.A()*Y+prob.Bu()*Ck       zeros(prob.n(),nwk)         zeros(prob.n(),nzk());  
                              (Xs*Bw+Bk*Dyw)'                           (Bw+prob.Bu()*Dk*Dyw)'        zeros(nwk,nwk)              zeros(nwk,nzk); 
                              Cz+Dzu*Dk*prob.Cy()                       Cz*Y+Dzu*Ck                   Dzw+Dzu*Dk*Dyw              zeros(nzk,nzk)]);

                if isobj(k)   
                    j = j+1;
                    switch opts.synthesis.lyapunovshape 
                        case 1
                          LMIs{k} = LMIs{k} + [Xsdot                       zeros(prob.n())                zeros(prob.n(),nwk)       zeros(prob.n(),nzk);
                                               zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)       zeros(prob.n(),nzk); 
                                               zeros(nwk,prob.n())         zeros(nwk,prob.n())            -sqgamma(j)*eye(nwk)      zeros(nwk,nzk); 
                                               zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)            -sqgamma(j)*eye(nzk)];
                        case 2
                          LMIs{k} = LMIs{k} + [Xsdot                       zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                               zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                               zeros(nwk,prob.n())         zeros(nwk,prob.n())            -sqgamma(j)*eye(nwk)     zeros(nwk,nzk); 
                                               zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -eye(nzk)];
                        case 3
                          LMIs{k} = LMIs{k} + [Xsdot                       zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                               zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                               zeros(nwk,prob.n())         zeros(nwk,prob.n())            -eye(nwk)                zeros(nwk,nzk); 
                                               zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)          -sqgamma(j)*eye(nzk)];
                        otherwise
                            prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
                    end
                else
                    LMIs{k} = LMIs{k} + [Xsdot                       zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                         zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                         zeros(nwk,prob.n())         zeros(nwk,prob.n())            -eye(nwk)                zeros(nwk,nzk); 
                                         zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -eye(nzk)];
                end
            end
            
            % stability constraint
            LMIs = [LMIs ; {-[Xs eye(prob.n()) ; eye(prob.n()) Y]}];
            
            % (D-)stability constraints, if any
            for k=1:length(prob.region)
                if ~isempty(prob.region(k).M) && ~isempty(prob.region(k).L)
                    LMIs = [LMIs ; 
                            {kron(prob.region(k).L,[Y             eye(prob.n()); 
                                                    eye(prob.n()) Xs           ]) 
                             + He(kron(prob.region(k).M,[prob.A()*Y+prob.Bu()*Ck   prob.A()+prob.Bu()*Dk*prob.Cy(); 
                                                         Ak                        Xs*prob.A()+Bk*prob.Cy()       ]))}];
                end
            end
            
        LMIs = cellfun(@(x) safe_degree_elevation(x,opts.synthesis.degreeelevations), LMIs, 'un', 0);
        LMIs = cellfun(@(x) safe_midpoint_refinement(x,opts.synthesis.knotinsertions), LMIs, 'un', 0);
        LMIs = cellfun(@(x) x - eye(size(x))*opts.synthesis.zerotol <= 0, LMIs, 'un', 0);
        problems = {LMIs}; 
        
        % 2.2 The "dual" version (Y parameter-varying, X "constant")
        % This is only solved whenever there is (at least) one parameter that has
        % nonzero but bounded rates. For all other cases, the dual
        % formulation equals the primal one, so that it doesn't make sense
        % to solve it twice.
        
        if solvedual
            % build the set of LMIs
            LMIs = cell(length(prob.specs),1); 
        
                % H-infinity performance
                j = 0; 
                for k=1:length(prob.specs) 
                    Bw = prob.Bw(); Bw = Bw(:,prob.specs(k).in); nwk = length(prob.specs(k).in);
                    Cz = prob.Cz(); Cz = Cz(prob.specs(k).out,:); nzk = length(prob.specs(k).out); 
                    Dzw = prob.Dzw(); Dzw = Dzw(prob.specs(k).out,prob.specs(k).in);
                    Dzu = prob.Dzu(); Dzu = Dzu(prob.specs(k).out,:);
                    Dyw = prob.Dyw(); Dyw = Dyw(:,prob.specs(k).in); 
                    LMIs{k} = He([X*prob.A()+Bk*prob.Cy()                   zeros(prob.n())                  zeros(prob.n(),nwk)       zeros(prob.n(),nzk);
                                  Ak'+prob.A()+prob.Bu()*Dk*prob.Cy()       prob.A()*Ys+prob.Bu()*Ck         zeros(prob.n(),nwk)       zeros(prob.n(),nzk);  
                                  (X*Bw+Bk*Dyw)'                            (Bw+prob.Bu()*Dk*Dyw)'           zeros(nwk,nwk)            zeros(nwk,nzk); 
                                  Cz+Dzu*Dk*prob.Cy()                       Cz*Ys+Dzu*Ck                     Dzw+Dzu*Dk*Dyw            zeros(nzk,nzk)]);

                    if isobj(k)   
                        j = j+1;
                        switch opts.synthesis.lyapunovshape 
                            case 1
                              LMIs{k} = LMIs{k} + [zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                                   zeros(prob.n())             -Ysdot                         zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                                   zeros(nwk,prob.n())         zeros(nwk,prob.n())            -sqgamma(j)*eye(nwk)     zeros(nwk,nzk); 
                                                   zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -sqgamma(j)*eye(nzk)];
                            case 2
                              LMIs{k} = LMIs{k} + [zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                                   zeros(prob.n())             -Ysdot                         zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                                   zeros(nwk,prob.n())         zeros(nwk,prob.n())            -sqgamma(j)*eye(nwk)     zeros(nwk,nzk); 
                                                   zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -eye(nzk)];
                            case 3
                              LMIs{k} = LMIs{k} + [zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                                   zeros(prob.n())             -Ysdot                         zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                                   zeros(nwk,prob.n())         zeros(nwk,prob.n())            -eye(nwk)                zeros(nwk,nzk); 
                                                   zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -sqgamma(j)*eye(nzk)];
                            otherwise
                                prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
                        end
                    else
                        LMIs{k} = LMIs{k} + [zeros(prob.n())             zeros(prob.n())                zeros(prob.n(),nwk)      zeros(prob.n(),nzk);
                                             zeros(prob.n())             -Ysdot                         zeros(prob.n(),nwk)      zeros(prob.n(),nzk); 
                                             zeros(nwk,prob.n())         zeros(nwk,prob.n())            -eye(nwk)                zeros(nwk,nzk); 
                                             zeros(nzk,prob.n())         zeros(nzk,prob.n())            zeros(nzk,nwk)           -eye(nzk)];
                    end
                end
            
                % stability constraint
                LMIs = [LMIs ; {-[X eye(prob.n()) ; eye(prob.n()) Ys]}];

                % (D-)stability constraints, if any
                for k=1:length(prob.region)
                    if ~isempty(prob.region(k).M) && ~isempty(prob.region(k).L)
                        LMIs = [LMIs ; 
                                {kron(prob.region(k).L,[Ys            eye(prob.n()); 
                                                        eye(prob.n()) X            ]) 
                                 + kron(prob.region(k).M,[prob.A()*Ys+prob.Bu()*Ck   prob.A()+prob.Bu()*Dk*prob.Cy(); 
                                                          Ak                         X*prob.A()+Bk*prob.Cy()       ])
                                 + kron(prob.region(k).M,[prob.A()*Ys+prob.Bu()*Ck   prob.A()+prob.Bu()*Dk*prob.Cy(); 
                                                          Ak                         X*prob.A()+Bk*prob.Cy()       ])'}];
                    end
                end
            
            LMIs = cellfun(@(x) x.safe_degree_elevation(opts.synthesis.degreeelevations), LMIs, 'un', 0);
            LMIs = cellfun(@(x) x.safe_midpoint_refinement(opts.synthesis.knotinsertions), LMIs, 'un', 0);
            LMIs = cellfun(@(x) x - eye(size(x))*opts.synthesis.zerotol <= 0, LMIs, 'un', 0);
            problems = [problems, {LMIs}]; 
        end
        
    % 3. Parse the LMIs for the solvers
    switch lower(opts.synthesis.gammasolver)
        case 'yalmip'
            sdpvars = ct_synthesis_YALMIP(problems,obj,opti,decvars,opts);
        case 'cvx'
            prob.watchdog.error('Parsing of the LMIs for CVX is not supported (yet). Use YALMIP instead.'); 
        case 'lmilab'
            prob.watchdog.error('Parsing of the LMIs for LMILAB is not supported (yet). Use YALMIP instead.'); 
        otherwise
            prob.watchdog.error('LMI parser for synthesis problem not defined.'); 
    end
    
    % 4. Reconstruct the controller 
    K = ct_controller_reconstruction(sdpvars,prob);
    
    % 5. Return gamma or gamma^2 (depending on Lyapunov shaping formulation)
    for k=1:length(K)
        K(k).sqgamma = sdpvars(k).sqgamma;
    end
    
    % 6. Sort according to performance (lowest first)
    [~,i] = sort([decvars.obj]);
    K = K(i); 

end

%% CONTINUOUS TIME - SYNTHESIS 

function sdpvars = ct_synthesis_YALMIP(problems,obj,opti,decvars,opts)

    % solve both problems
    sols = cell(length(problems),1);
    for k=1:length(problems)
        % build problem
        opti.minimize(obj)
        opti.subject_to(); % reset constraints
        opti.subject_to(problems{k}')
        opti.solver('yalmip',struct('yalmip_options',opts.yalmip));
        
        % solve
        sols{k} = opti.solve();
     
        % retrieve the SDP variables
        f = fieldnames(decvars);
        sdpvars = decvars; 
        for l=1:length(f)
            try
                sdpvars(k).(f{l}) = sols{k}.value(decvars.(f{l}));
            catch % decision variables not part of opti: X and Ys for "primal" formulation, Xs and Y for "dual" formulation
                sdpvars(k).(f{l}) = [];
            end
        end
    end

end

%% CONTINUOUS TIME - RECONSTRUCTION

function K = ct_controller_reconstruction(sdpvars,prob)

    K = struct(); 

    % 1. "Primal" formulation -> always
            if isnumeric(sdpvars(1).Xs) && any(isnan(sdpvars(1).Xs(:))) % infeasible or not properly solved
               K(1).K = [];
               K(1).fb = [];
            else
                if isnumeric(sdpvars(1).Xs) && isnumeric(sdpvars(1).Y) % constant Lyapunov matrices -> state-space form -> K = K1
                    AK = inv(eye(prob.n())-sdpvars(1).Xs*sdpvars(1).Y)*(sdpvars(1).Ak - sdpvars(1).Xs*(prob.A() - prob.Bu()*sdpvars(1).Dk*prob.Cy())*sdpvars(1).Y - sdpvars(1).Bk*prob.Cy()*sdpvars(1).Y - sdpvars(1).X*prob.Bu()*sdpvars(1).Ck); 
                    BK = inv(eye(prob.n())-sdpvars(1).Xs*sdpvars(1).Y)*(sdpvars(1).Bk - sdpvars(1).Xs*prob.Bu()*sdpvars(1).Dk);
                    CK = sdpvars(1).Ck - sdpvars(1).Dk*prob.Cy()*sdpvars(1).Y;
                    DK = sdpvars(1).Dk;
                    K(1).K = gssshinfcd.gsss(AK,BK,CK,DK,0,params);
                    K(1).fb = [];
                else % varying Lyapunov matrices -> LFT form -> K = lft(K1,fb1)
                    At = sdpvars(1).Ak - sdpvars(1).Xs*(prob.A() - prob.Bu()*sdpvars(1).Dk*prob.Cy())*sdpvars(1).Y - sdpvars(1).Bk*prob.Cy()*sdpvars(1).Y - sdpvars(1).Xs*prob.Bu()*sdpvars(1).Ck; 
                    Bt = sdpvars(1).Bk - sdpvars(1).Xs*prob.Bu()*sdpvars(1).Dk;
                    Ct = sdpvars(1).Ck - sdpvars(1).Dk*prob.Cy()*sdpvars(1).Y;
                    Dt = sdpvars(1).Dk;
                    AK = zeros(prob.n());
                    BK = [zeros(prob.n(),prob.ny()) eye(prob.n()) zeros(prob.n())];
                    CK = [zeros(prob.n()+prob.nu(),prob.n());eye(prob.n())];
                    DK = [Dt zeros(prob.nu(),prob.n()) Ct; 
                          Bt sdpvars(1).Xs*sdpvars(1).Y At;
                          zeros(prob.n(),2*prob.n()+prob.ny())];
                    K(1).K = gssshinfcd.gsss(AK,BK,CK,DK,0,prob.P.parameters);
                    K(1).fb = eye(2*prob.n());
                end
            end
            
    % 2. "Dual" formulation -> only if sdpvars has length > 1
        if length(sdpvars)>1
             if isnumeric(sdpvars(2).Ys) && any(isnan(sdpvars(2).Ys(:))) % infeasible or not properly solved
               K(2).K = [];
               K(2).fb = [];
            else
                if isnumeric(sdpvars(2).X) && isnumeric(sdpvars(2).Ys) % constant Lyapunov matrices -> state-space form -> K = K2
                    AK = (sdpvars(2).Ak - sdpvars(2).X*(prob.A() - prob.Bu()*sdpvars(2).Dk*prob.Cy())*sdpvars(2).Ys - sdpvars(2).Bk*prob.Cy()*sdpvars(2).Ys - sdpvars(2).X*prob.Bu()*sdpvars(2).Ck)*inv(eye(prob.n())-sdpvars(2).X*sdpvars(2).Ys); 
                    BK = sdpvars(2).Bk - sdpvars(2).X*prob.Bu()*sdpvars(2).Dk;
                    CK = (sdpvars(2).Ck - sdpvars(2).Dk*prob.Cy()*sdpvars(2).Ys)*inv(eye(prob.n())-sdpvars(2).X*sdpvars(2).Ys);
                    DK = sdpvars(2).Dk;
                    K(2).K = gssshinfcd.gsss(AK,BK,CK,DK,0,params);
                    K(2).fb = [];
                else % varying Lyapunov matrices -> LFT form -> K = lft(K2,fb2)
                    At = sdpvars(2).Ak - sdpvars(2).X*(prob.A() - prob.Bu()*sdpvars(2).Dk*prob.Cy())*sdpvars(2).Ys - sdpvars(2).Bk*prob.Cy()*sdpvars(2).Ys - sdpvars(2).X*prob.Bu()*sdpvars(2).Ck; 
                    Bt = sdpvars(2).Bk - sdpvars(2).X*prob.Bu()*sdpvars(2).Dk;
                    Ct = sdpvars(2).Ck - sdpvars(2).Dk*prob.Cy()*sdpvars(2).Ys;
                    Dt = sdpvars(2).Dk;
                    AK = zeros(prob.n());
                    BK = [zeros(prob.n(),prob.ny()) eye(prob.n()) zeros(prob.n())];
                    CK = [zeros(prob.n()+prob.nu(),prob.n());eye(prob.n())];
                    DK = [Dt zeros(prob.nu(),prob.n()) Ct; 
                          Bt zeros(prob.n()) At;
                          zeros(prob.n(),prob.n()+prob.ny()) sdpvars(2).X*sdpvars(2).Ys];
                    K(2).K = gssshinfcd.gsss(AK,BK,CK,DK,0,prob.P.parameters);
                    K(2).fb = eye(2*prob.n());
                end
            end
        end

end
