function K = solve(obj) 
% SOLVE Solves the control problem

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

    import gssshinfcd.util.*;

    if isempty(obj.opts)
        obj = obj.setoptions(); 
    end
    
    % 0. Start
    obj.watchdog.info(bold('THIS IS GSSSHINFCD.'));
    
    % 1. Formulate the problem as an SDP
    obj.watchdog.info('STARTED PREPROCESSING');
    obj.watchdog.startSubTask();
        obj.watchdog.info('Formulating the SDP...');
        obj.preprocessor = gssshinfcd.preprocessor(obj.problem, obj.opts);
    obj.watchdog.endSubTask();
    
    % 2. Solve the SDP
    obj.watchdog.info('STARTED SOLVING');
    obj.watchdog.startSubTask();
        obj.watchdog.info('Solving the SDP...');
        K = obj.preprocessor.solve();
    obj.watchdog.endSubTask();
    
    % 3. Postprocess the solution
    obj.watchdog.info('STARTED POSTPROCESSING...');
    obj.postprocessor = gssshinfcd.postprocessor(obj.problem, K);
    obj.watchdog.startSubTask();
        obj.watchdog.info('Compensating for algebraic loops...'); 
        obj.postprocessor = obj.postprocessor.compensateDyu();
    obj.watchdog.endSubTask(); 
    
    % 4. Return
    K = obj.postprocessor.K; 
    
    % 5. All done
    obj.watchdog.info('ALL DONE!'); 
    
end