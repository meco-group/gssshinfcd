function formatmsg = bold(msg)
% BOLD Puts a message to be displayed in bold
    
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

    if verLessThan('matlab','7.13') || (usejava('jvm') && ~feature('ShowFigureWindows'))
        formatmsg = msg; 
    else
        formatmsg = ['<strong>' msg '</strong>'];
    end
end