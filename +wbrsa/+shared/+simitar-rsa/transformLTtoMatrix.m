%   This file is part of Simitar
%
%   Simitar is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%   Simitar is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Simitar.  If not, see <http://www.gnu.org/licenses/>.
%

function [M] = transformLTtoMatrix(D,k,includeDiagonal)

idx = 1;
M = zeros(k,k);

if includeDiagonal
  for ik = 1:k
    len = k - ik + 1;
    range = idx:(idx+len-1);
    M(ik:k,ik) = D(range);
    idx = idx + len;
  end
else
  for ik = 1:(k-1)
    len = k - ik;
    range = idx:(idx+len-1);
    M((ik+1):k,ik) = D(range);
    idx = idx + len;
  end
  
end
