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

function [D] = transformMatrixToLT(M,includeDiagonal)

k = size(M,1);
idx = 1;

if includeDiagonal
  m = k*(k+1)/2;
  D = zeros(m,1);
  
  for ik = 1:k
    len = k - ik + 1;
    range = idx:(idx+len-1);
    D(range) = M(ik:k,ik);
    idx = idx + len;
  end
  
else
  m = k*(k-1)/2;
  D = zeros(m,1);
  
  for ik = 1:(k-1)
    len = k - ik;
    range = idx:(idx+len-1);
    D(range) = M((ik+1):k,ik);
    idx = idx + len;
  end
  
end
