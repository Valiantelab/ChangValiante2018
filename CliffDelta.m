function d = CliffDelta(X,Y)
% Calculates Cliff's Delta function, a non-parametric effect magnitude
% test. See: http://revistas.javeriana.edu.co/index.php/revPsycho/article/viewFile/643/1092
% for implementation details. 

% calculate length of vetors. 
lx = length(X);
ly = length(Y);

% comparison matrix. First dimension represnt elements in X, the second elements in Y
% Values calculated as follows:
% mat(i,j) = 1 if X(i) > Y(j), zero if they are equal, and -1 if X(i) < Y(j)
mat = zeros(lx, ly);	

% perform all the comparisons. 
for i = 1:lx
	for j = 1:ly
		if X(i) > Y(j)
			mat(i,j) = 1;
		elseif Y(j) > X(i)
			mat(i,j) = -1;
		end
	end
end

% calculate delta. 
d = sum(mat(:)) / (lx * ly)