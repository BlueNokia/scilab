function [x] = gausselimPP(A,b)
	[nA,mA] = size(A)
	[nb,mb] = size(b)
	if nA<>mA then
		error('gausselim - Matrix A must be square');
		abort;
	elseif mA<>nb then
		error('gausselim - incompatible dimensions between A and b');
		abort;
	end;
	a = [A b];
	n = nA;
	// Augmented matrix
	// Matrix size
	//Forward elimination with partial pivoting
	for k=1:n-1
		kpivot = k; amax = abs(a(k,k));
		//Pivoting
		for i=k+1:n
			if abs(a(i,k))>amax then
				kpivot = i; amax = a(k,i);
			end;
		end;
		temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
		for i=k+1:n
			//Forward elimination
			for j=k+1:n+1
				a(i,j)=a(i,j)-a(k,j)*a(i,k)/a(k,k);
			end;
		end;
	end;
		//Backward substitution
	x(n) = a(n,n+1)/a(n,n);
	for i = n-1:-1:1
		sumk=0
		for k=i+1:n
			sumk=sumk+a(i,k)*x(k);
		end;
		x(i)=(a(i,n+1)-sumk)/a(i,i);
	end;
endfunction;
