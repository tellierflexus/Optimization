function simplex(A, b, C, Xbi)
    format rat %To display fraction instead of decimal
    %Xbi = indices of basics variables
    % C is a column vector
    
    step = 1;

    while 1
        disp("Step " + step);
        disp("------------");
        B_1A = inv(A(:,Xbi))*A;
        Xb = inv(A(:,Xbi))*b; %Calculate value of basics variable
        Rc = C' - C(Xbi)'*inv(A(:,Xbi))*A; %We calculate the reduced cost
        f_v = -C(Xbi)'*Xb;
        Smpx_tab = [B_1A, Xb; Rc, f_v];
        disp(Smpx_tab);
        %Ok tricky part we gonna eliminate reduced cost corresponding to non
        %feasible solution
        [~, Nbi] = find(B_1A((Xb==0),:)<0); % We first check for Basic variables equals to 0, then if the corresponding lign in basic direction is non negative
        Rc(Nbi) = Inf; %We set the reduced cost of non feasible direction to inf
        disp("Non feasible direction indices: " + Nbi); %Matlab display this only if Nbi is not empty, don't know why but its working
        [v, i] = min(Rc);
        if v >= 0 % We check if it exists a negative reduced cost (actually if the minimum reduced cost is negative)
           disp("Problem Solved")
           break;
        end
        disp("Variable x" + i + " will enter the basis")
        alphas = Xb./B_1A(:,i); %We calculate the step following the direction to reach 0 for each basic variables
        alphas(alphas<0) = Inf;% We eliminate negative step (i.e the basics variables increase if we follow the direction)
        [~, alpha_i] = min(alphas); %The variable which as the smaller step is choosen to exit the basis
        disp("Variable x" + Xbi(alpha_i) + " will leave the basis")
        Xbi(alpha_i) = i;
        step = step + 1;
    end
    [~, col] = size(A);
    x = zeros(col, 1);
    x(Xbi) = Xb;
    disp("x : ");
    disp(x);
    disp("Value of the function : " + -f_v);
    disp("Problem Solved")
end

