function [Ti,g]=LM_calc_norm_str(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)

e0=8.85418782*1e-12; % Permittivity κενού-αέρα
m0=4*pi*1e-7; % Permeability κενού-αέρα

for k=1:1:num_files %  With increasing frequency!
    %k  Pointer to inform in the command window how often we are
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % Synthesis of Z 'for the present frequency - (ord x ord)  
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Synthesis of Y 'for the present frequency - (ord x ord)  
    end 
    
    if k==1 % QR decomposition is applied to the first frequency because this is how the algorithm starts and because it is 0 Hz where the product under consideration cannot be normalized normally - Attention !!! Manually can be applied to the second frequency and so on.

        S=Ydis*Zdis; % Product Y'Z' - (ord x ord)
        [Tidis,ldis]=eig(S); % QR decomposition - All tables are (ord x ord)
        gdis=sqrt(diag(ldis)); %Finding  propagation constant - (1 x ord)

        g(k,:)=gdis; % Save  propagation constant - (num_files x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Save the vector table - (num_files x ord^2)
        end

    else % The algorithm is applied to the other frequencies Levenberg-Marquardt (LM)
        S=Ydis*Zdis; % Product Y'Z' - (ord x ord)
        norm=-((2*pi*freq(k)).^2)*e0*m0; % Norm which normalizes the product Y'Z 'so that the eigenvalues are easier to calculate - Attention !!! The eigenvalues that will be calculated need to be reversed with the same norm while the eigenvectors remain unchanged since they are multiplication by a fixed number.
        %S=S./norm; % Normalization of the product Y'Z' - (ord x ord) 
        S=(S./norm)-O;
        
        gdis=g(k-1,:); % Introduction of the propagation constant of the previous frequency - Seeding!! - (1 x ord)
        ldis=gdis.^2; % Calculation of eigenvalues - (1 x ord)
        %ldis=ldis./norm; % Normalization of eigenvalues
        ldis=(ldis./norm)-1;
        
        for o=1:ord
            Tidis(o,:)=Ti(k-1,(o-1)*ord+1:o*ord); % Insert the table of eigenvectors of the previous frequency - Seeding!! - (ord x ord)
        end
               
        for o=1:ord
            x0=[real(Tidis(:,o));imag(Tidis(:,o));real(ldis(o));imag(ldis(o))]; % The vector x contains first the real and then the imaginary part of an eigenvector, and then the real and imaginary part of the corresponding eigenvalue - (2*ord+2 x 1)
            Q=[real(S) -imag(S);imag(S) real(S)]; % Auxiliary table for the implementation of the object function (2*ord x 2*ord)
            obj=@(x)[(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord));sum((x(1:ord).^2)-(x(ord+1:2*ord).^2))-1;sum(x(1:ord).*x(ord+1:2*ord))-0]; % Object function - Contains the problem of bidding plus two limiting equations for the measure of eigenvectors! Returns a real vector (2*ord+2 x 1)
                        
            options=optimset('Algorithm',{'levenberg-marquardt',0.0001},'Display','none','ScaleProblem','Jacobian'); % Various parameters to solve the object-function problem
            %options=optimset('Algorithm',{'levenberg-marquardt',1e-4},'Display','final','ScaleProblem','Jacobian','TolFun',1e-10,'TolX',1e-10); 
            
            x=fsolve(obj,x0,options); %Find the root of the object function - Returns the vector x that achieves the zero of the object function - (2*ord+2 x 1)
            
            Tidis(:,o)=x(1:ord)+1i*x(ord+1:2*ord); % Export the table of eigenvectors - (ord x ord)
            ldis(o)=x(2*ord+1)+1i*x(2*ord+2); % Extraction of the corresponding eigenvalue - (ord x ord)
        end
        
        %ldis=ldis.*norm; % Reversal of normalization - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % Calculation and storage propagation constant - (1 x ord)
               
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Save the vector table - (num_files x ord^2)
        end
        
    end
end
            
        