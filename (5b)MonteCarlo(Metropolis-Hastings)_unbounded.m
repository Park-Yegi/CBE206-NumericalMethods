% HOMEWORK 5b
% 20160253 Park Yegi

function fx = CBE206_hw5b_20160253(NumMolecules)
    fileID = fopen('MFI.txt', 'r');
    MFI = fscanf(fileID, '%g %g %g', [3 inf]);
    MFI = MFI';
    fclose(fileID);
    NumZeolites = size(MFI, 1);
   
    N = 200000;  % total number of Monte Carlo iteration
    throw_away = 20000;
    EpsilonMtoM = 148;
    SigmaMtoM = 3.73;
    EpsilonMtoZ = 115;
    SigmaMtoZ = 3.47;
    kB = 1;
    Temp = 300; % Temperature in Kelvin
    Rcut = 6.5;
    
    % 1. Random configuration of molecules
    Lx = 20.09;
    Ly = 19.738;
    Lz = 13.142;
    
    for k=1:NumMolecules
        x(k) = Lx*rand;
        y(k) = Ly*rand;
        z(k) = Lz*rand;
    end
    
    
    % 2. Compute initial potential energy
    UTotal = 0; % store total potential energy in UTotal
    for k=1:NumMolecules 
        for kk=k+1:NumMolecules  % avoid double counting
            deltaX = x(k) - x(kk);
            deltaY = y(k) - y(kk);
            deltaZ = z(k) - z(kk);
            r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);  % distance
            if (r <= Rcut)
                UPair = 4 * EpsilonMtoM * (SigmaMtoM^12/r^12 - SigmaMtoM^6/r^6);
                UTotal = UTotal + UPair; % update energy
            end
        end
        
        for kk=1:NumZeolites
            for nx = -1:1
                for ny = -1:1
                    for nz = -1:1
                        xlattice = MFI(kk,1) + nx*Lx;
                        deltaX = x(k) - xlattice;
                        ylattice = MFI(kk,2) + ny*Ly;
                        deltaY = y(k) - ylattice;   
                        zlattice = MFI(kk,3) + nz*Lz;
                        deltaZ = z(k) - zlattice;
                        r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);
                        if r <= Rcut
                            UPair = 4 * EpsilonMtoZ * (SigmaMtoZ^12/r^12 - SigmaMtoZ^6/r^6);
                            UTotal = UTotal + UPair; % update energy
                        end
                    end
                end
            end
        end
    end
    
    % 3. Monte Carlo loop
    for k=1:N
        % 3a. Select a molecule at random
        RandIndex = randi([1, NumMolecules]);
        
        % 3b. Compute its energy
        UOld = 0;
        for kk=1:NumMolecules
            if (kk ~= RandIndex)
                deltaX = x(RandIndex) - x(kk);
                deltaY = y(RandIndex) - y(kk);
                deltaZ = z(RandIndex) - z(kk);
                r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);  % distance
                if (r <= Rcut)
                    UPair = 4 * EpsilonMtoM * (SigmaMtoM^12/r^12 - SigmaMtoM^6/r^6);
                    UOld = UOld + UPair; % update energy
                end
            end
        end
        
        for kk=1:NumZeolites
            for nx = -1:1
                for ny = -1:1
                    for nz = -1:1
                        xlattice = MFI(kk,1) + nx*Lx;
                        deltaX = x(RandIndex) - xlattice;
                        ylattice = MFI(kk,2) + ny*Ly;
                        deltaY = y(RandIndex) - ylattice;
                        zlattice = MFI(kk,3) + nz*Lz;
                        deltaZ = z(RandIndex) - zlattice;
                        r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);
                        if r <= Rcut
                            UPair = 4 * EpsilonMtoZ * (SigmaMtoZ^12/r^12 - SigmaMtoZ^6/r^6);
                            UOld = UOld + UPair; % update energy
                        end
                    end
                end
            end
        end
        
        % 3c. Move a molecule at random
        xnew = x(RandIndex) + (8.0*rand - 4.0);
        ynew = y(RandIndex) + (8.0*rand - 4.0);
        znew = z(RandIndex) + (8.0*rand - 4.0);
        
        % if it moves outside the box, make it come out the other side
        if (xnew < 0)
            xnew = xnew + Lx;
        elseif (xnew > Lx)
            xnew = xnew - Lx;
        end
        if (ynew < 0)
            ynew = ynew + Ly;
        elseif (ynew > Ly)
            ynew = ynew - Ly;
        end
        if (znew < 0)
            znew = znew + Lz;
        elseif (znew > Lz)
            znew = znew - Lz;
        end
        
        % 3d. Compute the new energy
        UNew = 0;
        for kk=1:NumMolecules
            if (kk ~= RandIndex)
                deltaX = xnew - x(kk);
                deltaY = ynew - y(kk);
                deltaZ = znew - z(kk);
                r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);  % distance
                if (r <= Rcut)
                    UPair = 4 * EpsilonMtoM * (SigmaMtoM^12/r^12 - SigmaMtoM^6/r^6);
                    UNew = UNew + UPair; % update energy
                end
            end
        end
       
        for kk=1:NumZeolites     
            for nx = -1:1
                for ny = -1:1
                    for nz = -1:1
                        xlattice = MFI(kk,1) + nx*Lx;
                        deltaX = xnew - xlattice;
                        ylattice = MFI(kk,2) + ny*Ly;
                        deltaY = ynew - ylattice;
                        zlattice = MFI(kk,3) + nz*Lz;
                        deltaZ = znew - zlattice;
                        r = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);
                        if r <= Rcut
                            UPair = 4 * EpsilonMtoZ * (SigmaMtoZ^12/r^12 - SigmaMtoZ^6/r^6);
                            UNew = UNew + UPair; % update energy
                        end
                    end
                end
            end
        end
        
        % 3e. Metropolis-Hastings algorithm
        Metropolis = min(1, exp(-(UNew - UOld)/(kB*Temp)));
        MetroRand = rand;
        
        if (MetroRand <= Metropolis) % if it is inside the box and we accept the move
            x(RandIndex) = xnew;
            y(RandIndex) = ynew;
            z(RandIndex) = znew;
            UTotal = UTotal + (UNew - UOld);
        end
        
        % 3f. Store the values of the current total energy inside an array
        EnergyArray(k) = UTotal;
        
    end % End of Monte Carlo loop
    
    fprintf('average_potential_energy = %f Kelvin\n', sum(EnergyArray(throw_away+1:end))/(N-throw_away));
end
