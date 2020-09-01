% Matt Ireland 2020
% Top88 + Better FE

clear all;
close all;
global numNodes numElements numNodesPerElement nodalDOF elementDOF freeDOF numGaussPoints
global nodalCoordinates elementConnectivity dee nodalBoundaryConditions nodalLoads
global length height NXE NYE X_origin Y_origin dhx dhy rawDOF rawNodalDOF
format long g
tic;

%% Inputs

length = 30;% length of the model (mm)
height = 20;% height (mm)

NXE = 2*length;% Number of rows in the x direction
NYE = 2*height;% Number of rows in the y direction

dhx = length/NXE;% Element size in the x direction
dhy = height/NYE;% Element size in the y direction

X_origin = 0;% Keep at 0, Global coordinate reference
Y_origin = 0;% Keep at 0, Global coordinate reference

Force = 1;% N (magnitude)

render = 1;
% 0 for no plot rendering
% 1 to render plots

report = 1;
% 0 for no report
% 1 to render report (plot rendering needs to be on)

materialMode = 1;
% 0 for isotropic
% 1 for transverse isotropic
% 2 for top88 element stiffness matrix

theta = 45;
% Material orientation, degrees CCW from X+

bound = 1;
% Boundary Conditions:
% 0 for fixed
% 1 for midplane pin + rollers
% 2 for MBB

loadMode = 1;
% Loading Conditions:
% 0 for distributed load at free end
% 1 for midplane point load
% 2 for approx parabolic load (mesh dependent)
% 3 MBB point load

loadDirection = 0;
% Loading Direction (LC 0 & 1):
% 0 for Y-
% 1 for X+

%% Optimization Inputs

optimization = 1;% 0 for fea only, 1 for optimization loop

volumeConstraint = 0.1;% Global volume constraint
penalty = 3;% penalization factor (3 for isotropic)

Emin = 1e-9;% Dummy modulus for non-dense elements
rmin = 2.5;% minimum filter radius

filterMethod = 2;
% 0 for no filtering
% 1 for sensitivity filtering
% 2 for density filtering

%% Element and Material Definitions

numNodesPerElement = 4;% number of nodes per element
nodalDOF = 2;% number of nodal degrees of freedom
elementDOF = numNodesPerElement*nodalDOF;% elemental degrees of freedom

% Generate the mesh
Q4_mesh

% material properties
thick = 1;% Beam thickness in mm

%~~~~~~~~~~~~ isotropic ~~~~~~~~~~~~~~~~~~~
if materialMode == 0
    E = 1;     % Elastic modulus in MPa
    vu = 0.3;       % Poisson's ratio
    dee = formdsig(E,vu);
end

%~~~~~~~~~~~~ transverse isotropic ~~~~~~~~
if materialMode == 1
    E1 = 1;        % MPa
    nu12 = 0.27;
    E2 = 0.58;      % MPa
    G12 = 0.20;    % MPa
    Q = ReducedStiffness(E1,nu12,E2,G12);
    dee = OffAxisStiffness(theta,Q);

    E = E1;
end

%~~~~~~~~~~~ top88 element stiffness ~~~~~~
if materialMode == 2
    nu = 0.3;
    E = 1;
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    KE = 1/(1-nu^2)/24 * ([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
end

%% Boundary conditions

% Initialise the matrix nodalBoundaryConditions to 1
nodalBoundaryConditions = ones(numNodes, nodalDOF);
rawNodalDOF = ones(numNodes,nodalDOF);

% Implement boundary conditions according to input toggle
for i=1:numNodes

    % fixed boundary condition
    if bound == 0
        if nodalCoordinates(i,1) == X_origin
            nodalBoundaryConditions(i,:) = [0 0];
        end
    end

    % mid plane pin constrains X&Y, rollers otherwise
    if bound == 1
        if nodalCoordinates(i,1) == X_origin
            nodalBoundaryConditions(i,:) = [0 1];
        end
        if nodalCoordinates(i,2) == height/2 && nodalCoordinates(i,1) == X_origin
            nodalBoundaryConditions(i,:) = [0 0];
        end
    end

    % MBB boundary condition
    if bound == 2
        % Rollers on left edge constrain horizontal
        if nodalCoordinates(i,1) == X_origin
            nodalBoundaryConditions(i,:) = [0 1];
        end
        % Single roller on bottom rightmost edge
        if nodalCoordinates(i,1) == length && nodalCoordinates(i,2) == Y_origin
            nodalBoundaryConditions(i,:) = [1 0];
        end
    end
end

% Count degrees of freedom after boundary conditions
freeDOF = 0;
rawDOF = 1;
for i = 1:numNodes
    for j = 1:nodalDOF
        if nodalBoundaryConditions(i,j) ~= 0
            freeDOF = freeDOF + 1;
            nodalBoundaryConditions(i,j) = freeDOF;
        end
        rawNodalDOF(i,j) = rawDOF;
        rawDOF = rawDOF + 1;
    end
end

%% Loading

% Initialise the matrix of nodal loads to 0
nodalLoads= zeros(numNodes, 2);

j = 1; %index variable for parabolic load
for i=1:numNodes

    % distributed load at free end
    if loadMode == 0
        if nodalCoordinates(i,1) == length
            if loadDirection == 0 % Y-
                nodalLoads(i,:) = [0 -Force/(NYE+1)];
            end
            if loadDirection == 1 % X+
                nodalLoads(i,:) = [(Force/(NYE+1)) 0];
            end
        end
    else
        % midplane point load
        if loadMode == 1
            if nodalCoordinates(i,1) == length && nodalCoordinates(i,2) == height/2
                if loadDirection == 0 % Y-
                    nodalLoads(i,:) = [0 -Force];
                end
                if loadDirection == 1 % X+
                    nodalLoads(i,:) = [Force 0];
                end
            end
        else
            % approx parabolic traction (mesh dependent) Y-
            if loadMode == 2
                if nodalCoordinates(i,1) == length
                    para4 = [0 -302.1 -395.8 -302.1 0];
                    para8 = [0 -84.26 -142.9 -178 -189.7 -178 -142.9 -84.26 0];
                    if NYE == 4
                        nodalLoads(i,:) = [0 para4(j)];
                    end
                    if NYE == 8
                        nodalLoads(i,:) = [0 para8(j)];
                    end
                    j = j + 1;
                end
            else
                % MBB point load
                if loadMode == 3
                    if nodalCoordinates(i,1) == 0 && nodalCoordinates(i,2) == height
                        nodalLoads(i,:) = [0 -Force];
                    end
                end
            end
        end
    end
end

%% Assemble the global force vector

fg = zeros(freeDOF,1);
for i=1: numNodes
    if nodalBoundaryConditions(i,1) ~= 0
        fg(nodalBoundaryConditions(i,1))= nodalLoads(i,1);
    end
    if nodalBoundaryConditions(i,2) ~= 0
        fg(nodalBoundaryConditions(i,2))= nodalLoads(i,2);
    end
end

% Form the matrix containing the abscissas and the weights of Gauss points
numGaussPoints = 4;
samp = gauss(numGaussPoints);

%% Optimization

if optimization == 1

    % Initialize optimization variables
    elementDensity = repmat(volumeConstraint,NYE,NXE);
    physicalDensity = elementDensity;
    complianceGlobal = [];
    volumeGlobal = [];
    changeGlobal = [];
    counter = 0;
    change = 1;

    % Filter Prep
    iH = ones(NXE*NYE*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    % Build Filter Structure
    for i1 = 1:NXE
        for j1 = 1:NYE
            e1 = (i1-1)* NYE + j1;
            for i2 = max(i1 - (ceil(rmin)-1),1) : min(i1 + (ceil(rmin)-1) ,NXE)
                for j2 = max(j1-(ceil(rmin)-1),1) : min(j1+(ceil(rmin)-1) ,NYE)
                    e2 = (i2-1) * NYE + j2;
                    k = k + 1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0, rmin-sqrt( (i1-i2)^2 + (j1-j2)^2) );
                end
            end
        end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);

    % Begin Optimization Loop
    while change > 0.01
        counter = counter + 1;

        % Initialize
        kk = zeros(freeDOF, freeDOF);% Initialise the global stffness matrix
        ke = zeros(elementDOF,elementDOF,numElements);% Initialise the element stiffness matrix to zero
        ka = zeros(elementDOF,elementDOF,numElements);% Initialise the static element stiffness matrix
        ce = zeros(elementDOF,elementDOF,numElements);% Initialize the element compliance matrix to zero
        elementDOFmat = zeros(numElements,8);% Initialize a dof connectivity matrix

        % write element density matrix into a column vector for iteration
        elDens = reshape(physicalDensity,[],1);

        % Numerical integration and assembly of the global stiffness matrix
        for i = 1:numElements
            % coordinates of the nodes of element i, and its steering vector
            [coord,g,h] = elem_q4(i);
            % Write raw nodal DOF vector into connectivity matrix
            elementDOFmat(i,:) = h;

            % Build element stiffness matrix
            for ig = 1:numGaussPoints
                wi = samp(ig,2);
                for jg = 1:numGaussPoints
                    wj = samp(jg,2);
                    % Derivative of shape functions in local coordinates
                    [der,fun] = fmlin(samp,ig,jg);
                    % Compute Jacobian matrix
                    jac = der * coord;
                    % Compute determinant of Jacobian matrix
                    d = det(jac);
                    % Derivative of shape functions in global coordinates
                    deriv = jac\der;
                    % Form matrix [B]
                    bee = formbee(deriv,numNodesPerElement,elementDOF);

                    % Integrate stiffness matrix
                    if counter == 1 % First iteration
                        if materialMode ~= 2 % khennane element stiffness
                            ke(:,:,i) = ke(:,:,i) + d * thick * wi * wj * bee' * dee * bee;
                            ka(:,:,i) = ka(:,:,i) + d * thick * wi * wj * bee' * dee * bee;
                        else % top88 element stiffness
                            ke(:,:,i) = KE ;
                        end
                    end

                    if counter ~= 1 % After first iteration
                        if materialMode ~= 2 % khennane element stiffness
                            ke(:,:,i) = ke(:,:,i) + d * thick * wi * wj * bee' * dee * bee;
                            ka(:,:,i) = ka(:,:,i) + d * thick * wi * wj * bee' * dee * bee;
                        else % top88 element stiffness
                            ke(:,:,i) = KE;
                        end
                    end
                end
            end

            ke(:,:,i) = ke(:,:,i) .* (Emin + elDens(i)^penalty * (E - Emin));

            % Add to global stiffness matrix
            for ii = 1:elementDOF
                if g(ii) ~= 0
                    for jj = 1:elementDOF
                        if g(jj) ~= 0
                            kk(g(ii),g(jj))= kk(g(ii),g(jj)) + ke(ii,jj,i);
                        end
                    end
                end
            end
        end

        % Solve for unknown displacements
        delta = kk\fg ;

        % Write out nodal displacements
        displacement = zeros(numNodes,2);% Initialize the displacement array
        dispForOpt = zeros(numNodes,1);% Initialize the optimization disp array

        j = 2;
        for i = 1:numNodes
            if nodalBoundaryConditions(i,1) == 0
                x_disp =0.;
            else
                x_disp = delta(nodalBoundaryConditions(i,1));
            end

            if nodalBoundaryConditions(i,2) == 0
                y_disp = 0.;
            else
                y_disp = delta(nodalBoundaryConditions(i,2));
            end
            displacement(i,:) = [x_disp  y_disp];

            % write out full displacement vector observing boundary conditions
            dispForOpt(j-1) = x_disp;
            dispForOpt(j) = y_disp;
            j = 2 * (i+1);
        end

        % Generate compliance data structures
        U = zeros(size(elementDOFmat)); UKE = zeros(size(elementDOFmat));
        for ii = 1:size(U,1)
            for jj = 1:size(U,2)
                U(ii,jj) = dispForOpt(elementDOFmat(ii,jj));
            end
            UKE(ii,:) = U(ii,:) * ka(:,:,ii);
        end

        v = UKE .* U;
        V = sum(v,2);

        % Generate elemental compliance values
        complianceElemental = reshape( (sum( UKE .* U ,2)) ,NYE,NXE);
        % Evaluate objective function
        complianceGlobal(counter) = sum( sum( (Emin + physicalDensity .^ penalty * (E - Emin)) .* complianceElemental));
        % Generate elemental compliance sensitivity
        complianceSensitivity = - penalty * (E - Emin) * physicalDensity .^(penalty - 1) .* complianceElemental;
        % Generate elemental volume sensitivity
        volumeSensitivity = ones(NYE,NXE);

        % Write out density sum as global volume value
        volumeGlobal(counter) = sum(physicalDensity(:));

        % Update sensitivity by filtering mode
        if filterMethod == 1
            complianceSensitivity(:) = ...
                H * (elementDensity(:) .* complianceSensitivity(:)) ./ Hs ./ max( 1e-3, elementDensity(:) );
        elseif filterMethod == 2
            complianceSensitivity(:) = H * (complianceSensitivity(:) ./ Hs);
            volumeSensitivity(:) = H * (volumeSensitivity(:) ./ Hs);
        end

        % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        lagrangeOne = 0; lagrangeTwo = 1e9; move = 0.01;
        while (lagrangeTwo-lagrangeOne)/(lagrangeOne+lagrangeTwo) > 1e-3
            lagrangeMid = 0.5*(lagrangeTwo+lagrangeOne);
            newDensity = max(0, ...
                max(elementDensity - move, ...
                min(1, ...
                min(elementDensity + move, ...
                elementDensity .* sqrt(- complianceSensitivity ./ volumeSensitivity / lagrangeMid)))));
            % choose physical density value by filter method
            if filterMethod == 1
                physicalDensity = newDensity;
            elseif filterMethod == 2
                physicalDensity(:) = (H * newDensity(:))./ Hs;
            end

            if sum(physicalDensity(:)) > volumeConstraint*NYE*NXE
                lagrangeOne = lagrangeMid;
            else
                lagrangeTwo = lagrangeMid;
            end
        end
        change = max(abs(newDensity(:) - elementDensity(:)));
        elementDensity = newDensity;
        changeGlobal(counter) = change;

        % Plot current density if rendering is enabled
        if render == 1
            figure(1)
            fig = gcf;
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            fontsize = 18;
            map = colormap('gray');
            map = flip(map);
            % Element density
            subplot(3,2,1)
            imagesc(flip(physicalDensity));
            title('Material Placement','FontSize',fontsize);
            colormap(map);
            colorbar('eastoutside');
            pbaspect([NXE NYE 1]);
            % Compliance Sensitivity
            subplot(3,2,2)
            imagesc(-flip(complianceSensitivity));
            title('Compliance Sensitivity','FontSize',fontsize);
            colorbar('eastoutside');
            colormap(map);
            pbaspect([NXE NYE 1]);
            % Density Histogram
            subplot(3,2,3)
            histogram(physicalDensity,100,'Normalization','probability');
            title('Element Density Occurence','FontSize',fontsize);
            xlabel('Percent Dense');
            ylabel('Normalized Count');
            pbaspect([2 1 1]);
            % Objective Function Convergence
            subplot(3,2,4)
            plot(complianceGlobal);
            title('Obj. Function Value','FontSize',fontsize);
            xlabel('Iteration');
            ylabel('Obj. Fun.');
            pbaspect([2 1 1]);
            % Constraint Plot
            subplot(3,2,5)
            plot(volumeGlobal./numElements);
            title('Constraint Value','FontSize',fontsize);
            xlabel('Iteration');
            ylabel('Constraint');
            pbaspect([2 1 1]);
            % Change
            subplot(3,2,6)
            plot(changeGlobal);
            title('Obj. Function Change','FontSize',fontsize);
            xlabel('Iteration');
            ylabel('Change');
            pbaspect([2 1 1]);
        end

        % Write out to command line
        fprintf(' It.%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n'...
            ,counter,complianceGlobal(counter),mean(physicalDensity(:)),change);

    end % End Optimization Loop

    % Save figure output for report
    t = datetime('now','Format','yyyyMMddHHmm'); s = char(t);
    saveas(gcf,[s '_matl.png']);

    % Save workspace
    save([s '.mat']);
end

%% FEA Only

if optimization == 0

    kk = zeros(freeDOF, freeDOF);% Initialise the global stffness matrix
    ke = zeros(elementDOF,elementDOF,numElements);% Initialise the element stiffness matrix to zero
    elementDOFmat = zeros(numElements,8);% Initialize a dof connectivity matrix

    % Numerical integration and assembly of the global stiffness matrix
    for i=1:numElements
        % coordinates of the nodes of element i, and its steering vector
        [coord,g] = elem_q4(i);
        % Write steering vector into connectivity matrix, write nans to zeros
        elementDOFmat(i,:) = g;

        % Build element stiffness matrix
        for ig=1: numGaussPoints
            wi = samp(ig,2);
            for jg = 1:numGaussPoints
                wj = samp(jg,2);
                % Derivative of shape functions in local coordinates
                [der,fun] = fmlin(samp,ig,jg);
                % Compute Jacobian matrix
                jac = der * coord;
                % Compute determinant of Jacobian matrix
                d = det(jac);
                % Derivative of shape functions in global coordinates
                deriv = jac\der;
                % Form matrix [B]
                bee = formbee(deriv,numNodesPerElement,elementDOF);
                % Integrate stiffness matrix
                ke(:,:,i) = ke(:,:,i) + d * thick * wi * wj * bee' * dee * bee;
            end
        end

        % Add to global stiffness matrix
        for ii=1:elementDOF
            if g(ii) ~= 0
                for jj=1:elementDOF
                    if g(jj) ~= 0
                        kk(g(ii),g(jj))= kk(g(ii),g(jj)) + ke(ii,jj,i);
                    end
                end
            end
        end
    end

    % Solve for unknown displacements
    delta = kk\fg ;

    % Write out nodal displacements
    displacement = zeros(numNodes,2);% Initialize the displacement array
    for i=1: numNodes
        if nodalBoundaryConditions(i,1) == 0
            x_disp =0.;
        else
            x_disp = delta(nodalBoundaryConditions(i,1));
        end

        if nodalBoundaryConditions(i,2) == 0
            y_disp = 0.;
        else
            y_disp = delta(nodalBoundaryConditions(i,2));
        end
        displacement(i,:) = [x_disp  y_disp];
    end
end

%% Stresses and strains at the centre of each element
numGaussPoints = 1;
samp = gauss(numGaussPoints);

% preallocate element stress array
SIGMA = zeros(numElements,3);

for i = 1:numElements
    % coordinates of the nodes of element i, and its steering vector
    [coord,g] = elem_q4(i);
    % Initialise element displacement to zero
    eld = zeros(elementDOF,1);

    for m = 1:elementDOF
        if g(m) == 0
            eld(m) = 0;
        else
            % Retrieve element displacement from the global displacement vector
            eld(m) = delta(g(m));
        end
    end

    for ig = 1:numGaussPoints
        wi = samp(ig,2);
        for jg = 1:numGaussPoints
            wj = samp(jg,2);
            % Derivative of shape functions in local coordinates
            [der,fun] = fmlin(samp,ig,jg);
            % Compute Jacobian matrix
            jac = der * coord;
            % Derivative of shape functions in global coordinates
            %             deriv = inv(jac) * der;
            deriv = jac \ der;
            % Form matrix [B]
            bee = formbee(deriv,numNodesPerElement,elementDOF);
            % Compute strains
            eps = bee * eld;
            % Compute stresses
            sigma = dee * eps;
        end
    end
    % Store stresses for all elements
    SIGMA(i,:) = sigma ;
end

%% Average stresses at nodes, Generate plottable outputs

[ZX, ZY, ZT, Z1, Z2]=stresses_at_nodes_Q4(SIGMA);
U1 = displacement(:,1); U2 = displacement(:,2); % define plottable displacement vectors

%% Plot Outputs
if render == 1

    plotOut = 2;

    for plotOut = 2:6
        f = figure(plotOut);
        renderType = plotOut - 2;

        % Format plot variable and title
        if renderType == 0
            plotVar = ZX;
            plotLabel = 'Sigma XX, MPa';
            plotFileName = [s '_SXX'];
            colormap(parula(10));
        else
            if renderType == 1
                plotVar = ZY;
                plotLabel = 'Sigma YY, MPa';
                plotFileName = [s '_SYY'];
                colormap(parula(10));
            else
                if renderType == 2
                    plotVar = ZT;
                    plotLabel = 'Sigma XY, MPa';
                    plotFileName = [s '_SXY'];
                    colormap(parula(10));
                else
                    if renderType == 3
                        plotVar = U1;
                        plotLabel = 'Displacement U1, mm';
                        plotFileName = [s '_U1'];
                        colormap(winter(10));
                    else
                        if renderType == 4
                            plotVar = U2;
                            plotLabel = 'Displacement U2, mm';
                            plotFileName = [s '_U2'];
                            colormap(winter(10));
                        end
                    end
                end
            end
        end

        % Write title
        title(plotLabel,'FontSize',18);

        % Find plot value bounds
        cmin = min(plotVar);
        cmax = max(plotVar);
        caxis([cmin cmax]);

        % Use patch function to interpolate elemental values
        patch('Faces', elementConnectivity, 'Vertices', nodalCoordinates, 'FaceVertexCData', plotVar, ...
            'Facecolor','interp','Marker','.','EdgeColor','none');

        % Find colorbar bounds and discretize to 10 intervals
        cbh = colorbar('h');
        yticks = cmin:(cmax - cmin)/10:cmax;
        set(cbh,'YTick',yticks)

        % Name demical places on colorbar ticks
        yticks = arrayfun(@(x) sprintf('%.1E',x),yticks,'un',0);
        set(cbh,'TickLabels',yticks,'FontSize',12);

        % Place colorbar on right side
        cbh.Location = 'eastoutside';

        % Axis labels and font sizes
        ax = gca;
        ax.FontSize = 16;
        xlabel('Position, mm','FontSize',16);
        ylabel('Position, mm','FontSize',16);

        % Maintain 1:1 aspect ratio
        pbaspect([length height 1]);

        % Generate an image filename
        saveas(gcf,plotFileName,'png');

    end

end

% Write execute time to command line
ttc = toc;
disp(['Time to complete: ' num2str(ttc)])

%% Generate Report

if report == 1
    % Create filename
    filename = [s '_toFEA'];

    % Writeout vars, simplify datatypes
    h = height; l = length; numel = numElements;

    % Import report API classes
    import mlreportgen.dom.*
    import mlreportgen.report.*

    % Add report container
    rpt = Report(filename,'pdf');

    % Add content to container
    %     Title Page
    titlepg = TitlePage;
    titlepg.Title = 'FEA & TO';
    add(rpt,titlepg);% add title

    %     Input Data Table
    tableStyle = { Width("100%"), ...
        Border("solid"), ...
        RowSep("solid"), ...
        ColSep("solid") };

    headerStyle = { BackgroundColor("LightBlue"), ...
        Bold(true) };

    footerStyle = { BackgroundColor("LightCyan"), ...
        ColSep("none"), ...
        WhiteSpace("preserve") };

    headerContent = {'Field', '', 'Choice'};
    bodyContent = {'Part length', 'mm', l ; ...
        'Part height', 'mm', h ; ...
        'Num. Elements', 'Count', numel ; ...
        'Force', '|N|' , Force ; ...
        'Theta' , 'Deg' , theta ; ...
        'Modulus' , 'MPa' , E ; ...
        ' ' , ' ' , ' ' ; ...
        'Volume Constraint' , ' ' , sprintf('%1.2f',volumeConstraint) ; ...
        'Obj. Fun. Val.' , ' ' ,  sprintf('%2.2e',complianceGlobal(end)) ; ...
        'Filter Radius' , 'mm' , rmin ; ...
        'Lower Modulus Limit' , 'MPa' , sprintf('%1.1e',Emin) ; ...
        'Iterations' , 'Count' , counter; ...
        ' ' , ' ' , ' ' ; ...
        'BC' , 'Mode' , bound ; ...
        'Material' , 'Mode' , materialMode ; ...
        'Loading' , 'Mode' , loadMode ; ...
        'Direction' , 'Mode' , loadDirection ; ...
        'Filter Method' , 'Mode' , filterMethod; ...
        ' ' , ' ' , ' ' ; ...
        'Time to Complete' , 'S' , num2str(ttc) ; ...
        };

    tableContent = [headerContent; bodyContent];

    add(rpt, Heading1("Input Settings Used"));% Add table heading

    tbl = Table(tableContent);
    tbl.Style = tableStyle;

    firstRow = tbl.Children(1);
    firstRow.Style = headerStyle;

    tbl.TableEntriesHAlign = "center";

    add(rpt,tbl);% Add table

    %     Input Chapter Heading
    chap1 = Chapter('Optimization Output');
    add(rpt,chap1);

    %     Input material placement figure
    matl = Image([s '_matl.png']);
    matl.Style = [matl.Style {ScaleToFit}];
    add(rpt,matl);

    %     Input Chapter Heading
    chap2 = Chapter('FEA Plots After Optimization');
    add(rpt,chap2);

    %     Input sxx fea figure
    img = Image([s '_SXX.png']);
    img.Style = [img.Style {ScaleToFit}];
    add(rpt,img);

    %     Input syy fea figure
    img = Image([s '_SYY.png']);
    img.Style = [img.Style {ScaleToFit}];
    add(rpt,img);

    %     Input sxy fea figure
    img = Image([s '_SXY.png']);
    img.Style = [img.Style {ScaleToFit}];
    add(rpt,img);

    %     Input u1 fea figure
    img = Image([s '_U1.png']);
    img.Style = [img.Style {ScaleToFit}];
    add(rpt,img);

    %     Input u2 fea figure
    img = Image([s '_U2.png']);
    img.Style = [img.Style {ScaleToFit}];
    add(rpt,img);

    % Close the report
    close(rpt);
    % Display the report
    rptview(rpt);
end

%% Functions Called
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Function: Mesh Generation
function Q4_mesh

% This module generates a mesh of linear quadrilateral elements

global numNodes numElements
global nodalCoordinates elementConnectivity
global NXE NYE X_origin Y_origin dhx dhy

numNodes = 0;
k = 0;

for i = 1:NXE
    for j = 1:NYE
        k = k + 1;
        n1 = j + (i-1)*(NYE + 1);
        nodalCoordinates(n1,:) = [(i-1)*dhx - X_origin    (j-1)*dhy - Y_origin ];
        n2 = j + i*(NYE+1);
        nodalCoordinates(n2,:) = [i*dhx - X_origin       (j-1)*dhy - Y_origin  ];
        n3 = n1 + 1;
        nodalCoordinates(n3,:) = [(i-1)*dhx - X_origin       j*dhy - Y_origin  ];
        n4 = n2 + 1;
        nodalCoordinates(n4,:) = [i*dhx- X_origin       j*dhy - Y_origin       ];
        numElements = k;
        elementConnectivity(numElements,:) = [n1  n2  n4  n3];
        numNodes = n4;
    end
end
end

%% Function: Nodal Coordinates
function[coord,g,h] = elem_q4(i)

% This function returns the coordinates of the nodes of
% element i and its steering vector g

global numNodesPerElement nodalDOF rawNodalDOF
global nodalCoordinates elementConnectivity nodalBoundaryConditions

l = 0;
coord = zeros(numNodesPerElement,nodalDOF);
for k = 1: numNodesPerElement
    for j = 1:nodalDOF
        coord(k,j) = nodalCoordinates(elementConnectivity(i,k),j);
        l = l+1;
        g(l) = nodalBoundaryConditions(elementConnectivity(i,k),j);
        h(l) = rawNodalDOF(elementConnectivity(i,k),j);
    end
end
end

%% Function: Shape Function
function[der,fun] = fmlin(samp,ig,jg)

% This function returns the vector of the shape function and their
% derivatives with respect to xi and eta

xi = samp(ig,1);
eta = samp(jg,1);

fun = 0.25*[(1.- xi - eta + xi*eta);...
    (1.+ xi - eta - xi*eta);...
    (1.+ xi + eta + xi*eta);...
    (1.- xi + eta - xi*eta)];

der = 0.25*[-(1-eta)    (1-eta)    (1+eta)   -(1+eta);...
    -(1-xi)     -(1+xi)    (1+xi)     (1-xi)];
end

%% Function: Strain-Displacement Relation
function[bee] = formbee(deriv,numNodesPerElement,elementDOF)
%
%  This function assembles the matrix [bee] from the
%  derivatives of the shape functions in global coordinates
%
bee = zeros(3,elementDOF);
for m = 1:numNodesPerElement
    k = 2*m;
    l = k-1;
    x = deriv(1,m);
    bee(1,l) = x;
    bee(3,k) = x;
    y = deriv(2,m);
    bee(2,k) = y;
    bee(3,l) = y;
end
end

%% Function: Plane Stress Stiffness Matrix
function[dee] = formdsig(E,vu)
%
% This function forms the elasticity matrix for a plane stress problem

c = E/(1.-vu*vu);
dee = c*[1 vu 0; vu 1 0; 0 0 .5*(1.-vu)];

end

%% Function: Gauss Quadrature
function[samp] = gauss(numGaussPoints)
%
% This function returns the abscissas and weights of the Gauss points for numGaussPoints equal up to 4
%

samp = zeros(numGaussPoints,2);
if numGaussPoints == 1
    samp = [0.  2];
elseif numGaussPoints == 2
    samp = [1./sqrt(3)   1.;...
        -1./sqrt(3)  1.];
elseif numGaussPoints == 3
    samp = [.2*sqrt(15.)   5./9; ...
        0.            8./9.;...
        -.2*sqrt(15.)   5./9];
elseif numGaussPoints == 4
    samp = [0.861136311594053       0.347854845137454; ...
        0.339981043584856       0.652145154862546; ...
        -0.339981043584856       0.652145154862546; ...
        -0.861136311594053       0.347854845137454];
end
end

%% Function: Average Nodal Stress
function[ZX, ZY, ZT, Z1, Z2] = stresses_at_nodes_Q4(SIGMA)
%
% This function averages the stresses at the nodes
%
global numNodes numElements numNodesPerElement elementConnectivity

for k = 1:numNodes
    sigx = 0; sigy = 0; tau = 0;
    ne = 0;
    for iel = 1:numElements
        for jel=1:numNodesPerElement
            if elementConnectivity(iel,jel) == k
                ne = ne+1;
                sigx = sigx + SIGMA(iel,1);
                sigy = sigy + SIGMA(iel,2);
                tau = tau + SIGMA(iel,3);
            end
        end
    end
    ZX(k,1) = sigx/ne;
    ZY(k,1) = sigy/ne;
    ZT(k,1)=tau/ne;
    Z1(k,1)= ((sigx+sigy)/2 + sqrt(((sigx+sigy)/2)^2 +tau^2))/ne;
    Z2(k,1)= ((sigx+sigy)/2 - sqrt(((sigx+sigy)/2)^2 +tau^2))/ne;
end
end

%% Function: Reduced Stiffness Matrix
function Q = ReducedStiffness(E1,nu12,E2,G12)

nu21 =  E2 / E1 * nu12;
Q11 =  E1 / (1 - nu12 * nu21);
Q12 = nu12 * E2 / (1 - nu12 * nu21);
Q22 =  E2 / (1 - nu12 * nu21);
Q66 = G12;
Q = [Q11 Q12 0;Q12 Q22 0; 0 0 Q66];

end

%% Function: Off-axis Stiffness Matrix
function QBar = OffAxisStiffness(theta,Q)

% Calculate the cosine and sine of the off-axis angle theta specified in degrees
m = cosd(theta);
n = sind(theta);
Q11 = Q(1,1); Q12=Q(1,2); Q22=Q(2,2); Q66=Q(3,3);

QBar11 = Q11 * m^4 + 2 * (Q12 + 2* Q66) * m^2 * n^2 + Q22 * n^4;
QBar12 = (Q11 + Q22 - 4* Q66) * n^2 * m^2 + Q12 * (n^4 + m^4);
QBar16 = (Q11 - Q12 - 2* Q66) * n * m^3 + (Q12 - Q22 + 2* Q66) * n^3 * m;
QBar22 = Q11 * n^4 + 2 * (Q12 + 2* Q66) * n^2 * m^2 + Q22 * m^4;
QBar26 = (Q11 - Q12 - 2* Q66) * n^3 * m + (Q12 - Q22 + 2*Q66) * n * m^3;
QBar66 = (Q11 + Q22 - 2* Q12 - 2* Q66) * n^2 * m^2 + Q66 * (n^4 + m^4);

QBar =[QBar11 QBar12 QBar16; QBar12 QBar22 QBar26; QBar16 QBar26 QBar66];

end
