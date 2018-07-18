function [Area, Perimeter, Lattice, HeatMap, AFit,PFit] = PoreMC(NumIterations, MetalProb, ChalcProb) %Analysis of Nanopore Growth
%% Monte Carlo simulation of pore growth treating S atoms as separate entities occupying the same lattice site
%Inputs:
%NumIterations is a guess of the number of iterations used for initializing the area, perimeter, and heatmap outputs and determining x-axis scale in video figure.
%MetalProb is a number between 0 and 1 which corresponds to the probability that any edge metal atom will be ejected
%ChalcProb is same as MetalProb, but for chalcogen atoms

%Example inputs: PoreMC(1000, 0.03, 0.03)

%Outputs:
%Area is an array. First column is iteration number, second column is number of vacant atomic columns.
%Perimeter is same as Area, but for the number of edge atoms
%Lattice is the final map of where the remaining atoms are located.
%HeatMap labels which iteration caused each ejection to occur.
%AFit is a quadratic fit of the Area
%PFit is a linear fit of the Perimeter

%% Create Hexagonal Lattice on a Cartesian grid & Initialize HeatMap
%Armchair direction is x-axis  and 60 degree rotations of it
%Zigzag direction is y-axis and 60 degree rotations of it
%Lattice value of 2 if 2 chalcogen atoms present, 1.5 if 1 chalcogen atom present, 1 = metal atom present, 0 = no atom

Lattice = zeros(100,100); %intialize lattice
for i = 1:50 %rows/2
    for j = 1:25 %columns/4
        Lattice(2*i,4*j-1:4*j)=[2,1]; %X = 2, M = 1
        Lattice(2*i-1,4*j-3:4*j-2)=[2,1]; %X = 2, M = 1
    end
end

HeatMap=zeros(100,100); %initialize HeatMap
for i = 1:100
    for j = 1:100
        if Lattice(i,j)>0
            HeatMap(i,j)=-1; %if Heatmap(i,j) is never modified later in the code, we know this atom survives until the end of the simulation
        else
            HeatMap(i,j)=NaN; %empty index (not part of lattice)
        end
    end
end

%% Begin Measuring Lattice Properties and Initialize Pore (Vacancy)
TotalAtoms = zeros(NumIterations,1); %initialize count of total atoms (note: propto area)
TotalEdgeAtoms = zeros(NumIterations,1); %initialize count of edge atoms (note: propto perimeter)
Area = zeros(NumIterations,2); %initialize area
Perimeter = zeros(NumIterations,2);

i=1; %iteration number
Lattice(50,51)=0; %initialize pore near middle
HeatMap(50,51)=i; %position 50, 51 removed after first iteration

TotalAtoms(1,1) = sum(sum(Lattice~=0)); %store atom count at t=0 --> pore = single atom vacancy
Area(i,1)=i; %iteration number
Area(i,2)=TotalAtoms(1)+1-TotalAtoms(i); %area after ith iteration

TotalEdgeAtoms(1,1) = sum(sum(isEdgeAtom(Lattice)~=0)); %store edge atom count at t=0 --> pore = single atom vacancy
Perimeter(i,1)=i; %iteration number
Perimeter(i,2)=TotalEdgeAtoms(i); %perimeter after ith iteration

%% MP4 Video Formatting
if MetalProb<.099
	Percent = ['0',num2str(100.*MetalProb)];
else
	Percent = num2str(100.*MetalProb);
end
cm = [0, 0, 0; 0, 0, 0; 84, 181, 182; 255, 200, 50; 255, 255, 0]; % Color Map for values 0, 0.5, 1, 1.5, 2 corresponding to "empty", null, metal, 1chalc, 2chalc
cm = cm./255; %convert to decimal
v=VideoWriter(['PoreGrowth_p', num2str(Percent),'.mp4'], 'MPEG-4'); %initialize video writer
v.Quality = 100; %set quality to max
v.FrameRate = 3; % adjust based on total number of frames and desired runtime
open(v) %open video writer

fig = figure; %open figure
subplot(1,2,1) %left subplot
P = pcolor(Lattice); %heat map
P.LineStyle = 'none'; %no lines
Text = text(50,105,['P = ', Percent, '%'],'FontSize',14,'HorizontalAlignment', 'center');
colormap(cm); %use colormap cm defined above
caxis([0 2]); %set color axis
colorbar off %remove colorbar
pbaspect([1 1 1]) %set 1:1 aspect ratio
set(gca,'nextplot','replacechildren','visible','off') %turn off axes
subplot(1,2,2) %right plot
pbaspect([1 1 1]) %set 1:1 aspect ratio
xlabel('Iterations') %label x axis
xlim([0 NumIterations]) %number of iterations
yyaxis left %left y-axis
plot(Perimeter(i,1), Perimeter(i,2),'b.'); %plot perimeter data up until ith iteration
ylabel('Perimeter') %label left y-axis
ylim([0 35]) %fix y-axis limits
yyaxis right %right y-axis
plot(Area(i,1), Area(i,2),'r.'); %plot area data up until ith iteration
ylabel('Area') %label right y-axis
ylim([0 1600]) %fix y-axis limits
f = getframe(fig); %get frame
writeVideo(v,f.cdata); %add frame data to video
i = i+1;

%% Iterate Pore Growth and Store Lattice Properties at Each Step
while TotalAtoms(1)+1-TotalAtoms(i-1)<1500 %pore occupies 30% of atomic sites
	[Lattice, HeatMap] = MonteCarloPoreGrowth(i, Lattice, HeatMap, MetalProb, ChalcProb); %Pore growth happens here. New Lattice and HeatMap replace old ones.
    TotalAtoms(i,1)=sum(sum(Lattice~=0)); %add new atom count data point to the list
    TotalEdgeAtoms(i,1) = sum(sum(isEdgeAtom(Lattice)~=0)); %add new edge atom count data point to the list
    
    Area(i,1)=i; %same area and perimeter steps as above in Begin Measuring Lattice Properties and Initialize Pore (Vacancy) section
    Area(i,2)=TotalAtoms(1)+1-TotalAtoms(i);
	Perimeter(i,1)=i;
    Perimeter(i,2)=TotalEdgeAtoms(i);
    
    subplot(1,2,1) %same plotting steps as above in MP4 Video Formatting section
    P = pcolor(Lattice);
    P.LineStyle = 'none';
    Text = text(50,105,['P = ', Percent, '%'],'FontSize',14,'HorizontalAlignment', 'center');
    colormap(cm);
    caxis([0 2]);
    colorbar off
    pbaspect([1 1 1])
    set(gca, 'visible', 'off');
    subplot(1,2,2)
    pbaspect([1 1 1])
    xlabel('Iterations')
    xlim([0 NumIterations])
    yyaxis left
    plot(Perimeter(1:i,1), Perimeter(1:i,2),'b.');
    ylabel('Perimeter')
    ylim([0 150])
    yyaxis right
    plot(Area(1:i,1), Area(1:i,2),'r.');
    ylabel('Area')
    ylim([0 1600])
    f= getframe(fig);
    writeVideo(v,f.cdata);
    
    i=i+1; %increase iteration counter
end
close(v)

%% Clean up
HeatMap(HeatMap == -1) = NaN; %any atoms that survive until the end of the simulation should still have values of -1. To visualize the pore better, we can remove them now
Area = Area(1:max(Area(:,1)),:); %remove any extra 0's at the end if NumIterations overestimated
Perimeter = Perimeter(1:max(Perimeter(:,1)),:); %remove any extra 0's at the end if NumIterations overestimated

%% Fit Area and Perimeter Outputs
if ~isempty(Area(Area(:,2)>50,1)) %if Area is ever bigger than 50 (there tends to be a lot of statistical noise below 50 and fitting is improved if we ignore that)
    AFit=fit(Area(Area(:,2)>50,1),Area(Area(:,2)>50,2),'poly2'); %only use data from Area > 50 to do a quadratic fit
    PFit=fit(Perimeter(Area(:,2)>50,1),Perimeter(Area(:,2)>50,2),'poly1'); %only use data from Area > 50 to do a linear fit
else
    AFit.p1=0;
	PFit.p1=0;
end
end


%% isEdgeAtom Function
%Edge of the pore
function Map = isEdgeAtom(Lattice) %0 if not EdgeAtom, 1 if two-neighbor edge atom, 2 if one-neighbor edge atom, 3 if zero-neighbor edge atom
[row, col] = size(Lattice);
Map = zeros(row,col); %initialize output Boolean map

for i = 2:row-1 %ignore perimeter row and col of Lattice because they of course won't be saturated
    for j = 2:col-1
        if Lattice(i,j)>0 && sum(sum(Lattice(i-1:i+1,j-1:j+1)~=0))==3 %if there is an atom there and it has two nearest neighbors
            Map(i,j)=1; %if two-neighbor edge atom (missing 1 neighbor)
        elseif Lattice(i,j)>0 && sum(sum(Lattice(i-1:i+1,j-1:j+1)~=0))==2 %if there is an atom there and it has one nearest neighbor
            Map(i,j)=2; %if one-neighbor edge atom (missing 2 neighbors)
        elseif Lattice(i,j)>0 && sum(sum(Lattice(i-1:i+1,j-1:j+1)~=0))==1 %if there is an atom there and it has zero nearest neighbors
            Map(i,j)=3; %if zero-neighbor edge atom (missing 3 neighbors)
        end
    end
end
end

%% MonteCarloPoreGrowth Function
function [Lattice, HeatMap] = MonteCarloPoreGrowth(IterationNumber, Lattice, HeatMap, MetalProb, ChalcProb)
[row, col] = size(Lattice); %generalize function by determining size of input lattice
EdgeLattice = isEdgeAtom(Lattice); %1 if 2-neighbor edge atom. 2 if 1-neighbor edge atom. 3 if 0-neighbor edge atom. 0 if not an edge atom

for i = 2:row-1
    for j = 2:col-1
        r1 = rand(1);
        r2 = rand(1);
        
        %Island
        %Not perfect -- only applies if it is a single-atom island and not if there is a floating cluster
        %Could be improved with a simple message passing algorithm, but may significantly slow down the code
        if EdgeLattice(i,j)==3 %0-neighbor edge atom
            Lattice(i,j)=0; %disconnected floating atoms are always displaced (100% probability)
            HeatMap(i,j)=IterationNumber; %site becomes vacancy on iteration number i
        
        %Metal
        elseif EdgeLattice(i,j)==1 && Lattice(i,j)==1 %2-neighbor Metal
            if r1 < MetalProb
                Lattice(i,j)=0;
                HeatMap(i,j)=IterationNumber;
            end
        elseif EdgeLattice(i,j)==2 && Lattice(i,j)==1 %1-neighbor Metal
            %2-neighbor and 1-neighbor are separated in case one wants to treat them differently
            if r1 < (MetalProb)
                Lattice(i,j)=0;
                HeatMap(i,j)=IterationNumber;
            end
            
        %Dichalcogen
        elseif EdgeLattice(i,j)==1 && Lattice(i,j)==2 %2-neighbor dichalcogen
            if (r1 < ChalcProb || r2 < ChalcProb) %if at least one of the atoms is removed
                Lattice(i,j)=1.5;
                if (r1 < ChalcProb && r2 < ChalcProb) %if both atoms are removed
                    Lattice(i,j)=0;
                    HeatMap(i,j)=IterationNumber; %only add iteration number to this index of the HeatMap once atomic column is fully vacant
                end
            end
        elseif EdgeLattice(i,j)==2 && Lattice(i,j)==2 %1-neighbor dichalcogen
            %2-neighbor and 1-neighbor are separated in case one wants to treat them differently
            if (r1 < ChalcProb || r2 < ChalcProb) %if at least one of the atoms is removed
                Lattice(i,j)=1.5;
                if (r1 < ChalcProb && r2 < ChalcProb) %if both atoms are removed
                    Lattice(i,j)=0;
                    HeatMap(i,j)=IterationNumber; %only add iteration number to this index of the HeatMap once atomic column is fully vacant
                end
            end
            
        %Monochalcogen
        elseif EdgeLattice(i,j)==1 && Lattice(i,j)==1.5 %2-neighbor monochalcogen
            if r1 < ChalcProb
                Lattice(i,j)=0;
                HeatMap(i,j)=IterationNumber;
            end
        elseif EdgeLattice(i,j)==2 && Lattice(i,j)==1.5 %1-neighbor monohalcogen
            %2-neighbor and 1-neighbor are separated in case one wants to treat them differently
            if r1 < (ChalcProb)
                Lattice(i,j)=0;
                HeatMap(i,j)=IterationNumber;
            end
        end
    end
end

end
