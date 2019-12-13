%% Circular Cross-Section

% The definition of D in the pipe friction equation is "Hydraulic Diameter"
% Hydraulic diameter is defined as:
%   D = 4*cross_sectional_area / perimeter
%  This document proves the minimum value of D for a fixed CSA is for a
%  circle. With a fixed CSA we need to reduce perimeter
% The least surface area occurs for regular polygons. "Star" shaped ones
% tend to have lots of perimeter and little area, so we shall consider
% regular polygons as the best of the categories of 2D shapes
% Polygon with sides N and side length L, all angles in degrees

%% Area of polygon's triangles
% If you take each polygon and draw lines from each corner to the centre
% you have split it into N identical triangles with bottom length L
% For clarity of orientation, the length L is oriented downwards
% Either side of side L are two identical angles, B, that are both equal to
% half the interior angle of the polygon

    %interior_angle = 180 - (360/N);
    %B = interior_angle/2;

% This leaves the top angle as

    %A = 180 - 2*B;  % which simplifies to 360/N

% The two sides heading towards the centre of the polygon both have equal
% side lengths x, using Law of sines

    %x = L sind(A)/sind(B);

% Using Side-angle-side we get the area of a triangle as

    %area = x*L*sind(A) ==> L^2*sind^2(A)/sind(B);

% Total CSA is the angle of each triangle * N (number of triangles)

    %CSA = N*area;

% Perimeter is length of sides * number of sides

    %P = N*L;

% Thus we can put this into hydraulic diameter equation (ln 5)

    %D = 4*N*area/(N*L) ==> 4*area/L;
    %   ==> 4*L*sind^2(A)/sind(B);
  
% To minimise D we need to maximise sind(B) and minimise sind(A)
% This is done by minimising A and maximising B
    % A = 360/N which is decreased as N -> infinity
    % B = 90 - (180/N) which is increased as N -> infinity
% Thus as number of sides tends to infinity, hydraulic diameter for fixed
% CSA decreases.
% A regular polygon with infinite sides is a circle

% Thus a circle is the most efficient shape for a pipe

%% Graphs
% Visualisation of hydraulic diameter for side length 1 for N = 1:360
% Note, N=1,2 don't create a real shape

clear all
clc
% L = 1

for N = 1:360
    interior_angle(N) = 180 - (360/N);
    B(N) = interior_angle(N)/2;
    A(N) = 360/N;
    x(N) = sind(A(N))/sind(B(N));
    D(N) = 4*(sind(A(N))^2)/sind(B(N));
end

figure
plot(D)
title('Hydraulic Diameter vs Number of Sides')
ylabel('Hydraulic Diameter')
xlabel('N sides')
