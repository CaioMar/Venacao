# Venacao
Analise de dados mestrado

Outdated!!!
Check MAIN() function and do README update

	Leaf Venation Program Documentation and Manual

1 - Seed (for the pseudo-random number generator);
2 - lstepgrowth (How much is the leaf boundary growing in each step? ~ growthrate*characterisc length);
3 - Complementary Venation step multiplier (How many simulation steps are spent generating the second venation? stepmultiplier*first venation max number of steps);
4 - Growthstepmax (Maximum number of simulation steps employed in the generation of the first venation pattern);
5 - Veindist (length of the edge connecting the new venation node to its parent node);
6 - Birthdist (A new auxin node is considered valid if no existing node is inside the circle with radius equal to birthdist centered at the new auxin node);
7 - Killdist (if a vein node is placed inside the circle with radius killdist centered at an auxin node, that auxin node is removed from the leaf);
8 - Maxtries (Maximum number of failed attempts to place an auxin within a grid cell considered acceptable);
9 - Squarehex (accepts 2 inputs: 0 - square leaf boundary and 1 - hexagonal leaf boundary);
10 - denst (Density of auxin nodes in the leaf);
11 - numdiv (the square leaf boundary is partitioned in numdiv*numdiv cells);
12 - squareside (initial length of a side of the square leaf);
13 - finallength (maximum final square leaf side length accepted - it CAN be less);
14 - maxsteps2 (Number of steps in the initiation of the sencond venation);
15 - bdmultiplier (Multiplies the birthdist used in the first venation and uses the resultant birthdist in the second venation generation);
16 - kdmultiplier (Multiplies the killdist used in the first venation and uses the resultant killdist in the second venation generation);
17 - denstmultiplier (Multiplies the denst used in the first venation and uses the resultant denst in the second venation generation);
18 - veindistmultiplier (Multiplies the veindist used in the first venation and uses the resultant veindist in the second venation generation);
19 - auxincutoff (number of auxin sources left before stopping increase in the second venation);
20 - murrayexp (murray exponent that defines relations between the parent and children vein diameters);
21 - printvein (bool that defines how the output is going to be displayed 0 - simple display: all vein orders have the same color; 1 - veins of different orders are displayed in different colors
22 - background color control (1 - white 0 - green)

there are more parameters that need to be inputed. Will update the description when possible

Good tested parameters

Classic configuration (2 seeds - 1 placed at the middle of the base of square boundary and 1 at the top): 

tested parameters
./leafvenationsim.exe 6 0.13 3 65 0.148 0.68 0.79 3 0 190 65 0.1 50 4 0.72 0.72 0.72 0.78 100 2.42 1 1 0 2 1.7 0 0 3 4

Compiling the program
g++ leafvenationsim.cc -lm -lGLU -lglut -lGL mtrand.cpp gl2ps.c leafvenationsim2.hh -o leafvenationsim.exe

Dependent libraries:
Glew
Glut
OpenGL

Works on Ubuntu 16.04 and 14.04 (tested)
