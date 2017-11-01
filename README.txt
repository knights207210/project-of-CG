README
Author: Xu Han <xuha2442@colorado.edu>
Time Consuming: 7hours

---------------------------------------------------------------------------
Brife description
---------------------------------------------------------------------------

This openGL program is used to display a 3D scene with textures under lighting in the 3 modes: Overhead orthogonal; Overhead perspective; and First person navigation
---------------------------------------------------------------------------
The code is programed under OS X.
---------------------------------------------------------------------------

Run 'make' to build the code. 

---------------------------------------------------------------------------
Running the Program
---------------------------------------------------------------------------
 Key bindings:
 *  0      mode 0: Display the 3D scene in orthogonal mode
 *  1      mode 1: Display the 3D scene in projection mode
 *  2      mode 2: Display the 3D scene with first person navigation
 *  l          Toggles lighting
 *  a/A        Decrease/increase ambient light
 *  d/D        Decrease/increase diffuse light
 *  s/S        Decrease/increase specular light
 *  e/E        Decrease/increase emitted light
 *  n/N        Decrease/increase shininess
 *  F1         Toggle smooth/flat shading
 *  F2         Toggle local viewer mode
 *  F3         Toggle light distance (1/5)
 *  m          Toggles light movement
 *  []         Lower/rise light
 *  </>        Change the positions of light source mannually while the light 
        	movement is stopped
 *  +/-        Change field of view of perspective
 *  u/U        Make the first person navigation up and down 
 *  b/B        Change the scale of the scene in first person navigation(b: larger; B: 		smaller)
 *  x          Toggle axes
 *  arrows     Change view angle(<-/->: rotate around Y axis; up/down: rotate around X axis) when in mode 0 and mode 1
               Change the position and view angle(<-/->, make the first person navigation look left and right; up and down, 
               make the first person navigation move forward and backward)
 *  F5/F6      Zoom in and out for orthogonal
 *  r          Reset view angle
 *  ESC        Exit 

---------------------------------------------------------------------------
Source Files and Directory Structure
---------------------------------------------------------------------------

Usage: openGL

README              -- this file

Makefile            -- makefile for the project

hw6.c        —  the source code

