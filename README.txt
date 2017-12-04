README
Author: Xu Han <xuha2442@colorado.edu>

---------------------------------------------------------------------------
Brife description
---------------------------------------------------------------------------

This openGL program is for CSCI5229’s course project. It has two modes: the projection mode and first person navigation mode. Currently, users could make the following amazing things happen:

1. users could observe the water simulation in Forbidden Forest(press ‘w’/‘W’ in first person mode);
2. users could start fireworks using the wand(press ‘f’/‘F’ in first person mode);
3. users could go inside the room where basic scene and fire simulation is finished( press ‘q’/‘Q’ in first person mode or change the view anger in projection mode; when observe the fire simulation, users could press ‘c’ or ‘C’ to switch off the other light source);
4. users could change between a night and day mode;
5. 3D selection will happen when in the place to hit the ball(press ’t’/’T’ in first person mode). Users could click left mouse to select 3D objects: the flying broomstick—-help users to fly up; the ball—-the ball will fly around waiting to be hit; the bat—users could use the bat to hit the ball(when using the bat, press ‘F7’/‘F8’ to change the hitting angle, press ‘F4’ to hit).
6. light could move together with the first person view(press ‘i’/‘I’ in the first person mode to toggle the light).

What will make my project worthy of an A:
1. water wave simulation;
2. particle engines and light modeling of fire
3. 3D selection
4. ball’s moving trace calculation when hitting the ball
5. drawing the scene and light modeling of the scene.
---------------------------------------------------------------------------
The code is programed under OS X.
---------------------------------------------------------------------------

Run 'make' to build the code. 

---------------------------------------------------------------------------
Running the Program
---------------------------------------------------------------------------
 Key bindings:
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
 *  F4         hit the ball
 *  F7/F8      Change the angle of the bat when hitting the ball.
 *  F9/F10     Change the position of eye up and down in first person mode.
 *  m          Toggles light movement
 *  []         Lower/rise light
 *  </>        Change the positions of light source manually while the light 
        	movement is stopped
 *  +/-        Change field of view of perspective
 *  u/U        Make the first person navigation up and down 
 *  b/B        Change the scale of the scene in first person navigation(b: larger; B: 		smaller)
 *  c/C        Switch on/off the light source which goes around
               (for fire observation);
 *  f/F        When in the first person mode, use the wand to start the firework
 *  w/W        When in the first person mode, go to the place where the water
               simulation is
 *  t/T        When in the first person mode, go to the place where we can hit the
               ball
 *  q/Q        When in the first person mode, go to the place where has the fire
               simulation and indoor scene
 *  i/I        Toggle the light which moves together with stick
 *  k/K        Change between the day/night mode
 *  x          Toggle axes
 *  arrows     Change view angle(<-/->: rotate around Y axis; up/down: rotate around
               X axis) when in mode 1 Change the position and view angle(<-/->, make 
               the first person navigation look left and right; up and down, make the
               first person navigation move forward and backward)
 * LEFT MOUSE  When in the place to hit the ball, click left mouse to select 3D
               objects: the flying broomstick,the ball or the bat.
 *  r          Reset view angle
 *  ESC        Exit 

---------------------------------------------------------------------------
Source Files and Directory Structure
---------------------------------------------------------------------------

Usage: openGL

README              -- this file

Makefile            -- makefile for the project

hw7.cpp        —  the source code

