/*  Assignment 7: Xu Han
 *  Project
 *
 *
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
 *  F9/F10     Change the position of eye up and down in first person mode.
 *  m          Toggles light movement
 *  []         Lower/rise light
 *  </>        Change the positions of light source manually while the light 
          movement is stopped
 *  +/-        Change field of view of perspective
 *  u/U        Make the first person navigation up and down 
 *  b/B        Change the scale of the scene in first person navigation(b: larger; B:     smaller)
 *  f/F        When in the first person mode, use the wand to start the firework
 *  w/W        When in the first person mode, go to the place where the water simulation is
 *  q/Q        When in the first person mode, go to the place where has the fire simulation and indoor scene
 *  k/K        Change between the day/night mode
 *  x          Toggle axes
 *  arrows     Change view angle(<-/->: rotate around Y axis; up/down: rotate around X axis) when in mode 1
               Change the position and view angle(<-/->, make the first person navigation look left and right; up and down, 
               make the first person navigation move forward and backward)
 *  r          Reset view angle
 *  ESC        Exit 

 */
#include "CSCIx229.h"
#include <vector>
#include <iostream>
#define PI 3.14159265359
#define NUM_X_OSCILLATORS   150
#define NUM_Z_OSCILLATORS   150
#define NUM_OSCILLATORS     NUM_X_OSCILLATORS*NUM_Z_OSCILLATORS
#define OSCILLATOR_DISTANCE   0.05

#define OSCILLATOR_WEIGHT       0.0002
#define SQR(a) (a*a)

#define BILLBOARDING_NONE             0
#define BILLBOARDING_PERPTOVIEWDIR          1 
#define BILLBOARDING_PERPTOVIEWDIR_BUTVERTICAL    2  
#define RANDOM_FLOAT (((float)rand())/RAND_MAX)

int stick_light =0;
clock_t start;
double t_bat;
int hit =1;
int zh_bat=0;
int click_broom = 1;
int click_ball = 1;
int click_bat = 1;
int objID;
int width_win = 400;
int height_win = 400;
float mouseX =0;
float mouseY =0;
int zh=0;       // rotate around z
int NumOfEdges=50;   //to make up tower's circles and number of teeth on tower
float LowerHeight=1.0;    // tower's parameters
float HigherHeight=5.0;   // tower's parameters
float Radius=1.5;         // tower's parameters
float TowerTeethSize=0.5; // tower's parameters
float WallTeethSize=0.7;  // wall's teeth width
float WallHeight=3.5;     // wall's parameters
float HouseHeight=7;      // house's parameters 
float HouseCenterHeight=8;// house's parameters
float scale=0.2;  
float r=0.0;    
int plane=1; 
float speed=0.5;
int    sky[2];   //  Sky textures
int box=0;
int firework=0;
int fire;
int move_light = 1;

//initialized point view

float Ex = 0.0;
float Ey = 0.0;
float Ez = 0.0;

float s_atx=0.0;
float s_aty=0.0;
float s_atz=0.0;

float theta=0;

int axes=1;       //  Display axes
int mode=0;       //  Projection mode
int mode_project =1;
int move=1;       //  Move light
int th=160;         //  Azimuth of view angle
int ph=35;         //  Elevation of view angle
int fov=90;       //  Field of view (for perspective)
int light=1;      //  Lighting
double asp=1;     //  Aspect ratio
double dim=5.0;   //  Size of world
// Light values
int one       =   1;  // Unit value
int distance  =   5;  // Light distance
int inc       =  10;  // Ball increment
int smooth    =   1;  // Smooth/Flat shading
int local     =   0;  // Local Viewer Model
int emission  =   0;  // Emission intensity (%)
int ambient   =  0;  // Ambient intensity (%)
int diffuse   = 100;  // Diffuse intensity (%)
int specular  =   0;  // Specular intensity (%)
int shininess =   0;  // Shininess (power of two)
float shiny   =   1;  // Shininess (value)
int zh_l        =  90;  // Light azimuth
float ylight  =   0;  // Elevation of light
int numberOfTrees = 30;
int obj;

unsigned int texture[13];  //texture names

#define Broom       100               // This is the object ID for the SUN  
#define Ball     101               // This is the object ID for the EARTH
#define Bat     102 


/*
 *  Draw vertex in polar coordinates with normal
 */
static void Vertex(double th,double ph)
{
   double x = Sin(th)*Cos(ph);
   double y = Cos(th)*Cos(ph);
   double z =         Sin(ph);
   //  For a sphere at the origin, the position
   //  and normal vectors are the same
   glNormal3d(x,y,z);
   glVertex3d(x,y,z);
}
/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius (r)
 */
static void ball(double x,double y,double z,double r)
{
  float white[] = {1,1,1,1};
   float Emission[]  = {0.01*emission,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   int th,ph;
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball
   glColor3f(1,1,1);
   //  Bands of latitude
   for (ph=-90;ph<90;ph+=inc)
   {
      glBegin(GL_QUAD_STRIP);
      for (th=0;th<=360;th+=2*inc)
      {
         Vertex(th,ph);
         Vertex(th,ph+inc);
      }
      glEnd();
   }
   //  Undo transofrmations
   glPopMatrix();
}

/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius (r)
 */
static void ball_brown(double x,double y,double z,double r)
{
  float white[] = {1,1,1,1};
   float Emission[]  = {0.01*emission,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   int th,ph;
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball
  // glColor3f(1,1,1);
   glColor3f(0.543,0.270,0.074);
   //  Bands of latitude
   for (ph=-90;ph<90;ph+=inc)
   {
      glBegin(GL_QUAD_STRIP);
      for (th=0;th<=360;th+=2*inc)
      {
         Vertex(th,ph);
         Vertex(th,ph+inc);
      }
      glEnd();
   }
   glColor3f(1,1,1);
   //  Undo transofrmations
   glPopMatrix();
}
static void cube1(double x,double y,double z,
                 double dx,double dy,double dz,
                 double th)
{
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1,1,1);
   glBindTexture(GL_TEXTURE_2D,texture[11]);
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glRotated(th,0,1,0);
   glScaled(dx,dy,dz);
   //  Cube
   glBegin(GL_QUADS);
   //  Right
   //glColor3f(0,0.749,1);
   glNormal3f(0,0,1);
   glTexCoord2f(0,0);glVertex3f(-1,-1, 1);
   glTexCoord2f(1,0);glVertex3f(+1,-1, 1);
   glTexCoord2f(1,1);glVertex3f(+1,+1, 1);
   glTexCoord2f(0,1);glVertex3f(-1,+1, 1);
   //  Left
   glNormal3f( 0, 0,-1);
   glTexCoord2f(0,0);glVertex3f(+1,-1,-1);
    glTexCoord2f(1,0);glVertex3f(-1,-1,-1);
    glTexCoord2f(1,1);glVertex3f(-1,+1,-1);
    glTexCoord2f(0,1);glVertex3f(+1,+1,-1);
   //  Front

  // glColor3f(0,0.49,1);
   glNormal3f(+1, 0, 0);
   glTexCoord2f(0,0);glVertex3f(+1,-1,+1);
    glTexCoord2f(1,0);glVertex3f(+1,-1,-1);
    glTexCoord2f(1,1);glVertex3f(+1,+1,-1);
    glTexCoord2f(0,1);glVertex3f(+1,+1,+1);
   //  back
   //glColor3f(0.69,0.769,0.871);
   glNormal3f(-1, 0, 0);
  glTexCoord2f(0,0);glVertex3f(-1,-1,-1);
    glTexCoord2f(1,0);glVertex3f(-1,-1,+1);
    glTexCoord2f(1,1);glVertex3f(-1,+1,+1);
    glTexCoord2f(0,1);glVertex3f(-1,+1,-1);
   //  Top
   //glColor3f(0.541,0.169,0.886);
   glNormal3f( 0,+1, 0);
   glTexCoord2f(0,0);glVertex3f(-1,+1,+1);
    glTexCoord2f(1,0);glVertex3f(+1,+1,+1);
    glTexCoord2f(1,1);glVertex3f(+1,+1,-1);
    glTexCoord2f(0,1);glVertex3f(-1,+1,-1);
   //  Bottom
   glNormal3f( 0,-one, 0);
   glTexCoord2f(0,0);glVertex3f(-1,-1,-1);
    glTexCoord2f(1,0);glVertex3f(+1,-1,-1);
    glTexCoord2f(1,1);glVertex3f(+1,-1,+1);
    glTexCoord2f(0,1);glVertex3f(-1,-1,+1);
   //  End
   glEnd();
   //  Undo transformations
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
 }



void drawBroom(){
  float RGBA[] = {1,1,1,1};
  float Emission[]  = {0.0,0.0,0.01*emission,1.0};
  glMaterialf(GL_FRONT,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT,GL_SPECULAR,RGBA);
   glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
  obj = LoadOBJ("Broom.obj");
  glCallList(obj);
}

void drawBat(){
  float RGBA[] = {1,1,1,1};
  float Emission[]  = {0.0,0.0,0.01*emission,1.0};
  glMaterialf(GL_FRONT,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT,GL_SPECULAR,RGBA);
   glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
  obj = LoadOBJ("bat.obj");
  glCallList(obj);
}

 void RenderScene() 
{
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
  //glLoadIdentity();                 // Reset The matrix

  // We make our position a bit high and back to view the whole scene

    //    Position      View     Up Vector
  //gluLookAt(Ex, Ey, Ez,     0, 0, 0,     0, 1, 0);   // This determines where the camera's position and view is

  glInitNames();                    // This clears the name stack so we always start with 0 names.

  //GLUquadricObj *pObj = gluNewQuadric();        // Get a new Quadric off the stack

  //gluQuadricTexture(pObj, true);            // This turns on texture coordinates for our Quadrics

  // Bind the sun texture to the sun quadratic
  //glBindTexture(GL_TEXTURE_2D, texture[0]);     // Bind the Sun texture to the sun

  // Below we call glPushName().  We need to pass in an ID that we can check later that will
  // be associated with the polygons drawn next.  Here is how it works.  We call glPushName()
  // and pass an ID.  This pushes this ID onto the name stack.  Then we draw any primitives 
  // or shapes, then we call glPopName() which stops assigning polys to that object name.  
  // We now have a group of polygons that are given an ID.  Our ID SUN now refers to the 
  // sun Quadric we draw below.

 glPushName(Broom);                  // Push on our broom label 
  
  glPushMatrix();                   // Push on a new matrix scope
    //glTranslatef(0, 0, 0);              // Translate this sphere to the left
    glTranslatef(30,1.0,-5.0);
    //gluSphere(pObj, 0.5f, 20, 20);          // Draw the sunwith a radius of 0.5
    if(click_broom)
    drawBroom();
   // ball(0,0,0,0.1);
  glPopMatrix();                    // End the current scope of this matrix

  // Now that we drew the sun, we want to end our Object name.  We call glPopName()
  // to do that.  Now nothing else will be associated with the SUN ID.

  glPopName();                    // Stop assigning polygons to the SUN label (IMPORTANT)

  // Next, we want to bind the Earth texture to our Earth sphere
  //glBindTexture(GL_TEXTURE_2D, texture[1]);

  // Once again, we want to create a object ID for our earth, so we push on the EARTH ID.
  // Now, when we draw the next sphere, it will be associated with the EARTH ID.

  glPushName(Ball);                  // Push on our EARTH label (IMPORTANT)

  // Once again, we want to pop on a new matrix as not to affect any other spheres.
  // We rotate the sphere by its current rotation value FIRST before we translate it.
  // This makes it rotate around the origin, which is where the sun is.
  // Then we rotate it again about the Y-axis to make it spin around itself.

  glPushMatrix();                   // Push on a new matrix scope   
    //glTranslatef(-0.2, 0, 0);             // Translate this sphere to the left
    
    //gluSphere(pObj, 0.2f, 20, 20);          // Draw the sphere with a radius of 0.2 so it's smaller than the sun
    if(click_ball){
      glTranslatef(33,1.0,-5.0);
      //glColor3f(0.543,0.270,0.074);
    ball_brown(1,0,0,0.2);
    //glColor3f(1,1,1);
  
  }
    
  glPopMatrix();                    // End the current scope of this matrix

  // We are done assigning the EARTH object, so we need 
  // to stop assigning polygons to the current ID.

  glPopName();                    // Stop assigning polygons to the EARTH label (IMPORTANT)

  // Bind the pluto texture to the last sphere
  //glBindTexture(GL_TEXTURE_2D, texture[2]);

  // Finally, we want to be able to click on Pluto, so we need a pluto ID.

  glPushName(Bat);                  // Push on our PLUTO label (IMPORTANT)

  // Like we did with the earth, we rotate Pluto around the sun first,
  // then we translate it farther away from the sun.  Next, we rotate Pluto
  // around the Y axis to give it some spin.

  // Pass in our empty glBegin()/glEnd() statement because we are using Quadrics.
  // If we don't do this when using glLoadName(), it will grind to a hault on some cards.
  //glBegin(GL_LINES);
  //glEnd();

  glPushMatrix();                   // Push on a new matrix scope
    glTranslatef(25,1.0,-5.0);
    //gluSphere(pObj, 0.3f, 20, 20);          // Draw the sphere with a radius of 0.1 (smallest planet)
    //ball(2,0,0,0.1);
    if(click_bat)
    drawBat();
  glPopMatrix();                    // End the current scope of this matrix

  // We are finished with the PLUTO object ID, so we need to pop it off the name stack.

  glPopName();                    // Stop assigning polygons to our PLUTO label (IMPORTANT)

  //glutSwapBuffers();                 // Swap the backbuffers to the foreground

  //gluDeleteQuadric(pObj);               // Free the Quadric

}


// 3D selection
//This takes the cursor's x and y position and returns the closet object
//--------------------------------------------------------------------------------------------------------------------------------------

int RetrieveObjectID(int x, int y)
{
  int objectsFound = 0;               // This will hold the amount of objects clicked
  int viewportCoords[4] = {0};            // We need an array to hold our view port coordinates

  // This will hold the ID's of the objects we click on.
  // We make it an arbitrary number of 32 because openGL also stores other information
  // that we don't care about.  There is about 4 slots of info for every object ID taken up.
  unsigned int selectBuffer[32] = {0};        
                            
  // glSelectBuffer is what we register our selection buffer with.  The first parameter
  // is the size of our array.  The next parameter is the buffer to store the information found.
  // More information on the information that will be stored in selectBuffer is further below.

  glSelectBuffer(32, selectBuffer);         // Setup our selection buffer to accept object ID's

  // This function returns information about many things in OpenGL.  We pass in GL_VIEWPORT
  // to get the view port coordinates.  It saves it like a RECT with {top, left, bottom, right}

  glGetIntegerv(GL_VIEWPORT, viewportCoords);     // Get the current view port coordinates

  // Now we want to get out of our GL_MODELVIEW matrix and start effecting our
  // GL_PROJECTION matrix.  This allows us to check our X and Y coords against 3D space.

  glMatrixMode(GL_PROJECTION);            // We want to now effect our projection matrix
  
  glPushMatrix();                   // We push on a new matrix so we don't effect our 3D projection

    // This makes it so it doesn't change the frame buffer if we render into it, instead, 
    // a record of the names of primitives that would have been drawn if the render mode was
    // GL_RENDER are now stored in the selection array (selectBuffer).

    glRenderMode(GL_SELECT);            // Allows us to render the objects, but not change the frame buffer

    glLoadIdentity();               // Reset our projection matrix

    // gluPickMatrix allows us to create a projection matrix that is around our
    // cursor.  This basically only allows rendering in the region that we specify.
    // If an object is rendered into that region, then it saves that objects ID for us (The magic).
    // The first 2 parameters are the X and Y position to start from, then the next 2
    // are the width and height of the region from the starting point.  The last parameter is
    // of course our view port coordinates.  You will notice we subtract "y" from the
    // BOTTOM view port coordinate.  We do this to flip the Y coordinates around.  The 0 y
    // coordinate starts from the bottom, which is opposite to window's coordinates.
    // We also give a 2 by 2 region to look for an object in.  This can be changed to preference.

    gluPickMatrix(x, viewportCoords[3] - y, 20, 20, viewportCoords);

    // Next, we just call our normal gluPerspective() function, exactly as we did on startup.
    // This is to multiply the perspective matrix by the pick matrix we created up above. 

    gluPerspective(45.0f,(float)width_win/(float)height_win,0.1f,150.0f);
    
    glMatrixMode(GL_MODELVIEW);           // Go back into our model view matrix
  
    RenderScene();                  // Now we render into our selective mode to pinpoint clicked objects
    // If we return to our normal render mode from select mode, glRenderMode returns
    // the number of objects that were found in our specified region (specified in gluPickMatrix())

    objectsFound = glRenderMode(GL_RENDER);     // Return to render mode and get the number of objects found

    glMatrixMode(GL_PROJECTION);          // Put our projection matrix back to normal.
  glPopMatrix();                    // Stop effecting our projection matrix

  glMatrixMode(GL_MODELVIEW);             // Go back to our normal model view matrix

  // PHEW!  That was some stuff confusing stuff.  Now we are out of the clear and should have
  // an ID of the object we clicked on.  objectsFound should be at least 1 if we found an object.

  if (objectsFound > 0)
  {   
    // If we found more than one object, we need to check the depth values
    // of all the objects found.  The object with the LEAST depth value is
    // the closest object that we clicked on.  Depending on what you are doing,
    // you might want ALL the objects that you clicked on (if some objects were
    // behind the closest one), but for this tutorial we just care about the one
    // in front.  So, how do we get the depth value?  Well, The selectionBuffer
    // holds it.  For every object there is 4 values.  The first value is
    // "the number of names in the name stack at the time of the event, followed 
    // by the minimum and maximum depth values of all vertices that hit since the 
    // previous event, then followed by the name stack contents, bottom name first." - MSDN
    // The only ones we care about are the minimum depth value (the second value) and
    // the object ID that was passed into glLoadName() (This is the fourth value).
    // So, [0 - 3] is the first object's data, [4 - 7] is the second object's data, etc...
    // Be carefull though, because if you are displaying 2D text in front, it will
    // always find that as the lowest object.  So make sure you disable text when
    // rendering the screen for the object test.  I use a flag for RenderScene().
    // So, lets get the object with the lowest depth!   

    // Set the lowest depth to the first object to start it off.
    // 1 is the first object's minimum Z value.
    // We use an unsigned int so we don't get a warning with selectBuffer below.
    unsigned int lowestDepth = selectBuffer[1];

    // Set the selected object to the first object to start it off.
    // 3 is the first object's object ID we passed into glLoadName().
    int selectedObject = selectBuffer[3];

    // Go through all of the objects found, but start at the second one
    for(int i = 1; i < objectsFound; i++)
    {
      // Check if the current objects depth is lower than the current lowest
      // Notice we times i by 4 (4 values for each object) and add 1 for the depth.
      if(selectBuffer[(i * 4) + 1] < lowestDepth)
      {
        // Set the current lowest depth
        lowestDepth = selectBuffer[(i * 4) + 1];

        // Set the current object ID
        selectedObject = selectBuffer[(i * 4) + 3];
      }
    }
    
    // Return the selected object
    return selectedObject;
  }

  // We didn't click on any objects so return 0
  return 0;                     
}


//waterwave simulation
//-------------------------------------------------------------------------------------------------------------------
struct SOscillator
{
  GLfloat x,y,z;
  GLfloat nx,ny,nz;  //normal vector
  GLfloat UpSpeed;
  GLfloat newY;
  bool bIsExciter;
  //only in use, if bIsExciter is true:
  float ExciterAmplitude;  
  float ExciterFrequency;
};
//vertex data for the waves:
SOscillator * Oscillators;
int NumOscillators;  //size of the vertex array
std::vector <GLuint> IndexVect;  //we first put the indices into this vector, then copy them to the array below
GLuint * Indices;
int NumIndices;   //size of the index array

float g_timePassedSinceStart = 0.0f;  //note: this need not be the real time
bool  g_bExcitersInUse = true;


SF3dVector F3dVector ( GLfloat x, GLfloat y, GLfloat z )
{
  SF3dVector tmp;
  tmp.x = x;
  tmp.y = y;
  tmp.z = z;
  return tmp;
}
SF3dVector AddF3dVectors (SF3dVector* u, SF3dVector* v)
{
  SF3dVector result;
  result.x = u->x + v->x;
  result.y = u->y + v->y;
  result.z = u->z + v->z;
  return result;
}
void AddF3dVectorToVector ( SF3dVector * Dst, SF3dVector * V2)
{
  Dst->x += V2->x;
  Dst->y += V2->y;
  Dst->z += V2->z;
}

GLfloat GetF3dVectorLength( SF3dVector * v)
{
  return (GLfloat)(sqrt(v->x*v->x+v->y*v->y+v->z*v->z));  
}
SF3dVector CrossProduct (SF3dVector * u, SF3dVector * v)
{
  SF3dVector resVector;
  resVector.x = u->y*v->z - u->z*v->y;
  resVector.y = u->z*v->x - u->x*v->z;
  resVector.z = u->x*v->y - u->y*v->x;
  return resVector;
}
SF3dVector Normalize3dVector( SF3dVector v)
{
  SF3dVector res;
  float l = GetF3dVectorLength(&v);
  if (l == 0.0f) return F3dVector(0.0f,0.0f,0.0f);
  res.x = v.x / l;
  res.y = v.y / l;
  res.z = v.z / l;
  return res;
}
SF3dVector operator+ (SF3dVector v, SF3dVector u)
{
  SF3dVector res;
  res.x = v.x+u.x;
  res.y = v.y+u.y;
  res.z = v.z+u.z;
  return res;
}
SF3dVector operator- (SF3dVector v, SF3dVector u)
{
  SF3dVector res;
  res.x = v.x-u.x;
  res.y = v.y-u.y;
  res.z = v.z-u.z;
  return res;
}
SF3dVector operator* (SF3dVector v, float r)
{
  SF3dVector res;
  res.x = v.x*r;
  res.y = v.y*r;
  res.z = v.z*r;
  return res;
}
SF3dVector operator/ (SF3dVector v, float r)
{
  return v*(1/r);
}

//partical engine
//-----------------------------------------------------------------------------------------------------------------------------------------------
#define NULL_VECTOR F3dVector(0.0f,0.0f,0.0f)
class CCCParticleSystem;
class CCCParticle
{
private:
  //Position of the particle:  NOTE: This might be the global position or the local position (depending on the system's translating behavior)
  SF3dVector m_Position;
  //Moving the particle:
  //Particle's velocity:
  SF3dVector m_Velocity;
  //Particle's acceleration (per Sec):
  SF3dVector m_Acceleration;
  //Spinning the particle
  float    m_fSpinAngle; //radian measure
  //Particle's spin speed:
  float    m_fSpinSpeed;
  //Particle's spin acceleration:
  float      m_fSpinAcceleration;
  //Particle's alpha value (=transparency)
  float      m_fAlpha;
  float    m_fAlphaChange;  //how much is the alpha value changed per sec?
  //Particle's color:
  SF3dVector m_Color;   //x=r, y=g, z=b
  SF3dVector m_ColorChange;  //how to change to color per sec
  //Particle's size:  (the use of this value is dependent of m_bUseTexture in the parent!)
  float    m_fSize;
  float    m_fSizeChange; 
  //Handling the lifetime:
  float    m_fDieAge; //At what "age" will the particle die?
  float    m_fAge;    //Age of the particle (is updated   
  //Needed to access the system's values:
  CCCParticleSystem * m_ParentSystem;
  
public:
  bool     m_bIsAlive;  //Is the particle active or not? Must be visible for the System
  void Initialize(CCCParticleSystem * ParentSystem);
  void Update(float timePassed);  //called by UpdateSystem if particle is active
  void Render();  

};


/*******************
CCCParticleSystem
*******************/

class CCCParticleSystem
{
public:  //The values how to emit a particle must be public because the particle  
     //must be able to access them in the creation function.
//*************************************
// EMISSION VALUES
//*************************************

  //Position of the emitter:
  SF3dVector    m_EmitterPosition;
  //How far may the particles be created from the emitter?
  SF3dVector    m_MaxCreationDeviation;  //3 positive values. Declares the possible distance from the emitter
                       // Distance can be between -m_MaxCreationDeviation.? and +m_MaxCreationDeviation.?
  //Which direction are the particles emitted to?
  SF3dVector    m_StandardEmitDirection;
  SF3dVector    m_MaxEmitDirectionDeviation; //Works like m_MaxCreationDeviation

  //Which speed do they have when they are emitted?
  //->Somewhere between these speeds:
  float     m_fMinEmitSpeed;
  float     m_fMaxEmitSpeed;

  //How fast do they spin when being emitted? Speed here is angle speed (radian measure) per sec
  float     m_fMinEmitSpinSpeed;
  float     m_fMaxEmitSpinSpeed;
  //Spinning acceleration:
  float     m_fMinSpinAcceleration;
  float     m_fMaxSpinAcceleration;

  //The acceleration vector always has the same direction (normally (0/-1/0) for gravity):
  SF3dVector    m_AccelerationDirection;
  //...but not the same amount:
  float     m_fMinAcceleration;
  float     m_fMaxAcceleration;
  
  //How translucent are the particles when they are created?
  float     m_fMinEmitAlpha;
  float     m_fMaxEmitAlpha;
  //How translucent are the particles when they have reached their dying age?
  float     m_fMinDieAlpha;
  float     m_fMaxDieAlpha;

  //How big are the particles when they are created / when they die
  float     m_fMinEmitSize;
  float     m_fMaxEmitSize;
  float     m_fMinDieSize;
  float     m_fMaxDieSize;


  //The same with the color:
  SF3dVector    m_MinEmitColor;
  SF3dVector    m_MaxEmitColor;
  SF3dVector    m_MinDieColor;
  SF3dVector    m_MaxDieColor;

//*************************************
// OTHER PARTICLE INFORMATION
//*************************************

  //How long shall the particles live? Somewhere (randomly) between:
  float     m_fMinDieAge;
  float     m_fMaxDieAge;

  bool      m_bRecreateWhenDied;  //Set it true so a particle will be recreate itsself as soon
                      //as it died

//*************************************
// RENDERING PROPERTIES
//************************************* 

  int       m_iBillboarding;    //See the constants above
  
  //COGLTexture * m_Texture;        //Pointer to the texture (which is only an "alpha texture") 
  bool      m_bUseTexture;    //Set it false if you want to use GL_POINTS as particles!

  bool      m_bParticlesLeaveSystem;  //Switch it off if the particle's positions 
                       //shall be relative to the system's position (emitter position)
    
//*************************************
// STORING THE PARTICLES
//************************************* 
  //Particle array:
  CCCParticle    *m_pParticles;
  //Maximum number of particles (assigned when reserving mem for the particle array)
  int       m_iMaxParticles;
  //How many particles are currently in use?
  int       m_iParticlesInUse;
  //How many particles are created per second?
  //Note that this is an average value and if you set it too high, there won't be
  //dead particles that can be created unless the lifetime is very short and/or 
  //the array of particles (m_pParticles) is big 
  int       m_iParticlesCreatedPerSec;  //if bRecreateWhenDied is true, this is the ADDITIONAL number of created particles!
  float     m_fCreationVariance; //Set it 0 if the number of particles created per sec 
                     //should be the same each second. Otherwise use a positive value:
                     //Example: 1.0 affects that the NumParticlesCreatedPerSec varies 
                     //between m_iParticlesCreatedPerSec/2 and 1.5*m_iParticlesCreatedPerSec
//Do not set these values:
  float     m_fCurrentPointSize;  //required when rendering without particles
  //If Billboarding is set to NONE, the following vectors are (1,0,0) and (0,1,0).
  //If it is switched on, they are modified according to the viewdir/camera position (in Render of the System)
  SF3dVector    m_BillboardedX;
  SF3dVector    m_BillboardedY;     


//*************************************
// FUNCTIONS TO ASSIGN THESE MANY VECTORS MORE EASILY
//*************************************

  //Set the emitter position (you can pass a vector or x,y and z)
  void      SetEmitter(float x, float y, float z, float EmitterDeviationX,float EmitterDeviationY,float EmitterDeviationZ);
  void      SetEmitter(SF3dVector pos,SF3dVector dev);
  
  //Set the emission direction:
  void      SetEmissionDirection(float x, float y, float z,         //direction
                     float MaxDeviationX, float MaxDeviationY, float MaxDeviationZ);  //max deviation
  void      SetEmissionDirection(SF3dVector direction, SF3dVector Deviation);

  //Spin Speed
  void      SetSpinSpeed(float min, float max);
  
  //Acceleration
  void      SetAcceleration(float x, float y, float z, float Min, float Max);
  void      SetAcceleration(SF3dVector acc, float Min, float Max);

  //Color (at creation and dying age):
  void      SetCreationColor(float minr, float ming, float minb,
                     float maxr, float maxg, float maxb);
  void      SetCreationColor(SF3dVector min, SF3dVector max);

  void      SetDieColor   (float minr, float ming, float minb,
                     float maxr, float maxg, float maxb);
  void      SetDieColor   (SF3dVector min, SF3dVector max);
  //alpha:
  void      SetAlphaValues (float MinEmit, float MaxEmit, float MinDie, float MaxDie);
  //size:
  void      SetSizeValues (float EmitMin, float EmitMax, float DieMin, float DieMax);
    
//*************************************
// FUNCTIONS TO INITIALIZE THE SYSTEM
//*************************************

  CCCParticleSystem();                //constructor: sets default values

  bool      Initialize(int iNumParticles);    //reserves space for the particles

  bool      LoadTextureFromFile();

//*************************************
// FUNCTIONS TO UPDATE/RENDER THE SYSTEM
//*************************************

  void      UpdateSystem(float timePassed); //updates all particles alive
  void      Render();             //renders all particles alive

};

void CCCParticle::Initialize(CCCParticleSystem *ParentSystem)
{

  //Calculate the age, the particle will live:
  m_fDieAge = ParentSystem->m_fMinDieAge + 
           ((ParentSystem->m_fMaxDieAge - ParentSystem->m_fMinDieAge)*RANDOM_FLOAT);
  if (m_fDieAge == 0.0f) return;  //make sure there is no div 0
  m_fAge = 0.0f;

  //set the position:
  if (ParentSystem->m_bParticlesLeaveSystem)
  {
    //start with "global" coordinates (the current coordinates of the emitter position)
    m_Position = ParentSystem->m_EmitterPosition;
  }
  else
  {
    //In this case we assume a local coordinate system:
    m_Position = NULL_VECTOR;
  }
  //Add the deviation from the emitter position:
  m_Position.x += ParentSystem->m_MaxCreationDeviation.x*(RANDOM_FLOAT*2.0f-1.0f);
  m_Position.y += ParentSystem->m_MaxCreationDeviation.y*(RANDOM_FLOAT*2.0f-1.0f);
  m_Position.z += ParentSystem->m_MaxCreationDeviation.z*(RANDOM_FLOAT*2.0f-1.0f);
  //set the emission velocity
  m_Velocity.x = ParentSystem->m_StandardEmitDirection.x + ParentSystem->m_MaxEmitDirectionDeviation.x*(RANDOM_FLOAT*2.0f-1.0f);
  m_Velocity.y = ParentSystem->m_StandardEmitDirection.y + ParentSystem->m_MaxEmitDirectionDeviation.y*(RANDOM_FLOAT*2.0f-1.0f);
  m_Velocity.z = ParentSystem->m_StandardEmitDirection.z + ParentSystem->m_MaxEmitDirectionDeviation.z*(RANDOM_FLOAT*2.0f-1.0f);
  m_Velocity = m_Velocity*((ParentSystem->m_fMinEmitSpeed + 
                             (ParentSystem->m_fMaxEmitSpeed - ParentSystem->m_fMinEmitSpeed)*RANDOM_FLOAT));
  //m_Velocity = operator*(m_Velocity,(ParentSystem->m_fMinEmitSpeed + (ParentSystem->m_fMaxEmitSpeed - ParentSystem->m_fMinEmitSpeed)*RANDOM_FLOAT));
  //set the acceleration vector:
  m_Acceleration = ParentSystem->m_AccelerationDirection* 
                  (ParentSystem->m_fMinAcceleration + (ParentSystem->m_fMaxAcceleration-ParentSystem->m_fMinAcceleration)*RANDOM_FLOAT);
  //set the alpha / color values:
  m_Color = ParentSystem->m_MinEmitColor + 
       ((ParentSystem->m_MaxEmitColor-ParentSystem->m_MinEmitColor) * RANDOM_FLOAT);
  //calculate the "end color" (in order to get the ColorChange):
  SF3dVector EndColor = ParentSystem->m_MinDieColor + 
       ((ParentSystem->m_MaxDieColor-ParentSystem->m_MinDieColor) * RANDOM_FLOAT);
  m_ColorChange = (EndColor-m_Color) / m_fDieAge;

  m_fAlpha = ParentSystem->m_fMinEmitAlpha 
           + ((ParentSystem->m_fMaxEmitAlpha - ParentSystem->m_fMinEmitAlpha) * RANDOM_FLOAT);
  float fEndAlpha = ParentSystem->m_fMinDieAlpha 
           + ((ParentSystem->m_fMaxDieAlpha - ParentSystem->m_fMinDieAlpha) * RANDOM_FLOAT);
  m_fAlphaChange = (fEndAlpha - m_fAlpha) / m_fDieAge;

  //set the size values:
  m_fSize = ParentSystem->m_fMinEmitSize +
       ((ParentSystem->m_fMaxEmitSize - ParentSystem->m_fMinEmitSize) * RANDOM_FLOAT);
  float fEndSize = ParentSystem->m_fMinDieSize +
       ((ParentSystem->m_fMaxDieSize - ParentSystem->m_fMinDieSize) * RANDOM_FLOAT);
  m_fSizeChange = (fEndSize - m_fSize) / m_fDieAge;

  //spin values:
  m_fSpinAngle = 0.0f;
  m_fSpinSpeed = ParentSystem->m_fMinEmitSpinSpeed +
      ((ParentSystem->m_fMaxEmitSpinSpeed - ParentSystem->m_fMinEmitSpinSpeed) * RANDOM_FLOAT);
  m_fSpinAcceleration = ParentSystem->m_fMinSpinAcceleration +
      ((ParentSystem->m_fMaxSpinAcceleration - ParentSystem->m_fMinSpinAcceleration) * RANDOM_FLOAT);



  //Ok, we're done:
  m_bIsAlive = true;
  m_ParentSystem = ParentSystem;


}

void CCCParticle::Update(float timePassed)
{
  //Update all time-dependent values:
  m_fAge += timePassed;
  if (m_fAge >= m_fDieAge) 
  {
    if (m_ParentSystem->m_bRecreateWhenDied) 
    {
      Initialize(m_ParentSystem);
      Update(RANDOM_FLOAT * timePassed);  //see comment in UpdateSystem
    }
    else
    {
      m_fAge = 0.0f;
      m_bIsAlive = false;
      m_ParentSystem->m_iParticlesInUse--;
    }

    return;
  }

  m_fSize  += m_fSizeChange *timePassed;
  m_fAlpha += m_fAlphaChange*timePassed;
  m_Color = m_Color + m_ColorChange*timePassed;
  m_Velocity = m_Velocity + m_Acceleration*timePassed;
  //Note: exact would be: m_Position = 1/2*m_Acceleration*timePassed² + m_VelocityOLD*timePassed;
  //But this approach is ok, I think!
  m_Position = m_Position + (m_Velocity*timePassed);

  m_fSpinSpeed += m_fSpinAcceleration*timePassed;
  m_fSpinAngle += m_fSpinSpeed*timePassed;

  //That's all!
}

void CCCParticle::Render()
{
  if (!m_ParentSystem->m_bUseTexture) 
  {
    glPointSize(m_fSize*m_ParentSystem->m_fCurrentPointSize);
    float color[4];
    color[0] = m_Color.x;
    color[1] = m_Color.y;
    color[2] = m_Color.z;
    color[3] = m_fAlpha;

    glColor4fv(&color[0]);

    glBegin(GL_POINTS);
      glVertex3fv(&m_Position.x);
    glEnd();
  }
  else
  {
    //render using texture: (texture was already set active by the Render method of the particle system)
    
    float color[4];
    color[0] = m_Color.x;
    color[1] = m_Color.y;
    color[2] = m_Color.z;
    color[3] = m_fAlpha;
    glColor4fv(&color[0]);
    gluBuild2DMipmaps( GL_TEXTURE_2D, 4,5,5,GL_RGBA,GL_UNSIGNED_BYTE,color);

    SF3dVector RotatedX = m_ParentSystem->m_BillboardedX;
    SF3dVector RotatedY = m_ParentSystem->m_BillboardedY;


     //If spinning is switched on, rotate the particle now:
    if (m_fSpinAngle > 0.0f)
    {
      RotatedX = m_ParentSystem->m_BillboardedX * cos(m_fSpinAngle) 
               + m_ParentSystem->m_BillboardedY * sin(m_fSpinAngle);
      RotatedY = m_ParentSystem->m_BillboardedY * cos(m_fSpinAngle) 
               - m_ParentSystem->m_BillboardedX * sin(m_fSpinAngle);
    }
  
    
    //Render a quadrangle with the size m_fSize
    SF3dVector coords = m_Position - (RotatedX*(0.5f*m_fSize))
                     - (RotatedY*(0.5f*m_fSize));
    glBegin(GL_POLYGON);
      glVertex3fv(&coords.x);
      glTexCoord2f(0.0f,1.0f);
      coords = coords + RotatedY * m_fSize;
      glVertex3fv(&coords.x);
      glTexCoord2f(1.0f,1.0f);
      coords = coords + RotatedX * m_fSize;     
      glVertex3fv(&coords.x);
      glTexCoord2f(1.0f,0.0f);
      coords = coords - RotatedY * m_fSize;
      glVertex3fv(&coords.x);
      glTexCoord2f(0.0f,0.0f);
    glEnd();


  }

}

/*************************************

METHODS OF CCCParticleSytem class

**************************************/


CCCParticleSystem::CCCParticleSystem()
{
  //Set default values:
  
  //motion:
  this->m_EmitterPosition = NULL_VECTOR;
  this->m_MaxCreationDeviation = NULL_VECTOR;

  this->m_StandardEmitDirection = NULL_VECTOR;
  this->m_MaxEmitDirectionDeviation = NULL_VECTOR;
  this->m_fMaxEmitSpeed = 0.0f;
  this->m_fMinEmitSpeed = 0.0f;

  this->m_AccelerationDirection = NULL_VECTOR;
  this->m_fMaxAcceleration = 0.0f;
  this->m_fMinAcceleration = 0.0f;

  this->m_fMinEmitSpinSpeed = 0.0f;
  this->m_fMaxEmitSpinSpeed = 0.0f;
  
  this->m_fMaxSpinAcceleration = 0.0f;
  this->m_fMinSpinAcceleration = 0.0f;


  //look:
  this->m_fMaxEmitAlpha = 0.0f;
  this->m_fMinEmitAlpha = 0.0f;
  this->m_fMaxDieAlpha = 1.0f;
  this->m_fMinDieAlpha = 1.0f;

  this->m_MaxEmitColor = NULL_VECTOR;
  this->m_MinEmitColor = NULL_VECTOR;
  this->m_MaxDieColor = NULL_VECTOR;
  this->m_MinDieColor = NULL_VECTOR;

  //this->m_Texture = NULL;
  this->m_bUseTexture = false;
  this->m_iBillboarding = BILLBOARDING_NONE;

  //size:
  this->m_fMaxEmitSize = 0.0f;
  this->m_fMinEmitSize = 0.0f;
  this->m_fMaxDieSize = 0.0f;
  this->m_fMinDieSize = 0.0f;


  //behavior:
  this->m_bRecreateWhenDied = false;
  
  this->m_fMaxDieAge = 1.0f;
  this->m_fMinDieAge = 1.0f;

  this->m_iMaxParticles = 0;  //array is not yet created
  this->m_iParticlesInUse = 0;

  this->m_iParticlesCreatedPerSec = 0;
  this->m_fCreationVariance = 0.0f;
  this->m_bParticlesLeaveSystem = false;
  this->m_pParticles = NULL;

}
//*********************************************************
void CCCParticleSystem::SetEmitter(float x, float y, float z, float EmitterDeviationX,float EmitterDeviationY,float EmitterDeviationZ)
{
  SetEmitter(F3dVector(x,y,z),F3dVector(EmitterDeviationX,EmitterDeviationY,EmitterDeviationZ));
}

void CCCParticleSystem::SetEmitter(SF3dVector pos, SF3dVector dev)
{
  m_EmitterPosition = pos;
  m_MaxCreationDeviation = dev;
}

void CCCParticleSystem::SetEmissionDirection(float x, float y, float z,
                       float MaxDeviationX, float MaxDeviationY, float MaxDeviationZ)
{
  SetEmissionDirection(F3dVector(x,y,z),F3dVector(MaxDeviationX,MaxDeviationY,MaxDeviationZ));
}


void CCCParticleSystem::SetEmissionDirection(SF3dVector direction, SF3dVector Deviation)
{
  m_StandardEmitDirection = direction;
  m_MaxEmitDirectionDeviation = Deviation;
}




void CCCParticleSystem::SetSpinSpeed(float min, float max)
{
  m_fMinEmitSpinSpeed = min;
  m_fMaxEmitSpinSpeed = max;
}


void CCCParticleSystem::SetAcceleration(float x, float y, float z, float Min, float Max)
{
  SetAcceleration(F3dVector(x,y,z),Min,Max);
}

void CCCParticleSystem::SetAcceleration(SF3dVector acc, float Min, float Max)
{
  m_AccelerationDirection = acc;
  m_fMaxAcceleration = Max;
  m_fMinAcceleration = Min;
}

void CCCParticleSystem::SetCreationColor(float minr, float ming, float minb,
                       float maxr, float maxg, float maxb)
{
  SetCreationColor(F3dVector(minr,ming,minb),F3dVector(maxr,maxg,maxb));
}

void CCCParticleSystem::SetCreationColor(SF3dVector min, SF3dVector max)
{
  m_MinEmitColor = min;
  m_MaxEmitColor = max;
}


void CCCParticleSystem::SetDieColor (float minr, float ming, float minb,
                       float maxr, float maxg, float maxb)
{
  SetDieColor(F3dVector(minr,ming,minb),F3dVector(maxr,maxg,maxb));
}

void CCCParticleSystem::SetDieColor   (SF3dVector min, SF3dVector max)
{
  m_MinDieColor = min;
  m_MaxDieColor = max;
}

void CCCParticleSystem::SetAlphaValues (float MinEmit, float MaxEmit, float MinDie, float MaxDie)
{
  m_fMinEmitAlpha = MinEmit;
  m_fMaxEmitAlpha = MaxEmit;
  m_fMinDieAlpha = MinDie;
  m_fMaxDieAlpha = MaxDie;
}

void CCCParticleSystem::SetSizeValues (float EmitMin, float EmitMax, float DieMin, float DieMax)
{
  m_fMinEmitSize = EmitMin;
  m_fMaxEmitSize = EmitMax;
  m_fMinDieSize = DieMin;
  m_fMaxDieSize = DieMax;
}
//*********************************************************

bool CCCParticleSystem::Initialize(int iNumParticles)
{
  this->m_pParticles = new CCCParticle[iNumParticles];
  if (m_pParticles == NULL) 
  {
    return false;
    this->m_iMaxParticles = 0;
    this->m_iParticlesInUse = 0;
  }

  this->m_iMaxParticles = iNumParticles;
  this->m_iParticlesInUse = 0;

  //Set the status of each particle to DEAD
  for (int i = 0; i < iNumParticles; i++)
  {
    m_pParticles[i].m_bIsAlive = false;
  }

  return true;


  
}


void CCCParticleSystem::UpdateSystem(float timePassed)
{
  //We have to 
  //  -update the particles (= move the particles, change their alpha, color, speed values)
  //  -create new particles, if desired and there are "free" particles

  //First get the number of particles we want to create (randomly in a certain dimension (dependent of m_CreationVariance)
  
  int iParticlesToCreate = (int) ((float)m_iParticlesCreatedPerSec
                       *timePassed
                                       *(1.0f+m_fCreationVariance*(RANDOM_FLOAT-0.5f)));
  

  //loop through the particles and update / create them
  for (int i = 0; i < m_iMaxParticles; i++)
  {
    if (m_pParticles[i].m_bIsAlive)
    {
      m_pParticles[i].Update(timePassed);
    }

    //Should we create the particle?
    if (iParticlesToCreate > 0)
    {
      if (!m_pParticles[i].m_bIsAlive)
      {
        m_pParticles[i].Initialize(this);
        //Update the particle: This has an effect, as if the particle would have
        //been emitted some milliseconds ago. This is very useful on slow PCs:
        //Especially if you simulate something like rain, then you could see that 
        //many particles are emitted at the same time (same "UpdateSystem" call),
        //if you would not call this function:        
        m_pParticles[i].Update(RANDOM_FLOAT*timePassed);  
        iParticlesToCreate--;
      }
    }

  }
  
}

bool CCCParticleSystem::LoadTextureFromFile()
{ 
  glBindTexture( GL_TEXTURE_2D, fire);

  m_bUseTexture = true;


  return true;

}


void CCCParticleSystem::Render()
{
  //the calling method must switch on texturing!

  if (m_bUseTexture)
  {
    glBindTexture( GL_TEXTURE_2D, fire);
   // m_Texture->SetActive();
    //Calculate the "billboarding vectors" (the particles only store their positions, but we need quadrangles!)
    switch (m_iBillboarding)
    {
    case BILLBOARDING_NONE:
      {
        //independent from camera / view direction
        m_BillboardedX = F3dVector(1.0f,0.0f,0.0f);
        m_BillboardedY = F3dVector(0.0f,1.0f,0.0f);
        break;
      }
    case BILLBOARDING_PERPTOVIEWDIR:
      {
        //Retrieve the up and right vector from the modelview matrix:
        float fModelviewMatrix[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, fModelviewMatrix);

        //Assign the x-Vector for billboarding:
        m_BillboardedX = F3dVector(fModelviewMatrix[0], fModelviewMatrix[4], fModelviewMatrix[8]);

        //Assign the y-Vector for billboarding:
        m_BillboardedY = F3dVector(fModelviewMatrix[1], fModelviewMatrix[5], fModelviewMatrix[9]);
        break;
      }
    case BILLBOARDING_PERPTOVIEWDIR_BUTVERTICAL:
      {
        //Retrieve the right vector from the modelview matrix:
        float fModelviewMatrix[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, fModelviewMatrix);

        //Assign the x-Vector for billboarding:
        m_BillboardedX = F3dVector(fModelviewMatrix[0], fModelviewMatrix[4], fModelviewMatrix[8]);

        //Assign the y-Vector:
        m_BillboardedY = F3dVector(0.0f,1.0f,0.0f);       
        break;
      }
    }
  }
  else
  {
    glGetFloatv(GL_POINT_SIZE,&m_fCurrentPointSize);
  }
  for (int i = 0; i < m_iMaxParticles; i++)
  {
    if (m_pParticles[i].m_bIsAlive)
      m_pParticles[i].Render();
  }
}

CCCParticleSystem g_ParticleSystem1;
CCCParticleSystem g_ParticleSystem2;
CCCParticleSystem g_ParticleSystem4;
clock_t g_iLastRenderTime;
clock_t g_iLastRenderTime1;

void InitParticles()
{

//INIT SYSTEM 1 (POINTS, FIREWORK)
  g_ParticleSystem1.Initialize(800);  //particle system must not have more than 800 particles
  g_ParticleSystem1.m_iParticlesCreatedPerSec = 800;  //we create all particles in the first second of the system's life
  g_ParticleSystem1.m_fMinDieAge = 2.5f;      //but the particles live longer than one second
  g_ParticleSystem1.m_fMaxDieAge = 2.5f;      //-> this causes the system to "shoot" periodically
    
  g_ParticleSystem1.m_fCreationVariance = 1.0f;
  g_ParticleSystem1.m_bRecreateWhenDied = true;
  g_ParticleSystem1.SetCreationColor(1.0f,1.0f,1.0f,
                0.5f,0.5f,0.5f);
  g_ParticleSystem1.SetDieColor(0.0f,1.0f,0.0f,
                 0.0f,0.3f,0.0f);
  g_ParticleSystem1.SetAlphaValues(1.0f,1.0f,0.0f,0.0f);
  g_ParticleSystem1.SetEmitter(-10.0f,-0.0f,0.0f,
                0.02f,0.0f,0.02f);
  g_ParticleSystem1.SetAcceleration(F3dVector(0.0f,-1.0f,0.0f),0.83f,1.4f);
  g_ParticleSystem1.SetSizeValues(3.0f,3.0f,4.0f,4.0f);
  g_ParticleSystem1.m_fMaxEmitSpeed = 0.82f;
  g_ParticleSystem1.m_fMinEmitSpeed = 1.3f;
  g_ParticleSystem1.SetEmissionDirection(0.1f,2.0f,0.1f,
                    0.5f,0.5f,0.5f);

  g_ParticleSystem1.m_bParticlesLeaveSystem = true;


//INIT SYSTEM 2 (FIRE)

  g_ParticleSystem2.Initialize(300);
  g_ParticleSystem2.m_iParticlesCreatedPerSec = 300;
  g_ParticleSystem2.m_fCreationVariance = 0.0f;
  g_ParticleSystem2.m_bRecreateWhenDied = true;
  g_ParticleSystem2.m_fMinDieAge = 0.5f;
  g_ParticleSystem2.m_fMaxDieAge = 1.5f;
  g_ParticleSystem2.SetCreationColor(1.0f,0.0f,0.0f,
                  1.0f,0.5f,0.0f);
  g_ParticleSystem2.SetDieColor(1.0f,0.0f,0.0f,
                    1.0f,0.5f,0.0f);

  g_ParticleSystem2.SetAlphaValues(1.0f,1.0f,0.0f,0.0f);
  g_ParticleSystem2.SetEmitter(0.0f,0.0f,0.5f,
                0.3f,0.0f,0.3f);
  g_ParticleSystem2.SetAcceleration(F3dVector(0.0f,1.0f,0.0f),0.3f,0.4f);
  g_ParticleSystem2.SetSizeValues(0.04f,0.08f,0.06f,0.12f);
  g_ParticleSystem2.m_fMaxEmitSpeed = 0.1f;
  g_ParticleSystem2.m_fMinEmitSpeed = 0.2f;
  g_ParticleSystem2.SetEmissionDirection(0.0f,1.0f,0.0f,
                    0.08f,0.5f,0.08f);
  g_ParticleSystem2.m_bParticlesLeaveSystem = true;
  g_ParticleSystem2.SetSpinSpeed(-0.82*PI,0.82*PI);
  g_ParticleSystem2.m_iBillboarding = BILLBOARDING_PERPTOVIEWDIR;
  g_ParticleSystem2.LoadTextureFromFile();
}



//water simulation--draw the pool 


void CreatePool()
{
  
  NumOscillators = NUM_OSCILLATORS;
  Oscillators = new SOscillator[NumOscillators];
  IndexVect.clear();  //to be sure it is empty
  for (int xc = 0; xc < NUM_X_OSCILLATORS; xc++) 
    for (int zc = 0; zc < NUM_Z_OSCILLATORS; zc++) 
    {
      Oscillators[xc+zc*NUM_X_OSCILLATORS].x = OSCILLATOR_DISTANCE*float(xc);
      Oscillators[xc+zc*NUM_X_OSCILLATORS].y = 0.0f;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].z = OSCILLATOR_DISTANCE*float(zc);

      Oscillators[xc+zc*NUM_X_OSCILLATORS].nx = 0.0f;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].ny = 1.0f;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].nz = 0.0f;

      Oscillators[xc+zc*NUM_X_OSCILLATORS].UpSpeed = 0;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].bIsExciter = false;

      //create two triangles:
      if ((xc < NUM_X_OSCILLATORS-1) && (zc < NUM_Z_OSCILLATORS-1))
      {
        IndexVect.push_back(xc+zc*NUM_X_OSCILLATORS);
        IndexVect.push_back((xc+1)+zc*NUM_X_OSCILLATORS);
        IndexVect.push_back((xc+1)+(zc+1)*NUM_X_OSCILLATORS);

        IndexVect.push_back(xc+zc*NUM_X_OSCILLATORS);
        IndexVect.push_back((xc+1)+(zc+1)*NUM_X_OSCILLATORS);
        IndexVect.push_back(xc+(zc+1)*NUM_X_OSCILLATORS);
      }

    }

  //copy the indices:
  Indices = new GLuint[IndexVect.size()];  //allocate the required memory
  for (int i = 0; i < IndexVect.size(); i++)
  {
    Indices[i] = IndexVect[i];
  }

  Oscillators[100+30*NUM_X_OSCILLATORS].bIsExciter = true;
  Oscillators[100+30*NUM_X_OSCILLATORS].ExciterAmplitude = 0.5f;
  Oscillators[100+30*NUM_X_OSCILLATORS].ExciterFrequency = 50.0f;
  Oscillators[30+80*NUM_X_OSCILLATORS].bIsExciter = true;
  Oscillators[30+80*NUM_X_OSCILLATORS].ExciterAmplitude = 0.5f;
  Oscillators[30+80*NUM_X_OSCILLATORS].ExciterFrequency = 50.0f;
  NumIndices = IndexVect.size();
  IndexVect.clear();  //no longer needed, takes only memory
}

//water simulation--physicial simulation
void UpdateScene(bool bEndIsFree, float deltaTime, float time)
{
  for (int xc = 0; xc < NUM_X_OSCILLATORS; xc++) 
  {
    for (int zc = 0; zc < NUM_Z_OSCILLATORS; zc++) 
    {
      int ArrayPos = xc+zc*NUM_X_OSCILLATORS;

      //check, if oscillator is an exciter (these are not affected by other oscillators)
      if ((Oscillators[ArrayPos].bIsExciter) && g_bExcitersInUse)
      {
        Oscillators[ArrayPos].newY = Oscillators[ArrayPos].ExciterAmplitude*sin(time*Oscillators[ArrayPos].ExciterFrequency);
      }


      //check, if this oscillator is on an edge (=>end)
      if ((xc==0) || (xc==NUM_X_OSCILLATORS-1) || (zc==0) || (zc==NUM_Z_OSCILLATORS-1))
        ;//TBD: calculating oscillators at the edge (if the end is free)
      else
      {
        //calculate the new speed:
        

        //Change the speed (=accelerate) according to the oscillator's 4 direct neighbors:
        float AvgDifference = Oscillators[ArrayPos-1].y       //left neighbor
                   +Oscillators[ArrayPos+1].y       //right neighbor
                   +Oscillators[ArrayPos-NUM_X_OSCILLATORS].y  //upper neighbor
                   +Oscillators[ArrayPos+NUM_X_OSCILLATORS].y  //lower neighbor
                   -4*Oscillators[ArrayPos].y;        //subtract the pos of the current osc. 4 times  
        Oscillators[ArrayPos].UpSpeed += AvgDifference*deltaTime/OSCILLATOR_WEIGHT;

        //calculate the new position, but do not yet store it in "y" (this would affect the calculation of the other osc.s)
        Oscillators[ArrayPos].newY += Oscillators[ArrayPos].UpSpeed*deltaTime;
        
        
        
      }
    }   
  }

  //copy the new position to y:
  for ( int xc = 0; xc < NUM_X_OSCILLATORS; xc++) 
  {
    for (int zc = 0; zc < NUM_Z_OSCILLATORS; zc++) 
    {
      Oscillators[xc+zc*NUM_X_OSCILLATORS].y =Oscillators[xc+zc*NUM_X_OSCILLATORS].newY;
    }
  }
  //calculate new normal vectors (according to the oscillator's neighbors):
  for ( int xc = 0; xc < NUM_X_OSCILLATORS; xc++) 
  {
    for (int zc = 0; zc < NUM_Z_OSCILLATORS; zc++) 
    {
      ///
      //Calculating the normal:
    

      SF3dVector u,v,p1,p2; //u and v are direction vectors. p1 / p2: temporary used (storing the points)

      if (xc > 0) p1 = F3dVector(Oscillators[xc-1+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc-1+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc-1+zc*NUM_X_OSCILLATORS].z);
      else
            p1 = F3dVector(Oscillators[xc+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].z); 
      if (xc < NUM_X_OSCILLATORS-1) 
            p2 = F3dVector(Oscillators[xc+1+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+1+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+1+zc*NUM_X_OSCILLATORS].z);
      else
            p2 = F3dVector(Oscillators[xc+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].z); 
      u = p2-p1; //vector from the left neighbor to the right neighbor
      if (zc > 0) p1 = F3dVector(Oscillators[xc+(zc-1)*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+(zc-1)*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+(zc-1)*NUM_X_OSCILLATORS].z);
      else
            p1 = F3dVector(Oscillators[xc+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].z); 
      if (zc < NUM_Z_OSCILLATORS-1) 
            p2 = F3dVector(Oscillators[xc+(zc+1)*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+(zc+1)*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+(zc+1)*NUM_X_OSCILLATORS].z);
      else
            p2 = F3dVector(Oscillators[xc+zc*NUM_X_OSCILLATORS].x,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].y,
                     Oscillators[xc+zc*NUM_X_OSCILLATORS].z); 
      v = p2-p1; //vector from the upper neighbor to the lower neighbor
      //calculat the normal:
      SF3dVector normal = Normalize3dVector(CrossProduct(&u,&v));

      //assign the normal:
      Oscillators[xc+zc*NUM_X_OSCILLATORS].nx = normal.x;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].ny = normal.y;
      Oscillators[xc+zc*NUM_X_OSCILLATORS].nz = normal.z;
    }
  }

}
//draw the scene with pool
void DrawScene(void)
{
  float white[] = {1,1,1,1};
   float blue[] = {0.2,0.6,1,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,blue);
  glDrawElements( GL_TRIANGLES, //mode
            NumIndices,  //count, ie. how many indices
            GL_UNSIGNED_INT, //type of the index array
            Indices);;

}

//first person navigation
void FirstpersonNaviagtion(void)
{
      gluLookAt(Ex, Ey, Ez, Ex + 0.2*scale*10*Sin(zh), Ey+0.2*scale*50*Sin(theta), Ez - 0.2*scale*10*Cos(zh),0.0, 1.0, 0.0);
}




/* 
 *  Draw day sky box
 */
static void Sky(double D)
{
   glColor3f(1,1,1);
   glEnable(GL_TEXTURE_2D);

   //  Sides
   glBindTexture(GL_TEXTURE_2D,sky[0]);
   glBegin(GL_QUADS);
   glTexCoord2f(0.00,0.45); glVertex3f(-D,-0.7,-D);
   glTexCoord2f(0.25,0.45); glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.25,0.62); glVertex3f(+D,+D,-D);
   glTexCoord2f(0.00,0.62); glVertex3f(-D,+D,-D);

   glTexCoord2f(0.25,0.45); glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.50,0.45); glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.50,0.62); glVertex3f(+D,+D,+D);
   glTexCoord2f(0.25,0.62); glVertex3f(+D,+D,-D);

   glTexCoord2f(0.50,0.45); glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.75,0.45); glVertex3f(-D,-0.7,+D);
   glTexCoord2f(0.75,0.62); glVertex3f(-D,+D,+D);
   glTexCoord2f(0.50,0.62); glVertex3f(+D,+D,+D);

   glTexCoord2f(0.75,0.45); glVertex3f(-D,-0.7,+D);
   glTexCoord2f(1.00,0.45); glVertex3f(-D,-0.7,-D);
   glTexCoord2f(1.00,0.62); glVertex3f(-D,+D,-D);
   glTexCoord2f(0.75,0.62); glVertex3f(-D,+D,+D);
   glEnd();

   //  Top and bottom
   //glBindTexture(GL_TEXTURE_2D,sky[1]);
   glBegin(GL_QUADS);
   glTexCoord2f(0.25,0.62);glVertex3f(+D,+D,-D);
   glTexCoord2f(0.5,0.62);glVertex3f(+D,+D,+D);
   glTexCoord2f(0.5,0.87);glVertex3f(-D,+D,+D);
   glTexCoord2f(0.25,0.87);glVertex3f(-D,+D,-D);

   glTexCoord2f(0.5,0.12);glVertex3f(-D,-0.7,+D);
   glTexCoord2f(0.5,0.37);glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.25,0.12);glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.25,0.37); glVertex3f(-D,-0.7,-D);
   glEnd();

   glDisable(GL_TEXTURE_2D);
}

/* 
 *  Draw night sky box
 */
static void Sky1(double D)
{
   glColor3f(1,1,1);
   glEnable(GL_TEXTURE_2D);

   //  Sides
   glBindTexture(GL_TEXTURE_2D,sky[1]);
   glBegin(GL_QUADS);
   glTexCoord2f(0.00,0.37); glVertex3f(-D,-0.7,-D);
   glTexCoord2f(0.255,0.37); glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.255,0.615); glVertex3f(+D,+D,-D);
   glTexCoord2f(0.00,0.615); glVertex3f(-D,+D,-D);

   glTexCoord2f(0.255,0.37); glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.495,0.37); glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.495,0.615); glVertex3f(+D,+D,+D);
   glTexCoord2f(0.255,0.615); glVertex3f(+D,+D,-D);

   glTexCoord2f(0.495,0.37); glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.75,0.37); glVertex3f(-D,-0.7,+D);
   glTexCoord2f(0.75,0.615); glVertex3f(-D,+D,+D);
   glTexCoord2f(0.495,0.615); glVertex3f(+D,+D,+D);

   glTexCoord2f(0.75,0.37); glVertex3f(-D,-0.7,+D);
   glTexCoord2f(1.00,0.37); glVertex3f(-D,-0.7,-D);
   glTexCoord2f(1.00,0.615); glVertex3f(-D,+D,-D);
   glTexCoord2f(0.75,0.615); glVertex3f(-D,+D,+D);
   glEnd();

   //  Top and bottom
   //glBindTexture(GL_TEXTURE_2D,sky[1]);
   glBegin(GL_QUADS);
   glTexCoord2f(0.255,0.615);glVertex3f(+D,+D,-D);
   glTexCoord2f(0.495,0.615);glVertex3f(+D,+D,+D);
   glTexCoord2f(0.495,0.87);glVertex3f(-D,+D,+D);
   glTexCoord2f(0.255,0.87);glVertex3f(-D,+D,-D);

   glTexCoord2f(0.495,0.125);glVertex3f(-D,-0.7,+D);
   glTexCoord2f(0.495,0.37);glVertex3f(+D,-0.7,+D);
   glTexCoord2f(0.255,0.125);glVertex3f(+D,-0.7,-D);
   glTexCoord2f(0.255,0.37); glVertex3f(-D,-0.7,-D);
   glEnd();

   glDisable(GL_TEXTURE_2D);
}

void drawCourtTower(float R, float G, float B){
    
  float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  

  //draw lower part
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(R,G,B);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   
   glTexCoord2f(0.4,0); glVertex3f(0.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(2.0,0.0,0.0);
   glTexCoord2f(0,1.5); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(0.4,1.5);glVertex3f(0.0,12.0,0.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(3.5/4,0);glVertex3f(2.0,0.0,0.0);
   glTexCoord2f(4.5/4,0);glVertex3f(2.0,0.0,2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,0.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(2.0,0.0,2.0);
   glTexCoord2f(4.5/4,0);glVertex3f(0.0,0.0,2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(0.0,12.0,2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(3.5/4,0);glVertex3f(0.0,0.0,2.0);
   glTexCoord2f(4.5/4,0);glVertex3f(0.0,0.0,0.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(0.0,12.0,0.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(0.0,12.0,2.0);
   glEnd();
   glColor3f(1.0,1.0,1.0);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   //draw higher part
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(R,G,B);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_POLYGON);
   glNormal3f(1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,12.0,2.0);
   glTexCoord2f(1.0/4,0); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(1.0/4,1.0/4); glVertex3f(2.0,14.0,0.0);
   glTexCoord2f(0.8/4,1.0/4); glVertex3f(2.0,14.0,0.4);
   //glTexCoord2f(1*1.5,0); glVertex3f(2.0,12.0,2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,14.0,0.4);
   glTexCoord2f(1.0/4,0); glVertex3f(2.0,14.0,0.0);
   glTexCoord2f(1.0/4,1.0/4); glVertex3f(2.0,15.5,0.0);
   glTexCoord2f(0,1.0/4); glVertex3f(2.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(1.0,0.0,0.0);
   
   glTexCoord2f(0.8/4,1.0/4); glVertex3f(2.0,15.5,0.4);
   glTexCoord2f(0,0); glVertex3f(2.0,16.5,2.0);
   glTexCoord2f(1.0/4,0); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(0.0,1.0/4); glVertex3f(2.0,15.5,0.0);
   //glTexCoord2f(0,0); glVertex3f(2.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,12.0,2.0);
   glTexCoord2f(1.0/4,0); glVertex3f(0.0,12.0,0.0);
   glTexCoord2f(1.0/4,1.0/4); glVertex3f(0.0,14.0,0.0);
   glTexCoord2f(0,0.8/4); glVertex3f(0.0,14.0,0.4);
   //glTexCoord2f(1*1.5,0); glVertex3f(0.0,12.0,2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,14.0,0.4);
   glTexCoord2f(1.0/4,0); glVertex3f(0.0,14.0,0.0);
   glTexCoord2f(1.0/4,1.0/4); glVertex3f(0.0,15.5,0.0);
   glTexCoord2f(0,1.0/4); glVertex3f(0.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0.8/4,1.0/4); glVertex3f(0.0,15.5,0.4);
   glTexCoord2f(0,0); glVertex3f(0.0,16.5,2.0);
   glTexCoord2f(1.0/4,1.0/4); glVertex3f(0.0,16.5,0.0);
   glTexCoord2f(0,1.0/4); glVertex3f(0.0,15.5,0.0);
  // glTexCoord2f(0,0); glVertex3f(0.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0.0,0.0,-1.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,12.0,0.0);
   glTexCoord2f(0.4,0); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(0.4,0.7); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(0,0.7); glVertex3f(0.0,16.5,0.0);
  // glTexCoord2f(0,0); glVertex3f(0.0,15.5,0.4);
   glColor3f(1.0,1.0,1.0);
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   //draw the pointer
   glPushMatrix();
   //glEnable(GL_TEXTURE_2D);
   //glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(R,G,B);
   //glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_TRIANGLES);
   //float YtoH = (1.0-0.0) / (0.0-(18.5-16.5)) * (1.0-0.0);
   glNormal3f(0.0,2.0,-4.0);
   
   //glTexCoord2f(0,0); 
   glVertex3f(0.0,16.5,0.0);
   //glTexCoord2f(0,1*1.5); 
   glVertex3f(2.0,16.5,0.0);
   //glTexCoord2f(1*1.5,1*1.5); 
   glVertex3f(1.0,18.5,1.0);
   glEnd();

   glBegin(GL_TRIANGLES);
   glNormal3f(4.0,2.0,0.0);
   
   //glTexCoord2f(0,0); 
   glVertex3f(2.0,16.5,0.0);
   //glTexCoord2f(0,1*1.5); 
   glVertex3f(2.0,16.5,2.0);
   //glTexCoord2f(1*1.5,1*1.5); 
   glVertex3f(1.0,18.5,1.0);
   glEnd();

   glBegin(GL_TRIANGLES);
   glNormal3f(0.0,2.0,4.0);
   
   //glTexCoord2f(0,0); 
   glVertex3f(2.0,16.5,2.0);
   //glTexCoord2f(0,1*1.5); 
   glVertex3f(0.0,16.5,2.0);
   //glTexCoord2f(1*1.5,1*1.5); 
   glVertex3f(1.0,18.5,1.0);
   glEnd();

  glBegin(GL_TRIANGLES);
   glNormal3f(-4.0,2.0,0.0);
   
   //glTexCoord2f(0,0); 
   glVertex3f(0.0,16.5,2.0);
   //glTexCoord2f(0,1*1.5); 
   glVertex3f(0.0,16.5,0.0);
   //glTexCoord2f(1*1.5,1*1.5); 
   glVertex3f(1.0,18.5,1.0);
   glEnd();
   glColor3f(1.0,1.0,1.0);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   //pointer
   glPushMatrix();
   //glEnable(GL_TEXTURE_2D);
   //glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(R,G,B);
   //glBindTexture(GL_TEXTURE_2D,texture[2]);
   
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:
      int i=0;
      float x,z;
      
      //y is constant when the height is same 
     float YtoLowerHeight = (0.1-0.0) / (0.0-(21-18.5)) * (0.1-0.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glVertex3f((x)*0.1+1,18.5,(z)*0.1+1);
         //same x,z and NVect: 
         glVertex3f(1.0,20.0,1.0);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z); 
      glVertex3f((x)*0.1+1,18.5,(z)*0.1+1);
      //same x,z and NVect:
      glVertex3f(1.0,20.0,1.0);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
}

/*
r = torus ring radius
c = torus tube radius
rSeg, cSeg = number of segments/detail
*/

void drawTorus(double r, double c,int rSeg, int cSeg)
{
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  glFrontFace(GL_CW);

  glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture[1]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  const double TAU = 2 * PI;

  for (int i = 0; i < rSeg; i++) {
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j <= cSeg; j++) {
      for (int k = 0; k <= 1; k++) {
        double s = (i + k) % rSeg + 0.5;
        double t = j % (cSeg + 1);

        double x = (c + r * cos(s * TAU / rSeg)) * cos(t * TAU / cSeg);
        double y = (c + r * cos(s * TAU / rSeg)) * sin(t * TAU / cSeg);
        double z = r * sin(s * TAU / rSeg);

        double u = (i + k) / (float) rSeg;
        double v = t / (float) cSeg;

        glTexCoord2d(u, v);
        glNormal3f(2 * x, 2 * y, 2 * z);
        glVertex3d(2 * x, 2 * y, 2 * z);
      }
    }
    glColor3f(1.0,1.0,1.0);
    glEnd();
  }

  glFrontFace(GL_CCW);
  glPopMatrix();
  glDisable(GL_TEXTURE_2D);
}


void drawCourtGate(){
    float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  

   glPushMatrix();
      glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.688,0.766,0.867);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glBegin(GL_QUAD_STRIP);
   int i;
   float x,z;
      for ( i = 0; i < NumOfEdges; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges*360);
         z = Sin((float)i/(float)NumOfEdges*360);
         glNormal3f(x,0.0,z);
         glTexCoord2f((float)i/(float)NumOfEdges*2,LowerHeight/4); glVertex3f(x*0.1,0.0,z*0.1);
         glTexCoord2f((float)i/(float)NumOfEdges*2,HigherHeight/4); glVertex3f(x*0.1,HigherHeight,z*0.1);
      }
      x = Cos((float)i/(float)NumOfEdges*360);
      z = Sin((float)i/(float)NumOfEdges*360);
      glNormal3f(x,0.0,z);
      glTexCoord2f((float)i/(float)NumOfEdges*2,LowerHeight/4); glVertex3f(x*0.1,0.0,z*0.1);
      glTexCoord2f((float)i/(float)NumOfEdges*2,HigherHeight/4); glVertex3f(x*0.1,HigherHeight,z*0.1);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      glPushMatrix();
      glTranslatef(0.0,HigherHeight+0.2,0.0);
      drawTorus(0.05,0.2,16,8);
      glPopMatrix();

}
void drawCourt(){


  float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  
 glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.0,1.0,0.0);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glScalef(0.5,0.0,1.0);
   //draw the ground

      //glBegin(GL_TRIANGLES);
      int i;
      float x,z;
     
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

     glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glScalef(0.5,1.0,1.0);
   //draw the court wall
      glBegin(GL_QUADS);
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glTexCoord2f(0.0,0.0);glVertex3f(x*20,0.0,z*20);
         glTexCoord2f(0.0,1.0/4);glVertex3f(x*20,5.0,z*20);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glTexCoord2f(1.0/4,1.0/4);glVertex3f(x*20,5.0,z*20);
         glTexCoord2f(1.0/4,0.0);glVertex3f(x*20,0.0,z*20);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glTexCoord2f(0.0,0.0);glVertex3f(x*20,0.0,z*20);
      glTexCoord2f(0.0,1.0/4);glVertex3f(x*20,5.0,z*20);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glTexCoord2f(1.0/4,1.0/4);glVertex3f(x*20,5.0,z*20);
      glTexCoord2f(1.0/4,0.0);glVertex3f(x*20,0.0,z*20);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
      //draw court tower
      glPushMatrix();
      glTranslatef(10.0,0.0,0.0);
      glRotatef(270,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,0.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-10.0,0.0,1.5);
      glRotatef(90,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,0.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(1.5,0.0,20.0);
      glRotatef(180,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,1.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(0.0,0.0,-20.0);
      glRotatef(0,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,1.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(7.0,0.0,15.0);
      glRotatef(210,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,0.0,1.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-7.0,0.0,-14.5);
      glRotatef(20,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,0.0,1.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(6.3,0.0,-15.0);
      glRotatef(-40,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,1.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-5.9,0.0,15.8);
      glRotatef(140,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,1.0,0.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(9.5,0.0,7.5);
      glRotatef(220,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,0.0,1.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-9.5,0.0,-7.5);
      glRotatef(50,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(1.0,0.0,1.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(8.8,0.0,-7.5);
      glRotatef(-50,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,1.0,1.0);
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-8.5,0.0,8.4);
      glRotatef(-220,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower(0.0,1.0,1.0);
      glPopMatrix();

      //draw court gate

      glPushMatrix();
      glTranslatef(0,0.0,13.0);
      glRotatef(0,0,1,0);
      glScalef(1.2,1.2,1.2);
      drawCourtGate();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(0.0,0.0,-13.0);
      glRotatef(0,0,1,0);
      glScalef(1.2,1.2,1.2);
      drawCourtGate();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(2,0.0,13.0);
      glRotatef(0,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtGate();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-2.0,0.0,13.0);
      glRotatef(0,0,1,0);
      glScalef(1.0,1.0,1.0);
      drawCourtGate();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(2.0,0.0,-13.0);
      glRotatef(0,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtGate();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-2.0,0.0,-13.0);
      glRotatef(0,0,1,0);
      glScalef(1.0,1.0,1.0);
      drawCourtGate();
      //glColor3f(1.0,1.0,1.0);
      glPopMatrix();

      //draw marks
      glColor3f(1.0,1.0,1.0);
      glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1.0,1.0,1.0);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glScalef(0.5,1.0,1.0);
   //draw the ground
   glLineWidth(5);
      glBegin(GL_LINES);
      
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glVertex3f(x*16,0.01,z*16);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glVertex3f(x*16,0.01,z*16);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glVertex3f(x*16,0.01,z*16);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glVertex3f(x*16,0.01,z*16);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1.0,1.0,1.0);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   //draw the ground
      glBegin(GL_LINES);
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glVertex3f(x*3,0.01,z*3);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glVertex3f(x*3,0.01,z*3);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glVertex3f(x*3,0.01,z*3);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glVertex3f(x*3,0.01,z*3);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1.0,1.0,1.0);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   //draw the ground
      glBegin(GL_LINES);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(8,0.01,0);
      glVertex3f(-8.0,0.01,0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

}
/* draw branches
 *   for forbidden forest
 */

void drawBranches(float branchLength, float pointx, float pointy, float pointz){
    float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   float x,z;
   float YtoLowerHeight;  //y component for the NVects of the higher part
   
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[10]);
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:
      int i=0;
      
      //y is constant when the height is same 
      YtoLowerHeight = (0.5-0.0) / (0.0-branchLength) * (0.5-0.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f((float)i/(float)NumOfEdges,0); glVertex3f(x*0.5,0.0,z*0.5);
         //same x,z and NVect:
         glTexCoord2f((float)i/(float)NumOfEdges,pointy); glVertex3f(pointx,pointy,pointz);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f((float)i/(float)NumOfEdges,0); glVertex3f(x*0.5,0.0,z*0.5);
      //same x,z and NVect:
      glTexCoord2f((float)i/(float)NumOfEdges,pointy); glVertex3f(pointx,pointy,pointz);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
}

/* draw a tree
 *   for forbidden forest
 */
void drawTree(int treeHeight, float r1, float r2){
    float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   float x,z;
   float YtoLowerHeight;  //y component for the NVects of the higher part
   
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[10]);
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:
      int i=0;
      
      //y is constant when the height is same 
      YtoLowerHeight = (r1-r2) / (0.0-treeHeight) * (r1-r2);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f((float)i/(float)NumOfEdges*4,0); glVertex3f(x*r1,0.0,z*r1);
         //same x,z and NVect:
         glTexCoord2f((float)i/(float)NumOfEdges*4,treeHeight); glVertex3f(x*r2,treeHeight,z*r2);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f((float)i/(float)NumOfEdges*4,0); glVertex3f(x*r1,0.0,z*r1);
      //same x,z and NVect:
      glTexCoord2f((float)i/(float)NumOfEdges*4,treeHeight); glVertex3f(x*r2,treeHeight,z*r2);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      // draw braches
      //srand((unsigned)time(NULL));
      int n = 6; // number of branches
      int rand_height[6];
      for (int i=0; i<n; i++)         // generate random numbers
        rand_height[i]=(float)treeHeight/(float)(n+1)*(float)i;    //make the range of random number between 0 to 32767
      

      for (int j=0; j<n; j++){
      glPushMatrix();
       glTranslatef(0,rand_height[j],0);
       glScalef(1-0.02*n,1-0.02*n,1-0.02*n);
       drawBranches(3, Cos((float)j/(float)n * 360.0)*2, sqrt(9-4) ,Sin((float)j/(float)n * 360.0)*2);
       glPopMatrix();

       //draw the higher part of the tree
       glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[10]);
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:
      int i=0;
      
      //y is constant when the height is same 
      YtoLowerHeight = (r2-0.0) / (0.0-treeHeight*0.1) * (r2-0.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f((float)i/(float)NumOfEdges*2,0); glVertex3f(x*r2,treeHeight,z*r2);
         //same x,z and NVect:
         glTexCoord2f((float)i/(float)NumOfEdges*2,treeHeight*0.1); glVertex3f(0,treeHeight+treeHeight*0.1,0);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f((float)i/(float)NumOfEdges*2,0); glVertex3f(x*r2,treeHeight,z*r2);
      //same x,z and NVect:
      glTexCoord2f((float)i/(float)NumOfEdges*2,treeHeight*0.1); glVertex3f(0,treeHeight+treeHeight*0.1,0);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);


     }
}
void drawForest(int numberOfTrees){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
    // draw trees
      //srand((unsigned)time(NULL));
   int rand_x[numberOfTrees], rand_y[numberOfTrees];
   rand_x[0] = 7; rand_y[0] = -9;
   rand_x[1] = -3; rand_y[1] = 8;
   rand_x[2] = -0; rand_y[2] = -2;
   rand_x[3] = 4; rand_y[3] = 8;
   rand_x[4] = 3; rand_y[4] = 9;
   rand_x[5] = 0; rand_y[5] = 5;
   rand_x[6] = 2; rand_y[6] = 2;
   rand_x[7] = 7; rand_y[7] = 3;
   rand_x[8] = -7;rand_y[8] = -9;
   rand_x[9] = 0; rand_y[9] = 7;
   rand_x[10] = -3; rand_y[10] = 3;
   rand_x[11] = -9; rand_y[11] = 8;
   rand_x[12] = 0; rand_y[12] = 5;
   rand_x[13] = 9; rand_y[13] = 6;
   rand_x[14] = -6; rand_y[14] = 7;
   rand_x[15] = 1; rand_y[15] = 3;
   rand_x[16] = 0; rand_y[16] = 9;
   rand_x[17] = -5; rand_y[17] = -1;
   rand_x[18] = 2; rand_y[18] = 2;
   rand_x[19] = 6; rand_y[19] = 6;
   rand_x[20] = 9; rand_y[20] = 9;
   rand_x[21] = -4; rand_y[21] = -9;
   rand_x[22] = 5; rand_y[22] = 7;
   rand_x[23] = 7; rand_y[23] = -5;
   rand_x[24] = 1; rand_y[24] = 3;
   rand_x[25] = 9; rand_y[25] = 0;
   rand_x[26] = 5; rand_y[26] = -7;
   rand_x[27] = -7; rand_y[27] = 8;
   rand_x[28] = 6; rand_y[28] = 1;
   rand_x[29] = 3; rand_y[29] = -4;



      for (int j=0; j<numberOfTrees; j++){
      glPushMatrix();
       glTranslatef((float)(rand_x[j]),0.0,(float)(rand_y[j]));
       drawTree(rand_x[j]%4+10,1,0.5);
       glPopMatrix();
     }

   
}

/* draw the teeth of the wall
 *         the amount of teeth is TeethNumber
 */
void drawWallTeeth(int WallTeethNumber)
{
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
    glColor3f(1,1,1);
   glBindTexture(GL_TEXTURE_2D,texture[0]);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);

    glBegin(GL_QUADS);
    for(int i=0; i<WallTeethNumber; i++){
      glNormal3f(0.0,0.0,-1.0);
      glTexCoord2f(i*2.0*WallTeethSize/4,0.0); glVertex3f(i*2.0*WallTeethSize, 0.0, 0.0);
      glTexCoord2f(i*2.0*WallTeethSize/4,WallTeethSize/4);  glVertex3f(i*2.0*WallTeethSize, WallTeethSize, 0.0);
      glTexCoord2f((i*2.0+1.0)*WallTeethSize/4,WallTeethSize/4); glVertex3f((i*2.0+1.0)*WallTeethSize, WallTeethSize, 0.0);
      glTexCoord2f((i*2.0+1.0)*WallTeethSize/4,0.0); glVertex3f((i*2.0+1.0)*WallTeethSize, 0.0, 0.0);

    }
    glEnd();
     glPopMatrix();
  glDisable(GL_TEXTURE_2D);

    glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);
     glBegin(GL_QUADS);
    for(int i=0; i<WallTeethNumber; i++){
      glNormal3f(0.0,0.0,1.0);
      glTexCoord2f(i*2.0*WallTeethSize/4,0.0); glVertex3f(i*2.0*WallTeethSize, 0.0, 0.01);
      glTexCoord2f(i*2.0*WallTeethSize/4,WallTeethSize/4); glVertex3f(i*2.0*WallTeethSize, WallTeethSize, 0.01);
      glTexCoord2f((i*2.0+1.0)*WallTeethSize/4,WallTeethSize/4); glVertex3f((i*2.0+1.0)*WallTeethSize, WallTeethSize, 0.01);
      glTexCoord2f((i*2.0+1.0)*WallTeethSize/4,0.0); glVertex3f((i*2.0+1.0)*WallTeethSize, 0.0, 0.01);

    }
    glEnd();
     glPopMatrix();
  glDisable(GL_TEXTURE_2D);
}

/* draw the teeth of the tower
*/
void drawTowerTeeth()
{
  
     float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   int i=0;
   float x,z;
   while (i < NumOfEdges)
      {  
         glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


         glBegin(GL_QUADS);   
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x/4,0.0,z/4); glVertex3f(x,0.0,z);
         //same x,z and NVect:
         glTexCoord3f(x/4,TowerTeethSize/4,z/4); glVertex3f(x,TowerTeethSize,z);



         x = Cos((float)(i+3)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+3)/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x/4,TowerTeethSize/4,z/4); glVertex3f(x,TowerTeethSize,z);
         glTexCoord3f(x/4,0.0,z/4); glVertex3f(x,0.0,z);
         glEnd();
         glPopMatrix();
         glDisable(GL_TEXTURE_2D);

         glPushMatrix();
         glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
         glBegin(GL_QUADS);   
         x = 0.99*Cos((float)(i)/(float)NumOfEdges * 360.0);
         z = 0.99*Sin((float)(i)/(float)NumOfEdges * 360.0);
         glNormal3f(-x,0.0,-z);
         glTexCoord3f(x/4,0.0,z/4); glVertex3f(x,0.0,z);
         //same x,z and NVect:
         glTexCoord3f(x/4,TowerTeethSize/4,z/4); glVertex3f(x,TowerTeethSize,z);

         x = 0.99*Cos((float)(i+3)/(float)NumOfEdges * 360.0);
         z = 0.99*Sin((float)(i+3)/(float)NumOfEdges * 360.0);
         glNormal3f(-x,0.0,-z);
         glTexCoord3f(x/4,TowerTeethSize/4,z/4); glVertex3f(x,TowerTeethSize,z);
         glTexCoord3f(x/4,0.0,z/4); glVertex3f(x,0.0,z);
         glEnd();
         glPopMatrix();
         glDisable(GL_TEXTURE_2D);

         i+=14;
         
      }



}

/* draw the wall, together with the teeth
 *      with the length of the wall
 */
void drawWall(GLfloat length){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

 glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);
   

   glBegin(GL_QUADS);
   
   //glBindTexture(GL_TEXTURE_2D,texture[1]);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(0.0,0.0); glVertex3f(0.0, 0.0, 0.0);
   glTexCoord2f(0.0,WallHeight/4); glVertex3f(0.0, WallHeight, 0.0);
   glTexCoord2f(length/4,WallHeight/4); glVertex3f(length, WallHeight, 0.0);
   glTexCoord2f(length/4,0.0); glVertex3f(length, 0.0, 0.0);

   glEnd();
 glPopMatrix();
  glDisable(GL_TEXTURE_2D);
//  to see whether the length of teeth exceeds the wall length
   int i=(int)(length/WallTeethSize/2);
   if (i*WallTeethSize > length) 
      i=i-1;

   glPushMatrix();
      glTranslatef(0.0, WallHeight, 0.0);
      drawWallTeeth(i);
   glPopMatrix();

// draw the double wall
   glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(0.0,0.0); glVertex3f(0.0, 0.0, 1.0);
   glTexCoord2f(0.0,WallHeight/4); glVertex3f(0.0, WallHeight, 1.0);
   glTexCoord2f(length/4,WallHeight/4); glVertex3f(length, WallHeight, 1.0);
   glTexCoord2f(length/4,0.0); glVertex3f(length, 0.0, 1.0);

   glEnd();

glPopMatrix();
  glDisable(GL_TEXTURE_2D);
//  to see whether the length of teeth exceeds the wall length
    i=(int)(length/WallTeethSize/2);
   if (i*WallTeethSize > length) 
      i=i-1;

   glPushMatrix();

      
      glTranslatef(0.0, WallHeight, 1.0);
      glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);
      glBegin(GL_QUADS);
    for(int j=0; j<i; j++){
      glNormal3f(0.0,0.0,1.0);
      glTexCoord2f(j*2.0*WallTeethSize/4,0.0); glVertex3f(j*2.0*WallTeethSize, 0.0, 0.0);
      glTexCoord2f(j*2.0*WallTeethSize/4,WallTeethSize/4); glVertex3f(j*2.0*WallTeethSize, WallTeethSize, 0.0);
      glTexCoord2f((j*2.0+1.0)*WallTeethSize/4,WallTeethSize/4); glVertex3f((j*2.0+1.0)*WallTeethSize, WallTeethSize, 0.0);
      glTexCoord2f((j*2.0+1.0)*WallTeethSize/4,0.0); glVertex3f((j*2.0+1.0)*WallTeethSize, 0.0, 0.0);

    }
    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

     i=(int)(length/WallTeethSize/2);
   if (i*WallTeethSize > length) 
      i=i-1;

   glPushMatrix();
      
      glTranslatef(0.0, WallHeight, 1.0);
      glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
  glColor3f(0.465,0.531,0.598);
   glBindTexture(GL_TEXTURE_2D,texture[0]);
      glBegin(GL_QUADS);
    for(int j=0; j<i; j++){
      glNormal3f(0.0,0.0,-1.0);
      glTexCoord2f(j*2.0*WallTeethSize/4,0.0); glVertex3f(j*2.0*WallTeethSize, 0.0, -0.01);
      glTexCoord2f(j*2.0*WallTeethSize/4,WallTeethSize/4); glVertex3f(j*2.0*WallTeethSize, WallTeethSize, -0.01);
      glTexCoord2f((j*2.0+1.0)*WallTeethSize/4,WallTeethSize/4); glVertex3f((j*2.0+1.0)*WallTeethSize, WallTeethSize, -0.01);
      glTexCoord2f((j*2.0+1.0)*WallTeethSize/4,0.0); glVertex3f((j*2.0+1.0)*WallTeethSize, 0.0, -0.01);

    }
    glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
   
// draw the higher ground
   glPushMatrix();
    glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[4]);

   glBegin(GL_QUADS);
   
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(0.0,0.0);glVertex3f(0.0, WallHeight, 0.0);
   glTexCoord2f(0.0,1.0/4); glVertex3f(0.0, WallHeight, 1.0);
   glTexCoord2f(length/4,1.0/4); glVertex3f(length, WallHeight, 1.0);
   glTexCoord2f(length/4,0.0);glVertex3f(length, WallHeight, 0.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

// draw the lower ground
   glBegin(GL_QUADS);
   glNormal3f(0.0,-1.0,0.0);
   glVertex3f(0.0, 0.0, 0.0);
   glVertex3f(0.0, 0.0, 1.0);
   glVertex3f(length, 0.0, 1.0);
   glVertex3f(length, 0.0, 0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();

}

/* draw the tower, together with the teeh
 */

void drawTower(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   float x,z;
   float YtoLowerHeight;  //y component for the NVects of the higher part
   
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[2]);
   
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:
      int i=0;
      
      //y is constant when the height is same 
      YtoLowerHeight = (Radius-1.0) / (0.0-LowerHeight) * (Radius-1.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(x,LowerHeight,z);

        
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(x,LowerHeight,z);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      //Create the higher part of the tower
      glPushMatrix();
      glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.688,0.766,0.867);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   
    
   glBegin(GL_QUAD_STRIP);
      for ( i = 0; i < NumOfEdges; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges*360);
         z = Sin((float)i/(float)NumOfEdges*360);
         glNormal3f(x,0.0,z);
         glTexCoord2f((float)i/(float)NumOfEdges*2,LowerHeight/4); glVertex3f(x,LowerHeight,z);
         glTexCoord2f((float)i/(float)NumOfEdges*2,HigherHeight/4); glVertex3f(x,HigherHeight,z);
      }
      x = Cos((float)i/(float)NumOfEdges*360);
      z = Sin((float)i/(float)NumOfEdges*360);
      glNormal3f(x,0.0,z);
      glTexCoord2f((float)i/(float)NumOfEdges*2,LowerHeight/4); glVertex3f(x,LowerHeight,z);
      glTexCoord2f((float)i/(float)NumOfEdges*2,HigherHeight/4); glVertex3f(x,HigherHeight,z);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
   
      //draw the higherground
       glPushMatrix();
    glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[4]);
      glBegin(GL_TRIANGLES);
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = cos((float)i/(float)NumOfEdges * 360.0);
         z = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glTexCoord2f(0,0); glVertex3f(0.0,HigherHeight,0.0);
         glTexCoord2f(x/4,z/4); glVertex3f(x,HigherHeight,z);

         x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,1.0,0.0);
         glTexCoord2f(x/4,z/4); glVertex3f(x,HigherHeight,z);
      }
      x = cos((float)i/(float)NumOfEdges * 360.0);
      z = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glTexCoord2f(0,0); glVertex3f(0.0,HigherHeight,0.0);
      glTexCoord2f(x/4,z/4); glVertex3f(x,HigherHeight,z);
      x = cos(1.0/(float)NumOfEdges * 360.0);
      z = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,1.0,0.0);
      glTexCoord2f(x/4,z/4); glVertex3f(x,HigherHeight,z);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);



      //draw the ground
      glBegin(GL_TRIANGLES);
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = cos((float)i/(float)NumOfEdges * 360.0);
         z = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(0.0,0.0,0.0);
         glVertex3f(x*Radius,0.0,z*Radius);

         x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*Radius,0.0,z*Radius);
      }
      x = cos((float)i/(float)NumOfEdges * 360.0);
      z = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(x*Radius,0.0,z*Radius);
      x = cos(1.0/(float)NumOfEdges * 360.0);
      z = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*Radius,0.0,z*Radius);
      glEnd();



      // draw the teeth
      glPushMatrix();
      glTranslatef(0.0, HigherHeight, 0.0);
      drawTowerTeeth();
      glPopMatrix();
glColor3f(1.0,1.0,1.0);
}


/*
 * draw the board
 */

void drawBoard(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[7]);

   glBegin(GL_QUADS);
   
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,0.1,0.0);
   glTexCoord2f(2.0/4,0); glVertex3f(2.0,0.1,0.0);
   glTexCoord2f(2.0/4,-HigherHeight/4); glVertex3f(2.0,0.1,-HigherHeight);
   glTexCoord2f(0,-HigherHeight/4); glVertex3f(0.0,0.1,-HigherHeight);
   glEnd();

    glBegin(GL_QUADS);
   //glColor3f(0.543,0.270,0.074);
   glNormal3f(0.0,-1.0,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,0.09,0.0);
   glTexCoord2f(2.0/4,0); glVertex3f(2.0,0.09,0.0);
   glTexCoord2f(2.0/4,-HigherHeight/4); glVertex3f(2.0,0.09,-HigherHeight);
   glTexCoord2f(0,-HigherHeight/4); glVertex3f(0.0,0.09,-HigherHeight);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glBegin(GL_LINES);
   glVertex3f(0.0,0.0,-(HigherHeight-0.5));
   glVertex3f(0.0,HigherHeight/2+0.5,0.0);
   glEnd();

   glBegin(GL_LINES);
   glVertex3f(2.0,0.0,-(HigherHeight-0.5));
   glVertex3f(2.0,HigherHeight/2+0.5,0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();


}


/*
 *  Draw the chair
 */
void drawchair(){
      glPushMatrix();
      glScalef(1,1,0.5);
      cube1(0.6,0.15,0.55,0.1,0.08,0.25,0);
      cube1(0.6,0.2,0.55,0.15,0.03,0.35,0);
      
      cube1(0.6-0.15,0.4,0.55,0.03,0.21,0.35,0);
      glPopMatrix();
}

/*
 *  Draw the chair
 */
void drawcdesk(){
      glPushMatrix();
      glScalef(1,1,0.5);
      cube1(0.6,0.15,0.55,0.1,0.08,0.25,0);
      cube1(0.6,0.2,0.55,0.15,0.03,0.35,0);
      
      //cube1(0.6-0.15,0.4,0.55,0.03,0.21,0.35,0);
      glPopMatrix();
}
/*
 *  Draw the gate
 */
void drawGate(){
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.0,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);


 glPushMatrix();
   glTranslatef(8.0,0.55,0.1);
   //glRotatef(-45,0,1,0);
   glScalef(2.0,2.0,2.0);
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texture[12]);
   glutSolidTeapot(0.05);
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();



   glPushMatrix();
   glTranslatef(6.0,0.0,-1.0);
   glRotatef(-45,0,1,0);
   glScalef(2.0,2.0,2.0);

   drawchair();
   glPopMatrix();

  glPushMatrix();
   glTranslatef(7.5,0.0,-1.0);
   glRotatef(-45,0,1,0);
   glScalef(2.0,2.0,2.0);

   drawcdesk();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(11.0,0.0,0.0);
   glRotatef(225,0,1,0);
   glScalef(2.0,2.0,2.0);



   drawchair();
   glPopMatrix();

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(3.5,0.0,-2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(3.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(3.5,HigherHeight*1.2,2.0);
   glTexCoord2f(1*1.5,0); glVertex3f(3.5,0.0,2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(3.501,0.0,-2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(3.501,HigherHeight*1.2,-2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(3.501,HigherHeight*1.2,2.0);
   glTexCoord2f(1*1.5,0); glVertex3f(3.501,0.0,2.0);

   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(3.5,0.0,-2.0);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(3.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(3.5,0.0,-1.999);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-1.999);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-1.999);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(3.5,HigherHeight*1.2,-1.999);
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0);glVertex3f(4.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(4.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(4.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(4.5,0.0,-2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   
   glTexCoord2f(0,0);glVertex3f(4.501,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(4.501,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(4.501,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(4.501,0.0,-2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-4.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,-4.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-3.999);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,-3.999);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,-3.999);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-3.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(6.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(6.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(6.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.5,0.0,-2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(6.499,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(6.499,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(6.499,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.499,0.0,-2.0);
   
   
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);


   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,-2.0);
   glTexCoord2f(7.5/4,0);glVertex3f(7.5,0.0,-2.0);
   glTexCoord2f(7.5/4,HigherHeight*1.2/4);glVertex3f(7.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,-1.999);
   glTexCoord2f(7.5/4,0);glVertex3f(7.5,0.0,-1.999);
   glTexCoord2f(7.5/4,HigherHeight*1.2/4);glVertex3f(7.5,HigherHeight*1.2,-1.999);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,-1.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(7.5/4,HigherHeight/8);glVertex3f(7.5,HigherHeight/2,-2.0);
   glTexCoord2f(9.5/4,HigherHeight/8);glVertex3f(9.5,HigherHeight/2,-2.0);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(7.5/4,HigherHeight*1.2/4);glVertex3f(7.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(7.5/4,HigherHeight/8);glVertex3f(7.5,HigherHeight/2,-1.999);
   glTexCoord2f(9.5/4,HigherHeight/8);glVertex3f(9.5,HigherHeight/2,-1.999);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-1.999);
   glTexCoord2f(7.5/4,HigherHeight*1.2/4);glVertex3f(7.5,HigherHeight*1.2,-1.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-2.0);
   glTexCoord2f(9.5/4,0);glVertex3f(9.5,0.0,-2.0);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-2.0);
   glEnd();
  

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-1.999);
   glTexCoord2f(9.5/4,0);glVertex3f(9.5,0.0,-1.999);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-1.999);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-1.999);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0,0);glVertex3f(10.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(10.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(10.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(10.5,0.0,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0);glVertex3f(10.501,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(10.501,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(10.501,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(10.501,0.0,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-4.0);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-4.0);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-3.999);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-3.999);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-3.999);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-3.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0);glVertex3f(12.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(12.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(12.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(12.5,0.0,-2.0);
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0,0);glVertex3f(12.499,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(12.499,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(12.499,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(12.499,0.0,-2.0);
   
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-2.0);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,-2.0);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-1.999);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,-1.999);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,-1.999);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-1.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(1*1.5,0);glVertex3f(13.5,0.0,-2.0);
   glTexCoord2f(1*1.5,1*1.5);glVertex3f(13.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(0,1*1.5);glVertex3f(13.5,HigherHeight*1.2,2.0);
   glTexCoord2f(0,0); glVertex3f(13.5,0.0,2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(1*1.5,0);glVertex3f(13.499,0.0,-2.0);
   glTexCoord2f(1*1.5,1*1.5);glVertex3f(13.499,HigherHeight*1.2,-2.0);
   glTexCoord2f(0,1*1.5);glVertex3f(13.499,HigherHeight*1.2,2.0);
   glTexCoord2f(0,0); glVertex3f(13.499,0.0,2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,2.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,2.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,2.0);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,1.999);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,1.999);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,1.999);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,1.999);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
    glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0.0,0);glVertex3f(10.5,0.0,2.0);
   glTexCoord2f(0.0,1.5);glVertex3f(10.5,HigherHeight*1.2,2.0);
   glTexCoord2f(1/1.7,1.5);glVertex3f(10.5,HigherHeight*1.2,4.0);
   glTexCoord2f(1/1.7,0);glVertex3f(10.5,0.0,4.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0.0,0);glVertex3f(10.499,0.0,2.0);
   glTexCoord2f(0.0,1.5);glVertex3f(10.499,HigherHeight*1.2,2.0);
   glTexCoord2f(1/1.7,1.5);glVertex3f(10.499,HigherHeight*1.2,4.0);
   glTexCoord2f(1/1.7,0);glVertex3f(10.499,0.0,4.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,4.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,4.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,4.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,3.999);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,3.999);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,3.999);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,3.999);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0.0,0.0);glVertex3f(6.5,0.0,4.0);
   glTexCoord2f(0,1.5);glVertex3f(6.5,HigherHeight*1.2,4.0);
   glTexCoord2f(1/1.7,1.5);glVertex3f(6.5,HigherHeight*1.2,2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.5,0.0,2.0);
   
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0.0,0.0);glVertex3f(6.501,0.0,4.0);
   glTexCoord2f(0,1.5);glVertex3f(6.501,HigherHeight*1.2,4.0);
   glTexCoord2f(1/1.7,1.5);glVertex3f(6.501,HigherHeight*1.2,2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.501,0.0,2.0);
   
   
   
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,2.0);
   glTexCoord2f(3.5/4,0);glVertex3f(3.5,0.0,2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(3.5,HigherHeight*1.2,2.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,1.999);
   glTexCoord2f(3.5/4,0);glVertex3f(3.5,0.0,1.999);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(3.5,HigherHeight*1.2,1.999);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,1.999);
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   // draw the ground
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(3.5/4,-2.0/4); glVertex3f(3.5,0.0,-2.0);
   glTexCoord2f(13.5/4,-2.0/4); glVertex3f(13.5,0.0,-2.0);
   glTexCoord2f(13.5/4,2.0/4); glVertex3f(13.5,0.0,2.0);
   glTexCoord2f(3.5/4,2.0/4); glVertex3f(3.5,0.0,2.0);
   glEnd();


   glBegin(GL_POLYGON);

   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(10.5/4,4.0/4); glVertex3f(10.5,0.0,4.0);
   glTexCoord2f(10.5/4,2.0/4); glVertex3f(10.5,0.0,2.0);
   glTexCoord2f(6.5/4,2.0/4); glVertex3f(6.5,0.0,2.0);
   glTexCoord2f(6.5/4,4.0/4); glVertex3f(6.5,0.0,4.0);
 
 
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(4.5/4,-2.0/4); glVertex3f(4.5,0.0,-2.0);
   glTexCoord2f(4.5/4,-4.0/4); glVertex3f(4.5,0.0,-4.0);
   glTexCoord2f(6.5/4,-4.0/4); glVertex3f(6.5,0.0,-4.0);
   glTexCoord2f(6.5/4,-2.0/4); glVertex3f(6.5,0.0,-2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(10.5/4,-2.0/4); glVertex3f(10.5,0.0,-2.0);
   glTexCoord2f(10.5/4,-4.0/4); glVertex3f(10.5,0.0,-4.0);
   glTexCoord2f(12.5/4,-4.0/4); glVertex3f(12.5,0.0,-4.0);
   glTexCoord2f(12.5/4,-2.0/4); glVertex3f(12.5,0.0,-2.0);
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   // draw the higher ground
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(3.5/4,-2.0/4); glVertex3f(3.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(13.5/4,-2.0/4); glVertex3f(13.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(13.5/4,2.0/4); glVertex3f(13.5,HigherHeight*1.2,2.0);
   glTexCoord2f(3.5/4,2.0/4); glVertex3f(3.5,HigherHeight*1.2,2.0);
   glEnd();

   glBegin(GL_POLYGON);

   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(10.5/4,4.0/4); glVertex3f(10.5,HigherHeight*1.2,4.0);
   glTexCoord2f(10.5/4,2.0/4); glVertex3f(10.5,HigherHeight*1.2,2.0);
   glTexCoord2f(6.5/4,2.0/4); glVertex3f(6.5,HigherHeight*1.2,2.0);
   glTexCoord2f(6.5/4,4.0/4); glVertex3f(6.5,HigherHeight*1.2,4.0);
 
 
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(4.5/4,-2.0/4); glVertex3f(4.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(4.5/4,-4.0/4); glVertex3f(4.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(6.5/4,-4.0/4); glVertex3f(6.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(6.5/4,-2.0/4); glVertex3f(6.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(10.5/4,-2.0/4); glVertex3f(10.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(10.5/4,-4.0/4); glVertex3f(10.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(12.5/4,-4.0/4); glVertex3f(12.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(12.5/4,-2.0/4); glVertex3f(12.5,HigherHeight*1.2,-2.0);
   glEnd();
   glColor3f(1.0,1.0,1.0);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

void drawhousetower(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

 glPushMatrix();
      glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1.0,0.867,0.676);
   glBindTexture(GL_TEXTURE_2D,texture[6]);
   
      glBegin(GL_QUAD_STRIP);
      int i;
      float x,z;
      for ( i = 0; i < NumOfEdges; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges*360);
         z = Sin((float)i/(float)NumOfEdges*360);
         glNormal3f(x,0.0,z);
         glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight); glVertex3f(x*0.8,HouseHeight,z*0.8);
         glTexCoord2f((float)i/(float)NumOfEdges*4,0.0); glVertex3f(x*0.8,0.0,z*0.8);
      }
      x = Cos((float)i/(float)NumOfEdges*360);
      z = Sin((float)i/(float)NumOfEdges*360);
      glNormal3f(x,0.0,z);
      glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight); glVertex3f(x*0.8,HouseHeight,z*0.8);
      glTexCoord2f((float)i/(float)NumOfEdges*4,0.0); glVertex3f(x*0.8,0.0,z*0.8);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
   
   //draw higherpart
      glPushMatrix();
      glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[6]);
     
   glBegin(GL_QUAD_STRIP);
      for ( i = 0; i < NumOfEdges; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges*360);
         z = Sin((float)i/(float)NumOfEdges*360);
         glNormal3f(x,0.0,z);
         glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+2.0); glVertex3f(x,HouseHeight+2.0,z);
         glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight); glVertex3f(x,HouseHeight,z);
      }
      x = Cos((float)i/(float)NumOfEdges*360);
      z = Sin((float)i/(float)NumOfEdges*360);
      glNormal3f(x,0.0,z);
      glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+2.0); glVertex3f(x,HouseHeight+2.0,z);
      glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight); glVertex3f(x,HouseHeight,z);

      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

   //draw a point
       float YtoPoint = (1.0-0.0) / 2.0 * (1.0-0.0);
       glPushMatrix();
       glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.543,0.270,0.074);
   glBindTexture(GL_TEXTURE_2D,texture[1]);
    
   glBegin(GL_QUAD_STRIP);
      //Create the lower part of the tower:

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoPoint,z);
         glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+4.0); glVertex3f(0.0,HouseHeight+4.0,0.0);
         //same x,z and NVect:
         glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+2.0); glVertex3f(x,HouseHeight+2.0,z);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoPoint,z);
      glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+4.0); glVertex3f(0.0,HouseHeight+4.0,0.0);
      //same x,z and NVect:
      glTexCoord2f((float)i/(float)NumOfEdges*4,HouseHeight+2.0); glVertex3f(x,HouseHeight+2.0,z);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

}


void drawhouse(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.820,0.410,0.117);
   glBindTexture(GL_TEXTURE_2D,texture[5]);

   //draw the central part
   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,0.0,0.0);
   glTexCoord2f(0,HouseHeight/4); glVertex3f(0.0,HouseHeight,0.0);
   glTexCoord2f(3.5/4,HouseHeight/4); glVertex3f(3.5,HouseHeight,0.0);
   glTexCoord2f(3.5/4,0); glVertex3f(3.5,0.0,0.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(3.5/4,3.0/4); glVertex3f(3.5,3.0,0.0);
   glTexCoord2f(3.5/4,HouseHeight/4); glVertex3f(3.5,HouseHeight,0.0);
   glTexCoord2f(5.5/4,HouseHeight/4); glVertex3f(5.5,HouseHeight,0.0);
   glTexCoord2f(5.5/4,3.0/4); glVertex3f(5.5,3.0,0.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(5.5/4,0); glVertex3f(5.5,0.0,0.0);
   glTexCoord2f(5.5/4,HouseHeight/4); glVertex3f(5.5,HouseHeight,0.0);
   glTexCoord2f(9.0/4,HouseHeight/4); glVertex3f(9.0,HouseHeight,0.0);
   glTexCoord2f(9.0/4,0); glVertex3f(9.0,0.0,0.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,1.0);
   glTexCoord2f(0.0,0.0);glVertex3f(9.0,0.0,0.0);
   glTexCoord2f(0.0,1*1.5); glVertex3f(9.0,HouseHeight,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(9.0,HouseHeight,5.0);
   glTexCoord2f(1*1.5,0); glVertex3f(9.0,0.0,5.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(9.0/4,0); glVertex3f(9.0,0.0,5.0);
   glTexCoord2f(9.0/4,HouseHeight/4); glVertex3f(9.0,HouseHeight,5.0);
   glTexCoord2f(0,HouseHeight/4); glVertex3f(0.0,HouseHeight,5.0);
   glTexCoord2f(0,0);glVertex3f(0.0,0.0,5.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,0.0,5.0);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,HouseHeight,5.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(0.0,HouseHeight,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,0.0,0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);


   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[5]);
   glBegin(GL_QUADS);
   glNormal3f(0.0,-1.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,0.0,0.0);
   glTexCoord2f(9.0/4,0); glVertex3f(9.0,0.0,0.0);
   glTexCoord2f(9.0/4,5.0/4); glVertex3f(9.0,0.0,5.0);
   glTexCoord2f(0,5.0/4); glVertex3f(0.0,0.0,5.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,1.0,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,HouseHeight,0.0);
   glTexCoord2f(9.0/4,0); glVertex3f(9.0,HouseHeight,0.0);
   glTexCoord2f(9.0/4,5.0/4); glVertex3f(9.0,HouseHeight,5.0);
   glTexCoord2f(0,5.0/4); glVertex3f(0.0,HouseHeight,5.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   
   
}

//draw a flag
void drawflag(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   glLineWidth(1);
   glBegin(GL_LINES);

   
   glColor3f(0.543,0.270,0.074);
   glVertex3f(0.0,15.0,0.0);
   glVertex3f(0.0,16.5,0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();

   glBegin(GL_LINES);
   glVertex3f(0.0,16.5,0.0);
   glVertex3f(0.0,20.0,0.0);
   glEnd();

   
   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glColor3f(1.0,1.0,1.0);
   glColor3f(1.0,0.0,0.0);
   glVertex3f(0.0,19,0.0);
   glVertex3f(0.0,20.0,0.0);
   glVertex3f(2.0,20.0,0.0);
   glVertex3f(2.0,19.0,0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glColor3f(1.0,1.0,1.0);
   glColor3f(1.0,0.0,0.0);
   glVertex3f(0.0,19,-0.01);
   glVertex3f(0.0,20.0,-0.01);
   glVertex3f(2.0,20.0,-0.01);
   glVertex3f(2.0,19.0,-0.01);
   glColor3f(1.0,1.0,1.0);
   glEnd();


}

void drawWand(float pointx, float pointy, float pointz){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

   int i;
   float YtoPoint = (0.2-0.0) / sqrt(pointx*pointx+pointy*pointy+pointz*pointz) * (0.2-0.0);
   glBegin(GL_TRIANGLES);
      glColor3f(0.543,0.270,0.074);
      double wandx,wandz;
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         wandx = cos((float)i/(float)NumOfEdges * 360.0);
         wandz = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(wandx,YtoPoint,wandz);
         glVertex3f(pointx,pointy,pointz);
         glVertex3f(wandx*0.2,0.0,wandz*0.2);

         wandx = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         wandz = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(wandx,YtoPoint,wandz);
         glVertex3f(wandx*0.2,0.0,wandz*0.2);
      }
      wandx = cos((float)i/(float)NumOfEdges * 360.0);
      wandz = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(wandx,YtoPoint,wandz);
      glVertex3f(pointx,pointy,pointz);
      glVertex3f(wandx*0.2,0.0,wandz*0.2);
      wandx = cos(1.0/(float)NumOfEdges * 360.0);
      wandz = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(wandx,YtoPoint,wandz);
      glVertex3f(wandx*0.2,0.0,wandz*0.2);
      
      glColor3f(1.0,1.0,1.0);
      glEnd();

}

void drawGround(){
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.0,1.0,0.0);
   glBindTexture(GL_TEXTURE_2D,texture[9]);
   glBegin(GL_QUADS);
   glNormal3f(0.0,1.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(-500.0,-0.05,-500.0);
   glTexCoord2f(25*8,0); glVertex3f(500,-0.05,-500.0);
   glTexCoord2f(25*8,20*8); glVertex3f(500,-0.05,500);
   glTexCoord2f(0,20*8); glVertex3f(-500,-0.05,500);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,-1.0,0.0);
   glTexCoord2f(0,0); glVertex3f(-500.0,-0.1,-500.0);
   glTexCoord2f(25*8,0); glVertex3f(500,-0.1,-500.0);
   glTexCoord2f(25*8,20*8); glVertex3f(500,-0.1,500);
   glTexCoord2f(0,20*8); glVertex3f(-500.0,-0.1,500);
   //glColor3f(1.0,1.0,1.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

//I load an obj of bridge here but when I run the program it gets a little bit slow
//so I comment them out, but these code will work!
/*void drawBridge(){
  float RGBA[] = {1,1,1,1};
  float Emission[]  = {0.0,0.0,0.01*emission,1.0};
  glMaterialf(GL_FRONT,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT,GL_SPECULAR,RGBA);
   glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
  obj = LoadOBJ("woodenbridge.obj");
  glCallList(obj);
}*/

/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
   const double len=2.0;  //  Length of axes
   //  Erase the window and the depth buffer
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   //  Enable Z-buffering in OpenGL
   glEnable(GL_DEPTH_TEST);

   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective - set eye position
   if (mode_project==1)
   {
      double Ex = -2*dim*Sin(th)*Cos(ph);
      double Ey = +2*dim        *Sin(ph);
      double Ez = +2*dim*Cos(th)*Cos(ph);
      gluLookAt(Ex,Ey,Ez , 0,0,0 , 0,Cos(ph),0);
   }

   else if (mode_project==2){
      // set parameters for glulookat()
    fov = 40;
              
      FirstpersonNaviagtion();
   }


   //  Flat or smooth shading
   glShadeModel(smooth ? GL_SMOOTH : GL_FLAT);

   if(box) {
    Sky(3.5*dim);
  }
   else{
    Sky1(3.5*dim);
  }
   //  Light switch
   if (light)
   {
        //  Translate intensity to color vectors
        float Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0};
        float Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0};
        float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
        //  Light position
        float Position[]  = {distance*Cos(zh_l)-7.0,ylight,distance*Sin(zh_l),1.0};
        float Position_stick[]  = {Ex+0.5*Sin(zh),Ey-1.0+0.5*Sin(theta+10),Ez-0.5*Cos(zh),1.0};
       
        //float Position_fire[] = {2*Cos(zh_l)-7.0,ylight,2*Sin(zh_l),1.0};
        //  Draw light position as ball (still no lighting here)
        glColor3f(1,1,1);
        ball(Position[0],Position[1],Position[2] , 0.1);
        float Position_fire[] = {-6.9,-0.5,1.0,1.0};
   ball(Position_fire[0],Position_fire[1],Position_fire[2] , 0.01);   
   ball(Position_stick[0],Position_stick[1],Position_stick[2] , 0.01);

          //  OpenGL should normalize normal vectors
        glEnable(GL_NORMALIZE);
        //  Enable lighting
        glEnable(GL_LIGHTING);
        //  Location of viewer for specular calculations
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,local);
        //  glColor sets ambient and diffuse color materials
        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        //  Enable light 0
        glEnable(GL_LIGHT0);
        //  Set ambient, diffuse, specular components and position of light 0
        glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
        glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
        glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
        glLightfv(GL_LIGHT0,GL_POSITION,Position_fire);

       if(move_light){
       glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1,GL_AMBIENT ,Ambient);
        glLightfv(GL_LIGHT1,GL_DIFFUSE ,Diffuse);
        glLightfv(GL_LIGHT1,GL_SPECULAR,Specular);
        glLightfv(GL_LIGHT1,GL_POSITION,Position);
      }

      if(stick_light){
        glEnable(GL_LIGHT2);
        glLightfv(GL_LIGHT2,GL_AMBIENT ,Ambient);
        glLightfv(GL_LIGHT2,GL_DIFFUSE ,Diffuse);
        glLightfv(GL_LIGHT2,GL_SPECULAR,Specular);
        glLightfv(GL_LIGHT2,GL_POSITION,Position_stick);
      }
      else
        glDisable(GL_LIGHT2);
   }
   else
     glDisable(GL_LIGHTING);
   //  Set view angle
      
   glPushMatrix();
   glTranslatef(0.0,-0.5,0.0);


   //  Draw the stick
  if(mode_project== 2){
   // if(click_bat ==1 ){
  glPushMatrix();
   glTranslatef(Ex,Ey-1.0,Ez);
   glScalef(0.7*scale,0.7*scale,0.7*scale);
   drawWand(10*Sin(zh), 50*Sin(theta+10), -10*Cos(zh));
   glPopMatrix();
 //}//


}

//drawchair();

//partical engine
if(firework&&mode_project == 2){
  glPushMatrix();
 glTranslatef(Ex+10*Sin(zh)+13.8,Ey-1.0+50*Sin(theta+10)-7.0,Ez-10*Cos(zh)-0.0);
clock_t iNowTime = clock();
float timePassed = (float)(iNowTime- g_iLastRenderTime)/CLOCKS_PER_SEC;
g_ParticleSystem1.UpdateSystem(timePassed);
g_iLastRenderTime = iNowTime;
float zDist = Ez - g_ParticleSystem1.m_EmitterPosition.z;
  float xDist = Ex - g_ParticleSystem1.m_EmitterPosition.x;
  float CamDistToEmitter = sqrt(SQR(zDist)+SQR(xDist));
  if (CamDistToEmitter < 0.2f) //avoid too big particles
    CamDistToEmitter = 0.2f;
  glPointSize(1.0f/CamDistToEmitter);
  g_ParticleSystem1.Render();
  glPopMatrix();
}

glPushMatrix();
 glTranslatef(-7.0,0.0,0.5);
clock_t iNowTime1 = clock();
float timePassed1 = (float)(iNowTime1 - g_iLastRenderTime1)/CLOCKS_PER_SEC;
g_ParticleSystem2.UpdateSystem(timePassed1);
  g_ParticleSystem4.UpdateSystem(timePassed1);
  g_iLastRenderTime1 = iNowTime1;
  //glEnable(GL_TEXTURE_2D);
  float zDist = Ez - g_ParticleSystem1.m_EmitterPosition.z;
  float xDist = Ex - g_ParticleSystem1.m_EmitterPosition.x;
  float CamDistToEmitter = sqrt(SQR(zDist)+SQR(xDist));
  if (CamDistToEmitter < 0.2f) //avoid too big particles
    CamDistToEmitter = 0.2f;
  glPointSize(1.0f/CamDistToEmitter);

  glEnable(GL_TEXTURE_2D);
  g_ParticleSystem2.Render();
  g_ParticleSystem4.Render();
  glDisable(GL_TEXTURE_2D);
  glPopMatrix();


   //  Draw scene

 glScalef(2*scale,2*scale,2*scale);
 glTranslatef(-25,0,0);

 /*glPushMatrix();
 glScalef(0.45,0.45,0.45);
   glTranslatef(48.0,0.0,0.0);
   drawBridge();
   glPopMatrix();*/

 glPushMatrix();
 glScalef(0.8,0.8,0.8);
 glTranslatef(40.0,0.0,-5.0);
 drawCourt();
 glPopMatrix();

 //draw forest
 glPushMatrix();
 glTranslatef(25.0,0.0,25.0);
 drawForest(numberOfTrees);
 glPopMatrix();

 glPushMatrix();
 glTranslatef(35.0,0.0,25.0);
 drawForest(numberOfTrees);
 glPopMatrix();

 glPushMatrix();
 glTranslatef(35.0,0.0,35.0);
 drawForest(numberOfTrees);
 glPopMatrix();

 glPushMatrix();
 glTranslatef(25.0,0.0,35.0);
 drawForest(numberOfTrees);
 glPopMatrix();

 glPushMatrix();
 glTranslatef(25.0,0.0,37.0);
 glScalef(0.8,0.8,0.8);
 DrawScene();
 glPopMatrix();

 glPushMatrix();
 glTranslatef(-1.0,0.0,0.0);
 drawGate();
 glTranslatef(7.5,0.0,0.0);
 drawBoard(); 

 glPopMatrix();
 drawWall(2.5);

 glPushMatrix();
 glTranslatef(13.0,0.0,0.0);
 drawWall(2.5);
 glPopMatrix();
 
glPushMatrix();

      glTranslatef(15.0,0.0,0.0);
   
   glPushMatrix();
      glRotatef(270.0,0.0,1.0,0.0);
   drawWall(7.5);
   glPopMatrix();

     glTranslatef(0.0,0.0,8.0);
   glPushMatrix();
     glRotatef(270.0,0.0,0.1,0.0);
      drawWall(7.5);
   glPopMatrix();

     glTranslatef(0.0,0.0,8.0);
   glPushMatrix();
      glRotatef(225.0,0.0,1.0,0.0);
   drawWall(6.5);
   glPopMatrix();

     
      glTranslatef(-5.0,0.0,5.0);
   glPushMatrix();
      glRotatef(180.0,0.0,1.0,0.0);
   drawWall(4.5);
   
   glPopMatrix();
      
     glTranslatef(-5.0,0.0,0.0);
   
   glPushMatrix();
      glRotatef(135.0,0.0,1.0,0.0);
   drawWall(6.5);
   
   glPopMatrix();

   glTranslatef(-5.0,0.0,-5.0);
   
   glPushMatrix();
      glRotatef(90.0,0.0,1.0,0.0);
   drawWall(7.5);
   
   glPopMatrix();

   glTranslatef(0.0,0.0,-8.0);
    glPushMatrix();
      glRotatef(90.0,0.0,1.0,0.0);
   drawWall(7.5);
   
   glPopMatrix();

   glPopMatrix();

   

   // draw the tower
   glPushMatrix();
   glPushMatrix();
   glScalef(1.2,1.2,1.2);
   drawTower();
   glPopMatrix();
   
   glTranslatef(15.0,0.0,0.0);
   glPushMatrix();
   glScalef(1.2,1.2,1.2);
   drawTower();
   glPopMatrix();

   glTranslatef(0.0,0.0,8.0);
   drawTower();

   glTranslatef(0.0,0.0,8.0);
   drawTower();

   glTranslatef(-5.0,0.0,5.0);
   drawTower();

   glTranslatef(-5.0,0.0,0.0);
   drawTower();

   glTranslatef(-5.0,0.0,-5.0);
   drawTower();

   glTranslatef(0.0,0.0,-8.0);
   drawTower();
   glPopMatrix();

   //draw the house
   glPushMatrix();
   glTranslatef(3.0,0.0,9.0);
   drawhouse();
   glPopMatrix();
   
   //draw the house tower
   glPushMatrix();
   glTranslatef(3.0,0.0,9.0);
   drawhousetower();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(12.0,0.0,9.0);
   drawhousetower();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(3.0,0.0,14.0);
   drawhousetower();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(12.0,0.0,14.0);
   drawhousetower();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(7.5,0.0,15.0);
   glScalef(1.5,1.5,1.5);
   drawhousetower();
   glPopMatrix();

   //draw the flag
   glPushMatrix();
   glTranslatef(7.5,0.0,15.0);
   drawflag();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(-15.0,0.0,-15.0);
   drawGround();
   glPopMatrix();

   

    //glPushMatrix();
    //glTranslatef(30,1.0,-5.0);
   glColor3f(1,1,1);
    RenderScene();
    //glPopMatrix();

    objID = RetrieveObjectID(mouseX, mouseY);
    switch(objID){
      case Broom: 
        click_broom =0 ;
      //select the broomstick begin to fly
        Ex = 0.4;
        Ey = 3.0;
        Ez = -2.7;
        zh = 150;
        theta = -3;
      break;


      case Ball: 
      //select the ball, ball begin to fly
        click_ball =0;
    
      break;

      case Bat:
      //select the bat, begin to hit the ball
        click_bat =0;
      break;
      default:
      break;

    }

    if(click_ball ==0){
      if(hit == 1){
      float Position_ball[]  = {5*Cos(zh_l)+2.0,4.0,5*Sin(zh_l),1.0};
      
      glPushMatrix();

      glTranslatef(30,4.0,-5.0);
      //glColor3f(0.543,0.270,0.074);
      ball_brown(Position_ball[0],Position_ball[1],Position_ball[2],0.2);
      //glColor3f(1,1,1);
      glPopMatrix();
    }
      else{
        //start = clock();
        float Position_ball[]  = {8*t_bat*2*Cos(zh_bat)/(8*Sin(zh_bat)*Cos(zh)+3.0)-3.0,8*t_bat*2*Cos(zh_bat)/(8*Sin(zh_bat)*Sin(zh)-4.0)+4.0,t_bat*2};
      
      glPushMatrix();

      glTranslatef(30,4.0,-5.0);
      //glColor3f(0.543,0.270,0.074);
      ball_brown(Position_ball[0],Position_ball[1],Position_ball[2],0.2);
      //glColor3f(1,1,1);
      glPopMatrix();
      }
    }

    if(click_bat ==0){
      //float Position_ball[]  = {5*Cos(zh_l)+2.0,4.0,5*Sin(zh_l)b,1.0};
      
      glPushMatrix();
      //glRotatef(90,0,0,1);
      glTranslatef(26.6,8.0,-5.3);
      glRotatef(zh_bat,0,0,1);
      drawBat();
      glPopMatrix();
    }

   
   
   //  Draw axes - no lighting from here on
   glDisable(GL_LIGHTING);
   glColor3f(1,1,1);


   if (axes)
   {
      glBegin(GL_LINES);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(len+20.0,0.0,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,len+20.0,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,0.0,len+20.0);
      glEnd();
      //  Label axes
      glRasterPos3d(len+20.0,0.0,0.0);
      Print("X");
      glRasterPos3d(0.0,len+20.0,0.0);
      Print("Y");
      glRasterPos3d(0.0,0.0,len+20.0);
      Print("Z");
   }

   //  Display parameters
   glWindowPos2i(5,5);
   if (mode_project!=2)
   Print("Angle=%d,%d  Dim=%.1f FOV=%d Projection=%s Light=%s",
     th,ph,dim,fov,mode?"Perpective":"Orthogonal",light?"On":"Off");
   else if (mode_project ==2)
   Print("Angle=%d,%d  Dim=%.1f FOV=%d Projection=%s Light=%s",
     th,ph,dim,fov,"First person navigaiton",light?"On":"Off");

   if (light)
   {
      glWindowPos2i(5,45);
      Print("Model=%s LocalViewer=%s Distance=%d Elevation=%.1f",smooth?"Smooth":"Flat",local?"On":"Off",distance,ylight);
      glWindowPos2i(5,25);
      Print("Ambient=%d  Diffuse=%d Specular=%d Emission=%d Shininess=%.1f",ambient,diffuse,specular,emission,shiny);
   }

   //  Render the scene and make it visible
   ErrCheck("display");
   glFlush();
   glutSwapBuffers();
   glPopMatrix();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void idle()
{

   clock_t end = clock();
   t_bat = (end - start )/1000000.0;
   //  Elapsed time in seconds
   double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
   zh_l= fmod(90*t,360.0);
   float dtime = 0.004f;  //if you want to be exact, you would have to replace this by the real time passed since the last frame (and probably divide it by a certain number)
  g_timePassedSinceStart += dtime;

  if (g_timePassedSinceStart > 1.7f)
  {
    g_bExcitersInUse = false;  //stop the exciters
}


  UpdateScene(false,dtime,g_timePassedSinceStart);
 // Display();
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();

}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
   //  Right arrow key - increase angle by 5 degrees
   if (key == GLUT_KEY_RIGHT && mode_project!=2)
      th += 5;
   //  Left arrow key - decrease angle by 5 degrees
   else if (key == GLUT_KEY_LEFT && mode_project!=2)
      th -= 5;
   //  Up arrow key - increase elevation by 5 degrees
   else if (key == GLUT_KEY_UP && mode_project!=2)
      ph += 5;
   //  Down arrow key - decrease elevation by 5 degrees
   else if (key == GLUT_KEY_DOWN && mode_project!=2)
      ph -= 5;
   //  PageUp key - increase dim
   else if (key == GLUT_KEY_F5)
      dim += 0.1;
   //  PageDown key - decrease dim
   else if (key == GLUT_KEY_F6 && dim>1)
      dim -= 0.1;
   //  Smooth color model
   else if (key == GLUT_KEY_F1)
      smooth = 1-smooth;
   //  Local Viewer
   else if (key == GLUT_KEY_F2)
      local = 1-local;
   else if (key == GLUT_KEY_F3)
      distance = (distance==1) ? 5 : 1;
   else if(key == GLUT_KEY_LEFT && mode_project ==2){
      zh -= 5.0;                        //turn left for 2 degrees
   }
   else if (key== GLUT_KEY_RIGHT && mode_project ==2){
      zh += 5.0;                       //turn right for 2 degrees            
   }
// forward and backward
   else if(key == GLUT_KEY_UP && mode_project ==2){
      Ez -= scale*((float)Cos(zh) * speed);
      Ex += scale*((float)Sin(zh) * speed);
   }
   else if(key == GLUT_KEY_DOWN && mode_project==2){
      Ez += scale*((float)Cos(zh) * speed);
      Ex -= scale*((float)Sin(zh) * speed);
   }
   else if(key == GLUT_KEY_F9 && mode_project==2)
      Ey += scale*speed;
    else if(key == GLUT_KEY_F10 && mode_project==2)
      Ey -= scale*speed;
    else if(key == GLUT_KEY_F7 && mode_project==2)
    zh_bat += 5;

  else if(key == GLUT_KEY_F8 && mode_project==2)
    zh_bat -= 5;
  else if(key == GLUT_KEY_F4 && mode_project==2){
    hit = 1-hit;
    start = clock();
  }


   //  Keep angles to +/-360 degrees
   th %= 360;
   ph %= 360;
   zh %= 360;
   //  Update projection
   Project(mode_project?fov:0,asp,dim);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
   //  Exit on ESC
   if (ch == 27)
      exit(0);
   //  Reset view angle
   else if (ch == 'r' && mode_project!=2){
      th = 160;
      ph = 35;
      zh = 0;
      fov = 90;}
   else if (ch == 'r' && mode_project ==2){
      Ex = -2*dim*Sin(th)*Cos(ph);
      Ey = 0;
      Ez = +2*dim*Cos(th)*Cos(ph);
      zh=0;
   }
   //  Toggle axes
   else if (ch == 'x' || ch == 'X')
      axes = 1-axes;
   //  Toggle lighting
   else if (ch == 'l' || ch == 'L')
      light = 1-light;
   //  Switch projection mode
    else if ('1'<=ch && ch<='2')
      mode_project = ch-'0';
   //  Toggle light movement
   else if (ch == 'm' || ch == 'M')
      move = 1-move;
    else if (ch == 'f' || ch == 'F'){
      firework = 1-firework;
      Ex = 10;
      Ez = -5;
      theta = 5;
      zh = -90;
      stick_light = firework;
    }

   //  Move light
   else if (ch == '<')
      zh_l += 1;
   else if (ch == '>')
      zh_l -= 1;
   //  Change field of view angle
   else if (ch == '-' && ch>1)
      fov--;
   else if (ch == '+' && ch<179)
      fov++;
   //  Light elevation
   else if (ch=='[')
      ylight -= 0.1;
   else if (ch==']')
      ylight += 0.1;
   //  Ambient level
   else if (ch=='a' && ambient>0)
      ambient -= 5;
   else if (ch=='A' && ambient<100)
      ambient += 5;
   //  Diffuse level
   else if (ch=='d' && diffuse>0)
      diffuse -= 5;
   else if (ch=='D' && diffuse<100)
      diffuse += 5;
   //  Specular level
   else if (ch=='s' && specular>0)
      specular -= 5;
   else if (ch=='S' && specular<100)
      specular += 5;
   //  Emission level
   else if (ch=='e' && emission>0)
      emission -= 5;
   else if (ch=='E' && emission<100)
      emission += 5;
   //  Shininess level
   else if (ch=='n' && shininess>-1)
      shininess -= 1;
   else if (ch=='N' && shininess<7)
      shininess += 1;


   else if (ch == 'u' && mode_project==2)
      theta+=1;
   else if (ch=='U' && mode_project==2)
      theta-=1;
   else if (ch == 'b' && mode_project==2)
      scale+=0.01;
   else if (ch == 'B' && mode_project==2){
      if(scale-0.01>0)scale-=0.01; 
    }
   else if (ch == 'k' || ch == 'K'){
      box = 1-box; 
      if(box == 0){
        ambient = 0;
      }else
        ambient = 30;
   }
   else if ((ch == 'q' || ch == 'Q' )&& mode_project==2){

      Ex = -8;
      Ez = -11;
      zh = 180;
      theta = 0;
      firework = 0;

    }
   else if ((ch == 'w' || ch == 'W') && mode_project==2){
      Ex = 1.9;
      Ez = 13.7;
      zh = 190;
      Ey = -0.1;
      theta = 0;
      firework = 0;
  }
  else if(ch == 'c' || ch == 'C')
    move_light = 1- move_light;
  else if(ch == 'i' || ch == 'I')
    stick_light = 1- stick_light;

  else if(ch == 't' || ch == 'T'){
    Ex = 0.7;
    Ez =-3.7;
    Ey =0.0;
    theta = 0;
    zh = 120;
    click_broom = 1;
    click_bat = 1;
    click_ball = 1;
    
  
      //  Translate shininess power to value (-1 => 0)
  }

   shiny = shininess<0 ? 0 : pow(2.0,shininess);
   //  Tell GLUT it is necessary to redisplay the scene
   zh %=360;
   //  Reproject
   Project(mode_project?fov:0,asp,dim);
   //  Animate if requested
   glutIdleFunc(move?idle:NULL);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

void mouseCB(int button, int state, int x, int y)
{
    //mouseX = x;
    //mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        mouseX = x;
    mouseY = y;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        mouseX = x;
    mouseY = y;
  }
}



/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
  width_win = width;
  height_win = height;
   //  Ratio of the width to the height of the window
   asp = (height>0) ? (double)width/height : 1;
   //  Set the viewport to the entire window
   glViewport(0,0, width,height);
   //  Set projection
   Project(mode_project?fov:0,asp,dim);
}

/*
 *  Start up GLUT and tell it what to do
 */
int main(int argc,char* argv[])
{
   Ex = -2*dim*Sin(0)*Cos(0);
      Ey = 0;
      Ez = +2*dim*Cos(0)*Cos(0); 
   //  Initialize GLUT
   glutInit(&argc,argv);
   //  Request double buffered, true color window with Z buffering at 600x600
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(400,400);
   glutCreateWindow("Assignment 7: Xu Han");
   //compute the vertices and indices
  CreatePool();
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(  3,   //3 components per vertex (x,y,z)
            GL_FLOAT,
            sizeof(SOscillator),
            Oscillators);
  glNormalPointer(  GL_FLOAT,
            sizeof(SOscillator),
            &Oscillators[0].nx);  //Pointer to the first color*/
  glPointSize(2.0);
  glClearColor(0.0,0.0,0.0,0.0);


  glFrontFace(GL_CCW);   //Tell OGL which orientation shall be the front face
  glShadeModel(GL_SMOOTH);

  //initialize generation of random numbers:
  srand((unsigned)time(NULL));
  //InitParticles();
   

   //  Set callbacks
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   glutIdleFunc(idle);
   glutMouseFunc(mouseCB);
   //  Pass control to GLUT so it can interact with the user
   //  Load textures
   InitParticles();
   texture[0] = LoadTexBMP("1.bmp");
   texture[1] = LoadTexBMP("2.bmp");
   texture[2] = LoadTexBMP("3.bmp");
   texture[3] = LoadTexBMP("4.bmp");
   texture[4] = LoadTexBMP("5.bmp");
   texture[5] = LoadTexBMP("6.bmp");
   texture[6] = LoadTexBMP("7.bmp");
   texture[7] = LoadTexBMP("8.bmp");
   texture[8] = LoadTexBMP("9.bmp");
   texture[9] = LoadTexBMP("10.bmp");
   sky[0] = LoadTexBMP("sky0.bmp");
   sky[1] = LoadTexBMP("sky1.bmp");
   fire = LoadTexBMP("fire.bmp");
   texture[10] = LoadTexBMP("tree.bmp");
   texture[11] = LoadTexBMP("chair.bmp");
   texture[12] = LoadTexBMP("teapot.bmp");

   g_iLastRenderTime=clock();
   g_iLastRenderTime1=clock();

   ErrCheck("init");
   glutMainLoop();
   return 0;
}
