/*  Assignment 6: Xu Han
 *  Lighting and Texture
 *
 *
 *  Key bindings:
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
 *  </>        Change the positions of light source mannually while the light movement is stopped
 *  p          Toggles ortogonal/perspective projection
 *  +/-        Change field of view of perspective
 *  u/U        Make the first person navigation up and down 
 *  b/B        Change the scale of the scene(b: larger; B: smaller)
 *  x          Toggle axes
 *  arrows     Change view angle(<-/->: rotate around Y axis; up/down: rotate around X axis) when in mode 0 and mode 1
               Change the position and view angle(<-/->, make the first person navigation look left and right; up and down, 
               make the first person navigation move forward and backward)
 *  F5/F6      Zoom in and out for orthogonal
 *  r          Reset view angle
 *  ESC        Exit
 */
#include "CSCIx229.h"

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
int mode_project =0;
int move=1;       //  Move light
int th=0;         //  Azimuth of view angle
int ph=0;         //  Elevation of view angle
int fov=55;       //  Field of view (for perspective)
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
int ambient   =  30;  // Ambient intensity (%)
int diffuse   = 100;  // Diffuse intensity (%)
int specular  =   0;  // Specular intensity (%)
int shininess =   0;  // Shininess (power of two)
float shiny   =   1;  // Shininess (value)
int zh_l        =  90;  // Light azimuth
float ylight  =   0;  // Elevation of light

unsigned int texture[10];  //texture names

//first person navigation
void FirstpersonNaviagtion(void)
{
// the point where item lies

//gluLookAt(Ex, Ey, Ez, s_atx, s_aty, s_atz,0.0, 1.0, 0.0);
      gluLookAt(Ex, Ey, Ez, Ex + scale*10*Sin(zh), Ey+scale*50*Sin(theta), Ez - scale*10*Cos(zh),0.0, 1.0, 0.0);
}

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
   int th,ph;
   float yellow[] = {1.0,1.0,0.0,1.0};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball
   glColor3f(1,1,1);
   glMaterialf(GL_FRONT,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
   glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
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

/* draw branches
 *   for forbidden forest
 */

void drawBranches(float branchLength){
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
      YtoLowerHeight = (Radius-1.0) / (0.0-branchLength) * (Radius-1.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(x,branchLength,z);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(x,branchLength,z);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
}

/* draw a tree
 *   for forbidden forest
 */
void drawTree(float treeHeight, int numberOfBranches){
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
      YtoLowerHeight = (Radius-1.0) / (0.0-treeHeight) * (Radius-1.0);

      for ( i = 0; i < NumOfEdges; i++)    //create a circle
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(x,treeHeight,z);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(x,treeHeight,z);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);

      // draw braches
}

void drawLake(){
  
}
void drawForest(int numberOfTrees){

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
   /*glBegin(GL_QUADS);
      //Create the lower part of the tower:
      int i=0;
      
      //y is constant when the height is same 
      YtoLowerHeight = (Radius-1.0) / (0.0-LowerHeight) * (Radius-1.0);

      for ( i = 0; i < NumOfEdges-1; i++)    //create a circle
      {  
         x = cos((float)i/(float)NumOfEdges * 360.0);
         z = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);
         //same x,z and NVect:
         glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);

         x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);
         //same x,z and NVect:
         glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);
      }
     x = cos((float)i/(float)NumOfEdges * 360.0);
      z = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);
      //same x,z and NVect:
      glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);
      x = cos(0.0/(float)NumOfEdges * 360.0);
      z = sin(0.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);
      //same x,z and NVect:
      glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);
      glColor3f(1.0,1.0,1.0);
      glEnd();*/
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

         /*x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoLowerHeight,z);
         glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);
         //same x,z and NVect:
         glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);*/
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*Radius*3,z*Radius*2); glVertex3f(x*Radius,0.0,z*Radius);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(x,LowerHeight,z);
      /*x = cos(1.0/(float)NumOfEdges * 360.0);
      z = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord3f(x,LowerHeight,z); glVertex3f(x,LowerHeight,z);
      //same x,z and NVect:
      glTexCoord3f(x*Radius,0.0,z*Radius); glVertex3f(x*Radius,0.0,z*Radius);*/
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
   
      /*glBegin(GL_QUADS);
      for ( i = 0; i < NumOfEdges; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord2f(0,LowerHeight/4); glVertex3f(x,LowerHeight,z);
         glTexCoord2f(1,HigherHeight/4); glVertex3f(x,HigherHeight,z);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord2f(1,HigherHeight/4); glVertex3f(x,HigherHeight,z);
         glTexCoord2f(0,LowerHeight/4); glVertex3f(x,LowerHeight,z);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
      glTexCoord2f(0,LowerHeight/4); glVertex3f(x,LowerHeight,z);
      glTexCoord2f(1,HigherHeight/4); glVertex3f(x,HigherHeight,z);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
     glTexCoord2f(1,HigherHeight/4);  glVertex3f(x,HigherHeight,z);
      glTexCoord2f(0,LowerHeight/4); glVertex3f(x,LowerHeight,z);
      glEnd();*/
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
 *  Draw the gate
 */
void drawGate(){
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);

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
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(3.5,0.0,-2.0);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(3.5,HigherHeight*1.2,-2.0);
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
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(4.5/4,0);glVertex3f(4.5,0.0,-4.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,-4.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(4.5,HigherHeight*1.2,-4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0); glVertex3f(6.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(6.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(6.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.5,0.0,-2.0);
   
   
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
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(7.5/4,HigherHeight/8);glVertex3f(7.5,HigherHeight/2,-2.0);
   glTexCoord2f(9.5/4,HigherHeight/8);glVertex3f(9.5,HigherHeight/2,-2.0);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(7.5/4,HigherHeight*1.2/4);glVertex3f(7.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-2.0);
   glTexCoord2f(9.5/4,0);glVertex3f(9.5,0.0,-2.0);
   glTexCoord2f(9.5/4,HigherHeight*1.2/4);glVertex3f(9.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-2.0);
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
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,-4.0);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-4.0);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,-4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(0,0);glVertex3f(12.5,0.0,-4.0);
   glTexCoord2f(0,1*1.5);glVertex3f(12.5,HigherHeight*1.2,-4.0);
   glTexCoord2f(1/1.7,1*1.5);glVertex3f(12.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(12.5,0.0,-2.0);
   
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);

   glBegin(GL_QUADS);
   glNormal3f(0.0,-1.0,0.0);
   glTexCoord2f(12.5/4,0);glVertex3f(12.5,0.0,-2.0);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,-2.0);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(12.5/4,HigherHeight*1.2/4);glVertex3f(12.5,HigherHeight*1.2,-2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(1.0,0.0,0.0);
   glTexCoord2f(1*1.5,0);glVertex3f(13.5,0.0,-2.0);
   glTexCoord2f(1*1.5,1*1.5);glVertex3f(13.5,HigherHeight*1.2,-2.0);
   glTexCoord2f(0,1*1.5);glVertex3f(13.5,HigherHeight*1.2,2.0);
   glTexCoord2f(0,0); glVertex3f(13.5,0.0,2.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(13.5/4,0);glVertex3f(13.5,0.0,2.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,2.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,2.0);
   glTexCoord2f(13.5/4,HigherHeight*1.2/4);glVertex3f(13.5,HigherHeight*1.2,2.0);
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
   glNormal3f(0.0,0.0,1.0);
   glTexCoord2f(10.5/4,0);glVertex3f(10.5,0.0,4.0);
   glTexCoord2f(6.5/4,0);glVertex3f(6.5,0.0,4.0);
   glTexCoord2f(6.5/4,HigherHeight*1.2/4);glVertex3f(6.5,HigherHeight*1.2,4.0);
   glTexCoord2f(10.5/4,HigherHeight*1.2/4);glVertex3f(10.5,HigherHeight*1.2,4.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   glTexCoord2f(0.0,0.0);glVertex3f(6.5,0.0,4.0);
   glTexCoord2f(0,1.5);glVertex3f(6.5,HigherHeight*1.2,4.0);
   glTexCoord2f(1/1.7,1.5);glVertex3f(6.5,HigherHeight*1.2,2.0);
   glTexCoord2f(1/1.7,0);glVertex3f(6.5,0.0,2.0);
   
   
   
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
   /*glBegin(GL_QUADS);
   int i=0;
   float x,z;
   // draw lowerpart
   
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = cos((float)i/(float)NumOfEdges * 360.0);
         z = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x*0.8,0.0,z*0.8); glVertex3f(x*0.8,0.0,z*0.8);
         glTexCoord3f(x*0.8,HouseHeight,z*0.8); glVertex3f(x*0.8,HouseHeight,z*0.8);

         x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x*0.8,HouseHeight,z*0.8);glVertex3f(x*0.8,HouseHeight,z*0.8);
         glTexCoord3f(x*0.8,0.0,z*0.8);glVertex3f(x*0.8,0.0,z*0.8);
      }
      x = cos((float)i/(float)NumOfEdges * 360.0);
      z = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
      glTexCoord3f(x*0.8,0.0,z*0.8);glVertex3f(x*0.8,0.0,z*0.8);
      glTexCoord3f(x*0.8,HouseHeight,z*0.8);glVertex3f(x*0.8,HouseHeight,z*0.8);
      x = cos(1.0/(float)NumOfEdges * 360.0);
      z = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
      glTexCoord3f(x*0.8,HouseHeight,z*0.8);glVertex3f(x*0.8,HouseHeight,z*0.8);
      glTexCoord3f(x*0.8,0.0,z*0.8);glVertex3f(x*0.8,0.0,z*0.8);
   glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);*/
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
      /*glBegin(GL_QUADS);
      
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = cos((float)i/(float)NumOfEdges * 360.0);
         z = sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x,HouseHeight,z); glVertex3f(x,HouseHeight,z);
         glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);

         x = cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,0.0,z);
         glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);
         glTexCoord3f(x,HouseHeight,z); glVertex3f(x,HouseHeight,z);
      }
      x = cos((float)i/(float)NumOfEdges * 360.0);
      z = sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
      glTexCoord3f(x,HouseHeight,z);  glVertex3f(x,HouseHeight,z);
      glTexCoord3f(x,HouseHeight+2.0,z);  glVertex3f(x,HouseHeight+2.0,z);
      x = cos(1.0/(float)NumOfEdges * 360.0);
      z = sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,0.0,z);
      glTexCoord3f(x,HouseHeight+2.0,z);  glVertex3f(x,HouseHeight+2.0,z);
      glTexCoord3f(x,HouseHeight,z); glVertex3f(x,HouseHeight,z);*/
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
      /*glBegin(GL_TRIANGLES);
      
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoPoint,z);
         glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);
         glTexCoord3f(0.0,HouseHeight+4.0,0.0); glVertex3f(0.0,HouseHeight+4.0,0.0);
         

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(x,YtoPoint,z);
         glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);
         //glTexCoord3f(x,HouseHeight+2.0,z); 
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoPoint,z);
      glTexCoord3f(0.0,HouseHeight+4.0,0.0); glVertex3f(0.0,HouseHeight+4.0,0.0);
      glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoPoint,z);
      glTexCoord3f(x,HouseHeight+2.0,z); glVertex3f(x,HouseHeight+2.0,z);
      //glTexCoord3f(x,HouseHeight+2.0,z); */
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

   /*glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(1.0,1.0,0.0);
   glBindTexture(GL_TEXTURE_2D,texture[8]);

   glBegin(GL_TRIANGLES);
   glNormal3f(0.0,0.0,-1.0);

   glTexCoord3f(3.5/2,HouseHeight/2,0.0); glVertex3f(3.5,HouseHeight,0.0);
   glTexCoord3f(5.5/2,HouseHeight/2,0.0); glVertex3f(5.5,HouseHeight,0.0);
   glTexCoord3f(4.5/2,(HouseCenterHeight+1.0)/2,0.0); glVertex3f(4.5,HouseCenterHeight+1.0,0.0);
   glColor3f(1.0,1.0,1.0);
   glEnd();

    glBegin(GL_TRIANGLES);
   glNormal3f(0.0,0.0,1.0);
   glColor3f(1.0,1.0,0.0);
   glTexCoord3f(3.5/2,HouseHeight/2,0.01/2); glVertex3f(3.5,HouseHeight,0.01);
   glTexCoord3f(5.5/2,HouseHeight/2,0.01/2); glVertex3f(5.5,HouseHeight,0.01);
   glTexCoord3f(4.5/2,(HouseCenterHeight+1.0)/2,0.01/2); glVertex3f(4.5,HouseCenterHeight+1.0,0.01);
   glColor3f(1.0,1.0,1.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);*/
   
}

//draw a flag
void drawflag(){
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

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
   
   glTexCoord2f(0,0); glVertex3f(0.0,-0.05,0.0);
   glTexCoord2f(25/4,0); glVertex3f(50,-0.05,0.0);
   glTexCoord2f(25/4,20/4); glVertex3f(50,-0.05,50);
   glTexCoord2f(0,20/4); glVertex3f(0.0,-0.05,50);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,-1.0,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,-0.1,0.0);
   glTexCoord2f(25/4,0); glVertex3f(50,-0.1,0.0);
   glTexCoord2f(25/4,20/4); glVertex3f(50,-0.1,50);
   glTexCoord2f(0,20/4); glVertex3f(0.0,-0.1,50);
   //glColor3f(1.0,1.0,1.0);
   glEnd();
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

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
   else if (mode_project==0)
   {
         glRotated(ph,1,0,0);
         glRotated(th,0,1,0);
   }
   else if (mode_project==2){
      // set parameters for glulookat()
              
      FirstpersonNaviagtion();
   }


   //  Flat or smooth shading
   glShadeModel(smooth ? GL_SMOOTH : GL_FLAT);

   //  Light switch
   if (light)
   {
        //  Translate intensity to color vectors
        float Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0};
        float Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0};
        float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
        //  Light position
        float Position[]  = {distance*Cos(zh_l),ylight,distance*Sin(zh_l),1.0};
        //  Draw light position as ball (still no lighting here)
        glColor3f(1,1,1);
        ball(Position[0],Position[1],Position[2] , 0.1);
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
        glLightfv(GL_LIGHT0,GL_POSITION,Position);
   }
   else
     glDisable(GL_LIGHTING);
   //  Set view angle
   glPushMatrix();
   glTranslatef(0.0,-0.1,0.0);


   //  Draw the stick
  if(mode_project== 2){
  glPushMatrix();
   glTranslatef(Ex,Ey-1.0,Ez);
   glScalef(scale,scale,scale);
   drawWand(10*Sin(zh), 50*Sin(theta+10), -10*Cos(zh));
   glPopMatrix();
}

   //  Draw scene

 glScalef(scale,scale,scale);
 glPushMatrix();
 glTranslatef(-1.0,0.0,0.0);
 drawGate();
 glTranslatef(7.5,0.0,0.0);
 drawBoard(); 

 glPopMatrix();
 drawWall(3.5);

 glPushMatrix();
 glTranslatef(12.0,0.0,0.0);
 drawWall(3.5);
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
      Print("Ambient=%d  Diffuse=%d Specular=%d Emission=%d Shininess=%.0f",ambient,diffuse,specular,emission,shiny);
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
   //  Elapsed time in seconds
   double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
   zh_l= fmod(90*t,360.0);
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
      th = -200;
      ph = 30;
      zh = 0;}
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
    else if ('0'<=ch && ch<='2')
      mode_project = ch-'0';
   //  Toggle light movement
   else if (ch == 'm' || ch == 'M')
      move = 1-move;
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
      //  Translate shininess power to value (-1 => 0)
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

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
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
   glutCreateWindow("Assignment 6: Xu Han");
   //  Set callbacks
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   glutIdleFunc(idle);
   //  Pass control to GLUT so it can interact with the user
   //  Load textures
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
   ErrCheck("init");
   glutMainLoop();
   return 0;
}