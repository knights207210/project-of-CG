/*  Assignment 7: Xu Han
 *  Project
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
#include <vector>
#define PI 3.14159265359
#define NUM_X_OSCILLATORS   150
#define NUM_Z_OSCILLATORS   150
#define NUM_OSCILLATORS     NUM_X_OSCILLATORS*NUM_Z_OSCILLATORS
#define OSCILLATOR_DISTANCE   0.05

#define OSCILLATOR_WEIGHT       0.0002
#define SQR(a) (a*a)

#define BILLBOARDING_NONE             0
#define BILLBOARDING_PERPTOVIEWDIR          1  //align particles perpendicular to view direction
#define BILLBOARDING_PERPTOVIEWDIR_BUTVERTICAL    2  //like PERPToViewDir, but Particles are vertically aligned
#define RANDOM_FLOAT (((float)rand())/RAND_MAX)

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
int firework=1;

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
int ambient   =  30;  // Ambient intensity (%)
int diffuse   = 100;  // Diffuse intensity (%)
int specular  =   0;  // Specular intensity (%)
int shininess =   0;  // Shininess (power of two)
float shiny   =   1;  // Shininess (value)
int zh_l        =  90;  // Light azimuth
float ylight  =   0;  // Elevation of light
int numberOfTrees = 30;
int obj;

unsigned int texture[10];  //texture names

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
/*float operator* (SF3dVector v, SF3dVector u)  //Scalar product
{
  return v.x*u.x+v.y*u.y+v.z*u.z;
}*/



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

  bool      LoadTextureFromFile(char * Filename);

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

/*bool CCCParticleSystem::LoadTextureFromFile(char * Filename)
{
  //Create the texture pointer:
  m_Texture = new COGLTexture;

  if (m_Texture == NULL) return false;

  if (!m_Texture->LoadFromTGA(Filename,NULL,true)) return false;  //pass NULL as 2. param (only required if you want to combine rgb and alpha maps)

  m_Texture->SetActive();
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  
  
  m_bUseTexture = true;


  return true;

}*/


void CCCParticleSystem::Render()
{
  //the calling method must switch on texturing!

  if (m_bUseTexture)
  {
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
clock_t g_iLastRenderTime;

void InitParticles()
{

//INIT SYSTEM 2 (POINTS, FIREWORK)
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
  g_ParticleSystem1.SetEmitter(Ex+10*Sin(zh),Ey-1.0+50*Sin(theta+10),Ez-10*Cos(zh),
                0.02f,0.0f,0.02f);
  g_ParticleSystem1.SetAcceleration(F3dVector(0.0f,-1.0f,0.0f),0.83f,1.4f);
  g_ParticleSystem1.SetSizeValues(3.0f,3.0f,4.0f,4.0f);
  g_ParticleSystem1.m_fMaxEmitSpeed = 0.82f;
  g_ParticleSystem1.m_fMinEmitSpeed = 1.3f;
  g_ParticleSystem1.SetEmissionDirection(-1.0f,2.0f,0.0f,
                    0.5f,0.5f,0.5f);

  g_ParticleSystem1.m_bParticlesLeaveSystem = true;

}



///////////////////////////////////////////////////


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


void UpdateScene(bool bEndIsFree, float deltaTime, float time)
{
//********
// Here we do the physical calculations: 
// The oscillators are moved according to their neighbors.
// The parameter bEndIsFree indicates, whether the oscillators in the edges can move or not.
// The new position may be assigned not before all calculations are done!

// PLEASE NOTE: THESE ARE APPROXIMATIONS AND I KNOW THIS! (but is looks good, doesn't it?)

  //if we use two loops, it is a bit easier to understand what I do here.
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
      //Take the direction vectors 1.) from the left to the right neighbor 
      // and 2.) from the upper to the lower neighbor.
      //The vector orthogonal to these 

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
// the point where item lies

//gluLookAt(Ex, Ey, Ez, s_atx, s_aty, s_atz,0.0, 1.0, 0.0);
      gluLookAt(Ex, Ey, Ez, Ex + 0.2*scale*10*Sin(zh), Ey+0.2*scale*50*Sin(theta), Ez - 0.2*scale*10*Cos(zh),0.0, 1.0, 0.0);
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

/* 
 *  Draw sky box
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
 *  Draw sky box
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

void drawCourtTower(){
    
  float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  

  //draw lower part
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_QUADS);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,0.0,0.0);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,0.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,12.0,0.0);
   
   
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(2.0,0.0,0.0);
   glTexCoord2f(4.5/4,0);glVertex3f(2.0,0.0,2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,0.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
   glTexCoord2f(3.5/4,0);glVertex3f(2.0,0.0,2.0);
   glTexCoord2f(4.5/4,0);glVertex3f(0.0,0.0,2.0);
   glTexCoord2f(4.5/4,HigherHeight*1.2/4);glVertex3f(0.0,12.0,2.0);
   glTexCoord2f(3.5/4,HigherHeight*1.2/4);glVertex3f(2.0,12.0,2.0);
   glEnd();

   glBegin(GL_QUADS);
   glNormal3f(0.0,0.0,-1.0);
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
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,12.0,2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(2.0,14.0,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(2.0,14.0,0.4);
   glTexCoord2f(1*1.5,0); glVertex3f(2.0,12.0,2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,14.0,0.4);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,14.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(2.0,15.5,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(2.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,15.5,0.4);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,16.5,2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(2.0,15.5,0.0);
   glTexCoord2f(0,0); glVertex3f(2.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,12.0,2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,12.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(0.0,14.0,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,14.0,0.4);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,12.0,2.0);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,14.0,0.4);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,14.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(0.0,15.5,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,15.5,0.4);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,16.5,2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(0.0,16.5,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,15.5,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,15.5,0.4);
   glEnd();

   glBegin(GL_POLYGON);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,12.0,0.0);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,12.0,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(1*1.5,0); glVertex3f(0.0,16.5,0.0);
   glTexCoord2f(0,0); glVertex3f(0.0,15.5,0.4);
   glColor3f(1.0,1.0,1.0);
   glEnd();

   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   //draw the pointer
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);


   glBegin(GL_TRIANGLES);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,16.5,0.0);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(1.0,18.5,1.0);
   glEnd();

   glBegin(GL_TRIANGLES);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,16.5,0.0);
   glTexCoord2f(0,1*1.5); glVertex3f(2.0,16.5,2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(1.0,18.5,1.0);
   glEnd();

   glBegin(GL_TRIANGLES);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(2.0,16.5,2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,16.5,2.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(1.0,18.5,1.0);
   glEnd();

  glBegin(GL_TRIANGLES);
   glNormal3f(-1.0,0.0,0.0);
   
   glTexCoord2f(0,0); glVertex3f(0.0,16.5,2.0);
   glTexCoord2f(0,1*1.5); glVertex3f(0.0,16.5,0.0);
   glTexCoord2f(1*1.5,1*1.5); glVertex3f(1.0,18.5,1.0);
   glEnd();
   glColor3f(1.0,1.0,1.0);
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);

   //pointer
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[2]);
   
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
         glTexCoord2f(x*0.1*3,z*0.1*2); glVertex3f((x)*0.1+1,18.5,(z)*0.1+1);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(1.0,20.0,1.0);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*0.1*3,z*0.1*2); glVertex3f((x)*0.1+1,18.5,(z)*0.1+1);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(1.0,20.0,1.0);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
}

/*
r = torus ring radius
c = torus tube radius
rSeg, cSeg = number of segments/detail
*/

void drawTorus(double r, double c,int rSeg, int cSeg)
{
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  glFrontFace(GL_CW);

  glPushMatrix();
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture[1]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  //const double PI = 3.1415926535897932384626433832795;
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
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
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
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
  
 glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,mode?GL_REPLACE:GL_MODULATE);
   glColor3f(0.867,0.719,0.527);
   glBindTexture(GL_TEXTURE_2D,texture[3]);
   glScalef(0.5,0.0,1.0);
   //draw the ground
      glBegin(GL_TRIANGLES);
      int i;
      float x,z;
      for ( i = 0; i < NumOfEdges-1; i++)
      {  
         x = Cos((float)i/(float)NumOfEdges * 360.0);
         z = Sin((float)i/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(0.0,0.0,0.0);
         glVertex3f(x*20,0.0,z*20);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*20,0.0,z*20);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(x*20,0.0,z*20);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*20,0.0,z*20);
      glColor3f(1.0,1.0,1.0);
      glEnd();
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
         glVertex3f(x*20,0.0,z*20);
         glVertex3f(x*20,5.0,z*20);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*20,5.0,z*20);
         glVertex3f(x*20,0.0,z*20);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*20,0.0,z*20);
      glVertex3f(x*20,5.0,z*20);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*20,5.0,z*20);
      glVertex3f(x*20,0.0,z*20);
      glColor3f(1.0,1.0,1.0);
      glEnd();
      glPopMatrix();
      glDisable(GL_TEXTURE_2D);
      //draw court tower
      glPushMatrix();
      glTranslatef(10.0,0.0,0.0);
      glRotatef(270,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-10.0,0.0,1.5);
      glRotatef(90,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(1.5,0.0,20.0);
      glRotatef(180,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(0.0,0.0,-20.0);
      glRotatef(0,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(7.0,0.0,15.0);
      glRotatef(210,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-7.0,0.0,-14.5);
      glRotatef(20,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(6.3,0.0,-15.0);
      glRotatef(-40,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-5.9,0.0,15.8);
      glRotatef(140,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(9.5,0.0,7.5);
      glRotatef(220,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-9.5,0.0,-7.5);
      glRotatef(50,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(8.8,0.0,-7.5);
      glRotatef(-50,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
      glPopMatrix();

      glPushMatrix();
      glTranslatef(-8.5,0.0,8.4);
      glRotatef(-220,0,1,0);
      glScalef(0.8,0.8,0.8);
      drawCourtTower();
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
      glPopMatrix();

      //draw marks
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
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*16,0.01,z*16);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*16,0.01,z*16);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*16,0.01,z*16);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*16,0.01,z*16);
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
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*3,0.01,z*3);

         x = Cos((float)(i+1)/(float)NumOfEdges * 360.0);
         z = Sin((float)(i+1)/(float)NumOfEdges * 360.0);
         glNormal3f(0.0,-1.0,0.0);
         glVertex3f(x*3,0.01,z*3);
      }
      x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*3,0.01,z*3);
      x = Cos(1.0/(float)NumOfEdges * 360.0);
      z = Sin(1.0/(float)NumOfEdges * 360.0);
      glNormal3f(0.0,-1.0,0.0);
      glVertex3f(x*3,0.01,z*3);
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
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[2]);
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
         glTexCoord2f(x*0.5*3,z*0.5*2); glVertex3f(x*0.5,0.0,z*0.5);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(pointx,pointy,pointz);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*0.5*3,z*0.5*2); glVertex3f(x*0.5,0.0,z*0.5);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(pointx,pointy,pointz);
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
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[2]);
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
         glTexCoord2f(x*r1*3,z*r1*2); glVertex3f(x*r1,0.0,z*r1);
         //same x,z and NVect:
         glTexCoord2f(x*r2*3,z*r2*2); glVertex3f(x*r2,treeHeight,z*r2);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*r1*3,z*r1*2); glVertex3f(x*r1,0.0,z*r1);
      //same x,z and NVect:
      glTexCoord2f(x*r2*3,z*r2*2); glVertex3f(x*r2,treeHeight,z*r2);
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
   glColor3f(0.683,0.930,0.930);
   glBindTexture(GL_TEXTURE_2D,texture[2]);
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
         glTexCoord2f(x*r2*3,z*r2*2); glVertex3f(x*r2,treeHeight,z*r2);
         //same x,z and NVect:
         glTexCoord2f(x*3,z*2); glVertex3f(0,treeHeight+treeHeight*0.1,0);
      }
     x = Cos((float)i/(float)NumOfEdges * 360.0);
      z = Sin((float)i/(float)NumOfEdges * 360.0);
      glNormal3f(x,YtoLowerHeight,z);
      glTexCoord2f(x*r*3,z*r*2); glVertex3f(x*r2,treeHeight,z*r2);
      //same x,z and NVect:
      glTexCoord2f(x*3,z*2); glVertex3f(0,treeHeight+treeHeight*0.1,0);
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
   else if (mode_project==0)
   {
         glRotated(ph,1,0,0);
         glRotated(th,0,1,0);
   }
   else if (mode_project==2){
      // set parameters for glulookat()
    fov = 40;
              
      FirstpersonNaviagtion();
   }


   //  Flat or smooth shading
   glShadeModel(smooth ? GL_SMOOTH : GL_FLAT);

   if(box) 
    Sky(3.5*dim);
   else
    Sky1(3.5*dim);
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
   glTranslatef(0.0,-0.5,0.0);


   //  Draw the stick
  if(mode_project== 2){
  glPushMatrix();
   glTranslatef(Ex,Ey-1.0,Ez);
   glScalef(0.7*scale,0.7*scale,0.7*scale);
   drawWand(10*Sin(zh), 50*Sin(theta+10), -10*Cos(zh));
   glPopMatrix();


}



//partical engine
if(firework){
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
}
   //  Draw scene

 glScalef(2*scale,2*scale,2*scale);
 glTranslatef(-25,0,0);
 

 /*glPushMatrix();
 glScalef(0.45,0.45,0.45);
   glTranslatef(48.0,0.0,0.0);
   drawBridge();
   glPopMatrix();*/

 /*glPushMatrix();
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
   glPopMatrix();*/

   



   
   
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
   float dtime = 0.004f;  //if you want to be exact, you would have to replace this by the real time passed since the last frame (and probably divide it by a certain number)
  g_timePassedSinceStart += dtime;

  if (g_timePassedSinceStart > 1.7f)
  {
    g_bExcitersInUse = false;  //stop the exciters
  }
/*  //ENABLE THE FOLLOWING LINES FOR A RAIN EFFECT
  int randomNumber = rand();
  if (randomNumber < NUM_OSCILLATORS)
  {
    Oscillators[randomNumber].y = -0.05;
  }
  */


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
    else if (ch == 'f' || ch == 'F')
      firework = 1-firework;
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

   g_iLastRenderTime=clock();

   ErrCheck("init");
   glutMainLoop();
   return 0;
}
