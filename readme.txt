picture01.ppm 
 - Hidden surface removal.  
 - A green sphere in front of a red sphere.  

picture02.ppm
 - A green sphere in front of a red sphere.  
 - Another green sphere in front of another red sphere.  
 - Distortion from using a wide field of view. 

picture03.ppm
 - A magenta sphere under diffuse lighting. 
 - Two light sources.  

picture04.ppm
 - A magenta sphere under diffuse and specular lighting. 
 - Two light sources.  

picture05.ppm
 - Making a room with five infinite planes, partly occluding each other. 

picture06.ppm
 - Illustration of hard shadows. 
 - Two spheres casting shadows on a third larger sphere. 
 - Two light sources. 

picture07.gif
 - Illustration of soft shadows. 
 - Two spheres casting shadows on a third larger sphere. 
 - Each of the two light sources was splint into 27 smaller lights
   distributed over a small cube.  

picture08.gif
 - Buggy picture resulting from failure to use test t>0.0+epsilon when
   doing intersections. Ray from sphere to light source often intersects
   sphere, with effect that parts of the sphere shadow itself. 

picture09.ppm
 - Five planes making a room. 
 - A sphere inside casts shadows on the corners. 

picture10.ppm
 - A sphere being reflected into two mirrors. 
 - The mirrors are at 45 degree angles w.r.t. the z and x axies.  
   and parallel to the y axis. 
 - Two lights in front of the sphere: Lighting its near side. 
 - Two lights in behind of the sphere: Lighting its far side. 
 - Ray tracing recursion depth = 2. 

picture11.ppm
 - Two spheres being reflected on the surface of a third sphere. 
 - Ray tracing recursion depth = 2. 

picture12.alt.ppm
  - A sphere with a mirrored surface on a checkered floor in a room
    defined by six planes. 
  - Ray tracing recursion depth = 2. 
  - The one light source is visible as a specular highlight on the sphere.  

picture13.ppm
 - Two transparent spheres refracting light from a third sphere. 
 - Ray tracing recursion depth = 3. 

picture14.ppm
 - A transparent sphere refracting light from a plane with a
   checkerboard pattern.  
 - Ray tracing recursion depth = 3. 

picture15.alt.ppm
 - A transparent sphere refracting light from five sides of a room 
   constructed out of planes. 
 - Ray tracing recursion depth = 3. 
 - The checkerboard pattern on the walls was created by using a function 
   that computes coefficients of amibient reflection that depend on the
   point of intersection between the ray and the surface. The 
   procedure determines whether the intersection lies in a colored
   or gray square and returns the corresponding coefficients. 

picture16.ppm
  - A room defined by six planes with two spheres sitting on the floor. 
  - The two opposite ends of the room are mirror walls.  
  - Ray tracing recursion depth = 500. 
  - An infinity of reflections! 

picture17.ppm
  - A room defined by six planes with two spheres sitting on the floor. 
  - The two opposite ends of the room are mirror walls.  
  - The left and right sides of the room are mirror walls.  
  - Ray tracing recursion depth = 500. 
  - An infinity of reflections! 

