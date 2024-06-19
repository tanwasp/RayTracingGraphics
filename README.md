# Ray Tracing Project

## Description

This project is a ray tracer implementation which simulates the behavior of light to render 3D scenes onto a 2D image. It supports multiple geometric shapes, light sources, and material properties. The ray tracer outputs them in PPM (Portable Pixmap) format. It uses no additional libraries and is written in plain C++.

## Sample Output Images
![An infinite sea of pink spheres on a plane](images_png\picture17.png)

![A room defined by six planes with two spheres sitting on the floor.](images_png\picture16.png)

![Two transparent spheres refracting light from a third sphere](images_png\picture13.png)

Full sample scene and output image descriptions below



## PPM File Format

### Overview
The Portable Pixmap (PPM) format is a simple, ASCII-based image format. The structure of a PPM file is as follows:

1. **Magic Number**: Identifies the file type. For PPM, this is "P3".
2. **Whitespace**: Spaces, tabs, carriage returns, and line feeds.
3. **Width**: The width of the image in ASCII decimal format.
4. **Whitespace**: More whitespace.
5. **Height**: The height of the image in ASCII decimal format.
6. **Whitespace**: More whitespace.
7. **Max Color Value**: Maximum value for a color component in ASCII decimal format.
8. **Whitespace**: More whitespace.
9. **Pixel Data**: Width * height pixels, each represented by three ASCII decimal values for red, green, and blue components, respectively.


## PPM File Format

### Raw Format Variant
- **Magic Number**: "P6"
- **Pixel Data**: Stored as plain bytes without whitespace in the pixel area.
- **Max Color Value**: Must be <= 255.

## Scene File Format

### Structure
The scene file defines the camera, light sources, spheres, and planes in the scene.

### Camera
- **Origin**: 0.0 0.0 0.0
- **Projection Plane Z-coordinate**: -2.0
- **Raster Coordinates**: Min X and Max X: -1.0 1.0, Min Y and Max Y: -1.0 1.0
- **Background Color**: RGB values: 0.0 0.0 0.0
- **Ambient Light Color**: RGB values: 0.8 0.8 0.8
- **Max Ray Tracing Depth**: 1
- **Resolution**: Width and Height: 500 500

### Light Sources
- **Number of Light Sources**: 2
- **Light Source Details**:
  - **Position**: x y z
  - **Shading**: r g b
  - **Attenuation**: c0 c1 c2
  - **Example**: -4.5 4.5 6.0 1.0 1.0 1.0 1.0 0.0 0.0

### Spheres
- **Number of Spheres**: 1
- **Sphere Details**:
  - **Geometry**: x y z r
  - **Ambient Coefficients**: kar kag kab
  - **Diffuse Coefficients**: kdr kdg kdb
  - **Specular Coefficients**: ksr ksg ksb
  - **Reflectivity Coefficients**: krr krg krb
  - **Transmissivity Coefficients**: ktr ktg ktb
  - **Specularity**: specularity
  - **Index of Refraction**: refractionIndex
  - **Flags**: Reflective: isReflective, Transparent: isTransparent
  - **Example**: 0.0 0.0 -7.50 1.0 0.00 0.3 0.3 0.0 0.3 0.3 0.3 0.3 0.3 0.0 0.0 0.0 0.0 0.0 0.0 20.0 1.0 0 0

### Planes
- **Number of Planes**: 5
- **Plane Details**:
  - **Geometry**: x y z r
  - **Ambient Coefficients**: kar kag kab
  - **Diffuse Coefficients**: kdr kdg kdb
  - **Specular Coefficients**: ksr ksg ksb
  - **Reflectivity Coefficients**: krr krg krb
  - **Transmissivity Coefficients**: ktr ktg ktb
  - **Specularity**: specularity
  - **Index of Refraction**: refractionIndex
  - **Flags**: Reflective: isReflective, Transparent: isTransparent
  - **Checkered Pattern Directions**: p1x p1y p1z, p2x p2y p2z
  - **Example**: 0.0 0.0 1.0 -7.50 0.00 0.3 0.3 0.0 0.3 0.3 0.3 0.3 0.3 0.0 0.0 0.0 0.0 0.0 0.0 20.0 1.0 0 0

### Example scene
```
0.0  0.0  0.0

-2.0

-1.0  1.0

-1.0  1.0

0.0  0.0  0.0

1.0  1.0  1.0

2

500 500

2

 5.0   7.5  -5.0   1.0  1.0  1.0   1.0  0.0  0.0

-5.0   7.5  -5.0   1.0  1.0  1.0   1.0  0.0  0.0

3

 0.0  -7.0  -25.0  10.0    0.3  0.0  0.3    0.35  0.0  0.35   0.2  0.2  0.2   
 0.0  0.0  0.0             0.0  0.0  0.0    20.0  1.0  0   0  

 2.0  4.0  -11.5  1.0      0.0  0.3  0.3    0.0  0.35  0.35   0.2  0.2  0.2   
 0.0  0.0  0.0             0.0  0.0  0.0    20.0  1.0   0   0  

 -2.0  4.0  -11.5  1.0    0.3  0.3  0.0     0.35  0.35  0.0   0.2  0.2  0.2   
  0.0  0.0  0.0           0.0  0.0  0.0     20.0  1.0   0   0  
  ```

## Usage





## Sample Output Images
### picture01.ppm
Hidden surface removal.
A green sphere in front of a red sphere.

### picture02.ppm
A green sphere in front of a red sphere.
Another green sphere in front of another red sphere.
Distortion from using a wide field of view.

### picture03.ppm
A magenta sphere under diffuse lighting.
Two light sources.

### picture04.ppm
A magenta sphere under diffuse and specular lighting.
Two light sources.

### picture05.ppm
Making a room with five infinite planes, partly occluding each other.

### picture06.ppm
Illustration of hard shadows.
Two spheres casting shadows on a third larger sphere.
Two light sources.

### picture07.gif
Illustration of soft shadows.
Two spheres casting shadows on a third larger sphere.
Each of the two light sources was split into 27 smaller lights distributed over a small cube.

### picture08.gif
Buggy picture resulting from failure to use test t>0.0+epsilon when doing intersections. Ray from sphere to light source often intersects sphere, causing parts of the sphere to shadow itself.

### picture09.ppm
Five planes making a room.
A sphere inside casts shadows on the corners.

### picture10.ppm
A sphere being reflected into two mirrors.
The mirrors are at 45 degree angles w.r.t. the z and x axes, and parallel to the y axis.
Two lights in front of the sphere: lighting its near side.
Two lights behind the sphere: lighting its far side.
Ray tracing recursion depth = 2.

### picture11.ppm
Two spheres being reflected on the surface of a third sphere.
Ray tracing recursion depth = 2.

### picture12.alt.ppm
A sphere with a mirrored surface on a checkered floor in a room defined by six planes.
Ray tracing recursion depth = 2.
The one light source is visible as a specular highlight on the sphere.

### picture13.ppm
Two transparent spheres refracting light from a third sphere.
Ray tracing recursion depth = 3.

### picture14.ppm
A transparent sphere refracting light from a plane with a checkerboard pattern.
Ray tracing recursion depth = 3.

### picture15.alt.ppm
A transparent sphere refracting light from five sides of a room constructed out of planes.
Ray tracing recursion depth = 3.
The checkerboard pattern on the walls was created by using a function that computes coefficients of ambient reflection that depend on the point of intersection between the ray and the surface. The procedure determines whether the intersection lies in a colored or gray square and returns the corresponding coefficients.

### picture16.ppm
A room defined by six planes with two spheres sitting on the floor.
The two opposite ends of the room are mirror walls.
Ray tracing recursion depth = 500.
An infinity of reflections!

## Code Structure

### Header Files
- `ray-tracing.h`: Contains class and function declarations.

### Main Components

#### Globals and Utility Functions
- Various global variables for scene dimensions, resolution, etc.
- Utility functions for clamping and maximum value calculations.

#### Color Class
Handles color operations:
- Constructors for default, stream input, and specific color values.
- Methods to write to stream, add, subtract, multiply colors, and calculate magnitude and difference.

#### Vec Class
Represents a 3D vector:
- Constructors for default, stream input, specific coordinates, and copying.
- Methods for normalization, magnitude, dot product, addition, subtraction, and indexing.

#### Vect Class
Represents a 4-dimensional vector used for homogenous coordinates:
- Constructors and methods for element access, subtraction, dot product, norm calculation, and initialization from `Vec`.

#### Matrix Class
Represents a 4x4 matrix:
- Constructors and methods for row access, matrix-vector multiplication.

#### Ray Class
Represents a ray with an origin and end point:
- Methods for calculating direction and position along the ray.

#### Figure Class
Base class for all geometric shapes:
- Material properties and methods for loading them.
- Virtual methods for intersection calculation and normal vector calculation.

#### Light Class
Represents a light source:
- Position, shading color, and attenuation coefficients.
- Methods for reading coefficients, calculating attenuation, and distance squared.

#### Sphere Class
Derived from `Figure`, represents a sphere:
- Center and radius.
- Methods for intersection calculation and normal vector calculation.

#### Plane Class
Derived from `Figure`, represents a plane:
- Normal vector and distance from origin.
- Methods for intersection calculation and normal vector calculation.

#### Quadric Class
Derived from `Figure`, represents a quadric surface:
- Geometry matrix and clipping planes.
- Methods for intersection calculation, normal vector calculation, and clipping.

### Ray Tracing Functions

#### Scene Parsing
- `openFile()`: Opens a file for reading.
- `parseCameraSettings()`: Reads camera settings.
- `parseDimensions()`: Reads scene dimensions.
- `readColor()`: Reads a color from the file.
- `calculatePixels()`: Calculates pixel width and height.
- `parseObjects()`: Reads lights and shapes from the file.
- `sceneParsing()`: Parses the entire scene.

#### Image Generation
- `createImageStorage()`: Initializes the image storage.
- `initializePixelValues()`: Sets initial pixel values to black.
- `beginImageGeneration()`: Starts the image generation process.

#### Output Functions
- `writeP3Header()`: Writes the PPM header to the output file.
- `openOutputFile()`: Opens the output file for writing.

#### Reflection and Refraction
- `calculateReflectionDirection()`: Calculates the reflection direction.
- `calculateReflectedColor()`: Calculates the color from a reflected ray.