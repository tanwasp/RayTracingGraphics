#include <cstdlib>
#include <cmath>
#include <cassert>
#include <list>
#include <vector>

#include <iostream>

#include <cfloat>
#include <fstream>



using namespace std;

// importing header file

#include "ray-tracing.h"

// defining global variables and constants

double  pixH;
double pixW;


double minimumY, maximumY;
double minimumX,   maximumX;



double  pixStartX;
double pixStartY;

int ssRate = 1;

int horizontalRes, verticalRes;
double zCoor;


vector < vector<  Color > > image;

list<Figure*> shapeList;
list<Light*> lightList;


const double epsilon = 0.00000000001;





double camera1, camera2, camera3;


const double aaThreshold = 0.25;

int maxDepth;
Color backgroundColor;
Color ambient;






// Utility functions for clamping and getting maximum values



// Color Implementation

int clamp(int value, int minVal, int maxVal) { return (value < minVal) ? minVal : (value > maxVal ? maxVal : value); }
int getMax(int first, int second) { return (first > second) ? first : second; }

const double maximumTimes = 1000000;

Color::Color() : redComponent(0.0), greenComponent(0.0), blueComponent(0.0) {}

Color::Color(std::ifstream& input) {
    input >> redComponent >> greenComponent >> blueComponent;
}

const Color Default_black;


Color::Color(double red, double green, double blue) : redComponent(red), greenComponent(green), blueComponent(blue) {}

void Color::writeToStream(std::ostream& output) {
    const int colorScale = 255;
    int rOutput = clamp(static_cast<int>(redComponent * colorScale), 0, colorScale);
    int gOutput = clamp(static_cast<int>(greenComponent * colorScale), 0, colorScale);
    int bOutput = clamp(static_cast<int>(blueComponent * colorScale), 0, colorScale);
    output << rOutput << " " << gOutput << " " << bOutput << std::endl;
}

Color Color::add(const Color& other) const {
    return Color(redComponent + other.redComponent, greenComponent + other.greenComponent, blueComponent + other.blueComponent);
}

Color Color::subtract(const Color& other) const {
    return Color(redComponent - other.redComponent, greenComponent - other.greenComponent, blueComponent - other.blueComponent);
}

Color Color::multiply(const Color& other) const {
    return Color(redComponent * other.redComponent, greenComponent * other.greenComponent, blueComponent * other.blueComponent);
}

double Color::magnitude() const {
    return std::fabs(redComponent) + std::fabs(greenComponent) + std::fabs(blueComponent);
}

double Color::difference(const Color& other) const {
    return this->subtract(other).magnitude();
}

// Color Implementation End




// Vec Implementation




Vec::Vec() : components{ 0.0, 0.0, 0.0 } {}

Vec::Vec(std::ifstream& input) {
    input >> components[0] >> components[1] >> components[2];
}

Vec::Vec(double x, double y, double z) : components{ x, y, z } {}

Vec::Vec(const Vec& other) : components{ other.components[0], other.components[1], other.components[2] } {}

void Vec::normalize() {
    double len = magnitude();
    if (len == 0) throw std::runtime_error("Attempt to normalize a zero vector");
    components[0] /= len;
    components[1] /= len;
    components[2] /= len;
}

double Vec::magnitude() const {
    return std::sqrt(components[0] * components[0] + components[1] * components[1] + components[2] * components[2]);
}

double Vec::dotProduct(const Vec& other) const {
    return components[0] * other.components[0] + components[1] * other.components[1] + components[2] * other.components[2];
}

Vec Vec::add(const Vec& other) const {
    return Vec(components[0] + other.components[0], components[1] + other.components[1], components[2] + other.components[2]);
}

Vec Vec::subtract(const Vec& other) const {
    return Vec(components[0] - other.components[0], components[1] - other.components[1], components[2] - other.components[2]);
}

const double& Vec::operator[](int index) const {
    if (index < 0 || index >= 3) throw std::out_of_range("Index out of range");
    return components[index];
}

double& Vec::operator[](int index) {
    if (index < 0 || index >= 3) throw std::out_of_range("Index out of range");
    return components[index];
}

Vec::Vec(const Vect& v) : components{ v.getElement(0), v.getElement(1), v.getElement(2)} {}

Vec operator*(double scalar, const Vec& vector) {
    return Vec(scalar * vector.components[0], scalar * vector.components[1], scalar * vector.components[2]);
}


// Vec Implementation End


// Vect Implementation


Vect::Vect() : data(4, 0.0) {}  // Default constructor initializes a 4-dimensional zero vector.

Vect::Vect(const Vec& vec) {
    initializeFromVec(vec);
}

double Vect::getElement(int index) const {
    if (index < 0 || index >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[index];
}

void Vect::setElement(int index, double value) {
    if (index < 0 || index >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    data[index] = value;
}

Vect Vect::subtract(const Vect& other) const {
    Vect result;
    for (size_t i = 0; i < data.size(); ++i) {
        result.setElement(i, this->getElement(i) - other.getElement(i));
    }
    return result;
}

double Vect::dotProduct(const Vect& other) const {
    double result = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        result += this->getElement(i) * other.getElement(i);
    }
    return result;
}

double Vect::calculateNorm() const {
    double sum = 0.0;
    for (double elem : data) {
        sum += elem * elem;
    }
    return std::sqrt(sum);
}

void Vect::initializeFromVec(const Vec& vec) {
    for (int i = 0; i < 3; ++i) {
        setElement(i, vec[i]);  // Assumes Vec defines operator[].
    }
    setElement(3, 1.0);  // Homogeneous coordinate for 3D vectors.
}

Matrix::Matrix() : data(4, std::vector<double>(4, 0.0)) {}  // Default constructor initializing a 4x4 zero matrix.

Matrix::Matrix(std::ifstream& ifs) : data(4, std::vector<double>(4, 0.0)) {
    ifs >> data[0][0];  // a
    ifs >> data[1][1];  // b
    ifs >> data[2][2];  // c
    ifs >> data[0][1]; data[1][0] = data[0][1];  // d
    ifs >> data[1][2]; data[2][1] = data[1][2];  // e
    ifs >> data[0][2]; data[2][0] = data[0][2];  // f
    ifs >> data[0][3]; data[3][0] = data[0][3];  // g
    ifs >> data[1][3]; data[3][1] = data[1][3];  // h
    ifs >> data[2][3]; data[3][2] = data[2][3];  // j
    ifs >> data[3][3];  // k
}

const std::vector<double>& Matrix::getRow(int index) const {
    if (index < 0 || index >= data.size()) throw std::out_of_range("Row index out of range");
    return data[index];
}

std::vector<double>& Matrix::getRow(int index) {
    if (index < 0 ||  index >= data.size()) throw std::out_of_range("Row index out of range");
    return data[index ];
}

Vect Matrix::multiplyByVector(const Vect& vector) const {
    Vect result ; 
    for (int row = 0;  row < 4; row++) {
        double sum = 0.0;

        for (int col =  0; col < 4; col++) {

            sum += data[row][col] * vector.getElement(col);
        }
        result.setElement(row, sum);
    }
    return result;
}



Color operator*(double num, const Color& c)
{
    Color newColor(num * c.redComponent, num * c.greenComponent, num * c.blueComponent);

    return newColor;
}

Ray::Ray() : origin(), endPoint() {}  // Ensuring both vectors are initialized to their defaults

Ray::Ray(const Vec& origin, const Vec& endPoint) : origin(origin), endPoint(endPoint) {}

Vec Ray::computePosition(double t) const {
    Vec directionVector = calculateDirection();
    return origin.add(t * directionVector);  
}

Vec Ray::calculateDirection() const {
    return endPoint.subtract(origin);
}

const Vec& Ray::getOrigin() const {
    return origin;
}

const Vec& Ray::getEndPoint() const {
    return endPoint;
}

// Figure Implementation

Figure::Figure() : ambient(), diffuse(), specular(), reflectivity(), transmissivity(),
shininess(0.0), indexOfRefraction(1.0), reflectionFlag(0), transmissionFlag(0) {}

void Figure::loadMaterialProperties(std::ifstream& ifs) {
    ambient = Color(ifs);
    diffuse = Color(ifs);
    specular = Color(ifs);
    reflectivity = Color(ifs);
    transmissivity = Color(ifs);
    ifs >> shininess >> indexOfRefraction >> reflectionFlag >> transmissionFlag;
}

Color Figure::getAmbientColor() const { return ambient; }
Color Figure::getDiffuseColor() const { return diffuse; }
Color Figure::getSpecularColor() const { return specular; }
Color Figure::getReflectivity() const { return reflectivity; }
Color Figure::getTransmissivity() const { return transmissivity; }
double Figure::getShininess() const { return shininess; }
double Figure::getRefractionIndex() const { return indexOfRefraction; }

bool Figure::hasTransmission() const { return transmissionFlag == 1; }
bool Figure::hasReflection() const { return reflectionFlag == 1; }

// Figure Implementation End

// Light Implementation
//


Light::Light(std::ifstream& ifs) : position(ifs), shading(ifs) {
    readAttenuationCoefficients(ifs);
}

void Light::readAttenuationCoefficients(std::ifstream& ifs) {
    ifs >> constantAttenuation >> linearAttenuation >> quadraticAttenuation;
}

Vec Light::getPosition() const {
    return position;
}

Color Light::getShading() const {
    return shading;
}

double Light::calculateAttenuation(const Vec& intersection) const {
    double distanceSquared = computeDistanceSquared(intersection);
    double distance = std::sqrt(distanceSquared);

    if (linearAttenuation == 0.0 && quadraticAttenuation == 0.0)
        return 1.0 / constantAttenuation;

    return 1.0 / (constantAttenuation + linearAttenuation * distance + quadraticAttenuation * distanceSquared);
}

double Light::computeDistanceSquared(const Vec& point) const {
    Vec separation = point.subtract(position);
    return separation.magnitude();  // If magnitude() returns the distance squared, rename or adjust
}

// Sphere implementation
Sphere::Sphere(std::ifstream& ifs) : Figure(), center(ifs) {
    ifs >> radius;
    loadMaterialProperties(ifs);
}

double Sphere::calculateIntersection(const Ray& ray, double minT, double maximumTimes) const {
    Vec rayDirection = ray.getEndPoint().subtract(ray.getOrigin());
    Vec toRayOrigin = ray.getOrigin().subtract(center);

    double a = rayDirection.dotProduct(rayDirection);
    double b = 2.0 * rayDirection.dotProduct(toRayOrigin);
    double c = toRayOrigin.dotProduct(toRayOrigin) - (radius * radius);

    double discriminant = (b * b) - (4.0 * a * c);
    if (discriminant < 0.0) {
        return maximumTimes + 1.0; // No intersection
    }

    double sqrtDiscriminant = sqrt(discriminant);
    double invDenom = 0.5 / a;  // Inverse of denominator 2a

    // Try the first root
    double t1 = (-b - sqrtDiscriminant) * invDenom;
    if (t1 >= minT && t1 < maximumTimes) {
        return t1;
    }

    // Try the second root
    double t2 = (-b + sqrtDiscriminant) * invDenom;
    if (t2 >= minT && t2 < maximumTimes) {
        return t2;
    }

    return maximumTimes + 1.0; // No valid intersection in the range
}

std::pair<Vec, bool> Sphere::calculateNormal(const Vec& intersectionPoint, const Ray& ray) const {
    Vec normal = intersectionPoint.subtract(center);
    normal.normalize(); // Assuming normalize modifies the vector in place

    bool isFrontFacing = ray.calculateDirection().dotProduct(normal) < 0.0;
    if (!isFrontFacing) {
        normal = -1.0 * normal;
    }

    return { normal, isFrontFacing };
}


// Plane implementation


Plane::Plane(std::ifstream& ifs) : Figure(), normalVector(ifs) {
    ifs >> distanceFromOrigin;
    loadMaterialProperties(ifs);
    direction1 = Vec(ifs);
    direction2 = Vec(ifs);
}

double Plane::calculateIntersection(const Ray& ray, double minT, double maximumTimes) const {
    Vec rayDirection = ray.calculateDirection();
    double numerator = distanceFromOrigin - ray.getOrigin().dotProduct(normalVector);
    double denominator = rayDirection.dotProduct(normalVector);

    if (std::abs(denominator) < epsilon) // Using std::abs for clarity
        return maximumTimes + 1.0;

    double t = numerator / denominator;
    return (t >= minT && t <= maximumTimes) ? t : maximumTimes + 1.0;
}

std::pair<Vec, bool> Plane::calculateNormal(const Vec&, const Ray& ray) const {
    Vec rayDirection = ray.calculateDirection();
    Vec normalizedNormal = normalVector; // Make a copy of the normal vector
    normalizedNormal.normalize(); // Normalize the copy

    bool isFrontFacing = rayDirection.dotProduct(normalizedNormal) < 0;
    return { isFrontFacing ? normalizedNormal : -1.0 * normalizedNormal, isFrontFacing };
}


// Quadric implementation

Quadric::Quadric(std::ifstream& ifs) : Figure(), geometry(ifs) {
    int numClips;
    ifs >> numClips;
    for (int i = 0; i < numClips; i++) {
        clips.push_front(std::make_unique<Matrix>(ifs));
    }
    loadMaterialProperties(ifs);
}

double Quadric::calculateIntersection(const Ray& ray, double minT, double maximumTimes) const {
    Vect p0Vect(ray.getOrigin());
    Vect p1Vect(ray.getEndPoint());
    Vect pDiffVect = p1Vect.subtract(p0Vect);

    double a = pDiffVect.dotProduct(geometry.multiplyByVector(pDiffVect));
    double b = 2.0 * pDiffVect.dotProduct(geometry.multiplyByVector(p0Vect));
    double c = p0Vect.dotProduct(geometry.multiplyByVector(p0Vect));

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < epsilon) return maximumTimes + 1.0; // No real roots

    double sqrtDiscriminant = sqrt(discriminant);
    double invDenom = 0.5 / a;  // Inverse of denominator 2a

    double t1 = (-b - sqrtDiscriminant) * invDenom;
    if (t1 >= minT && t1 < maximumTimes && !clip(ray.computePosition(t1))) return t1;

    double t2 = (-b + sqrtDiscriminant) * invDenom;
    if (t2 >= minT && t2 < maximumTimes && !clip(ray.computePosition(t2))) return t2;

    return maximumTimes + 1.0;
}

std::pair<Vec, bool> Quadric::calculateNormal(const Vec& intersectionPoint, const Ray& ray) const {
    Vect nVect = geometry.multiplyByVector(Vect(intersectionPoint));
    Vec normal(nVect); // Create Vec from Vect
    normal.normalize();

    bool isFrontFacing = ray.calculateDirection().dotProduct(normal) < 0;
    return { isFrontFacing ? normal : -1.0 * normal, isFrontFacing };
}

bool Quadric::clip(const Vec& point) const {
    Vect position(point);
    for (const auto& clipMatrix : clips) {
        if (position.dotProduct((*clipMatrix).multiplyByVector(position)) > 0) return true;
    }
    return false;
}

// END OF HEADER CLASS METHODS








// BEGINNING OF IMPLEMENTATION OF RAY TRACING FUNCTIONS

//



void openFile(ifstream& ifs, const char* filename) {
    ifs.open(filename);
    if (!ifs.is_open()) {
        throw runtime_error("Failed to open file");
    }
}

void parseCameraSettings(ifstream& ifs) {
    ifs >> camera1 >> camera2 >> camera3 >> zCoor;
}

void parseDimensions(ifstream& ifs) {
    ifs >> minimumX >> maximumX >> minimumY >> maximumY;
}

Color readColor(ifstream& ifs) {
    double r, g, b;
    ifs >> r >> g >> b;
    return Color(r, g, b);
}

void calculatePixels() {
    pixW = (maximumX - minimumX) / horizontalRes;
    pixH = (maximumY - minimumY) / verticalRes;
    pixStartX = (minimumX + pixW / 2.0);
    pixStartY = (maximumY - pixH / 2.0);
}

void parseObjects(ifstream& ifs) {
    int numLights, numSpheres, numPlanes, numQuadrics;
    ifs >> numLights;
    for (int i = 0; i < numLights; ++i) lightList.push_front(new Light(ifs));
    ifs >> numSpheres;
    for (int i = 0; i < numSpheres; ++i) shapeList.push_front(new Sphere(ifs));
    ifs >> numPlanes;
    for (int i = 0; i < numPlanes; ++i) shapeList.push_front(new Plane(ifs));
    ifs >> numQuadrics;
    for (int i = 0; i < numQuadrics; ++i) shapeList.push_front(new Quadric(ifs));
}


void sceneParsing(char* sceneName) {
    assert(sceneName != nullptr);
    ifstream ifs;
    openFile(ifs, sceneName);
    parseCameraSettings(ifs);
    parseDimensions(ifs);
    backgroundColor = readColor(ifs);
    ambient = readColor(ifs);
    ifs >> maxDepth;
    ifs >> horizontalRes >> verticalRes;
    calculatePixels();
    parseObjects(ifs);
    ifs.close();
}



void createImageStorage() {
    image.reserve(verticalRes);
    for (int i = 0; i < verticalRes; ++i) {
        image.push_back(vector<Color>());
        image[i].reserve(horizontalRes);
    }
}

void initializePixelValues() {
    for (int i = 0; i < verticalRes; ++i) {
        for (int j = 0; j < horizontalRes; ++j) {
            image[i].push_back(Default_black);
        }
    }
}

void beginImageGeneration() {
    createImageStorage();
    initializePixelValues();
}

void writeP3Header(ofstream& ofs, const char* outputFile) {
    ofs << "P3" << "\n"
        << horizontalRes << " " << verticalRes << "\n"
        << 255 << "\n";
}



void openOutputFile(char* outputFile, ofstream& ofs) {
    assert(outputFile != nullptr);
    ofs.open(outputFile);
    if (!ofs.is_open()) {
        throw runtime_error("Failed to open output file.");
    }
    writeP3Header(ofs, outputFile);
}



// Calculates the reflection direction given the ray and surface normal
Vec calculateReflectionDirection(const Ray& ray, const Vec& normal)
{
    Vec incomingDirection = ray.calculateDirection();
    incomingDirection.normalize();
    double incidenceDotNormal = incomingDirection.dotProduct(normal);
    return (2 * incidenceDotNormal * normal).subtract(incomingDirection);
}

// Calculates the color from a reflected ray
Color calculateReflectedColor(Figure* obj, const Ray& reflectedRay, double depth)
{
    Color reflectionFactor = obj->getReflectivity();
    Color recursiveColor = RT_trace(reflectedRay, depth + 1);
    return reflectionFactor.multiply(recursiveColor);
}

// Reflects a ray based on the interaction with an object's surface
Color RT_reflect(Figure* obj, const Ray& ray, const Vec& intersection, const Vec& normal, double depth)
{
    if (depth > maxDepth || !obj->hasReflection()) {
        return Default_black;
    }

    Vec reflectionDirection = calculateReflectionDirection(ray, normal);
    Ray reflectedRay(intersection, intersection.add(reflectionDirection));
    return calculateReflectedColor(obj, reflectedRay, depth);
}







// Computes the index ratio based on whether the ray is entering or exiting the material
double calculateIndexRatio(Figure* obj, bool entering)
{
    double refractionIndex = obj->getRefractionIndex();
    return entering ? 1.0 / refractionIndex : refractionIndex;
}

// Calculates the direction of a transmitted ray considering refraction
Vec calculateTransmissionDirection(const Ray& ray, const Vec& normal, double indexRatio)
{
    Vec incomingDirection = ray.calculateDirection();  // Get the incoming direction
    incomingDirection.normalize();  // Normalize the copied vector

    double incidenceDotNormal = incomingDirection.dotProduct(normal);
    double k = 1.0 - indexRatio * indexRatio * (1.0 - incidenceDotNormal * incidenceDotNormal);

    if (k < 0.0) {
        return Vec();  // Indicates total internal reflection
    }

    return ((indexRatio * incidenceDotNormal - sqrt(k)) * normal).subtract(indexRatio * incomingDirection);
}

// Calculates the color from a transmitted ray
Color calculateTransmittedColor(Figure* obj, const Ray& transmittedRay, double depth)
{
    Color transmissivity = obj->getTransmissivity();
    Color recursiveColor = RT_trace(transmittedRay, depth + 1);
    return transmissivity.multiply(recursiveColor);
}

// Transmits a ray through an object considering refraction
Color RT_transmit(Figure* obj, const Ray& ray, const Vec& intersection, const Vec& normal, bool entering, double depth)
{
    if (depth > maxDepth || !obj->hasTransmission()) {
        return Default_black;
    }

    double indexRatio = calculateIndexRatio(obj, entering);
    Vec transmissionDirection = calculateTransmissionDirection(ray, normal, indexRatio);

    if (transmissionDirection.magnitude() == 0) {  // Handles total internal reflection
        return Default_black;
    }

    Ray transmittedRay(intersection, intersection.add(transmissionDirection));
    return calculateTransmittedColor(obj, transmittedRay, depth);
}



// Calculates the diffuse component of shading based on the light's interaction with an object
Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
    return (dotProduct > 0.0) ? dotProduct * light->getShading().multiply(obj->getDiffuseColor()) : Default_black;
}

// Helper to calculate reflected direction for specular lighting
Vec calculateReflectDirection(const Vec& normal, const Vec& lightDirection, double dotProduct)
{
    return (2.0 * dotProduct * normal).subtract(lightDirection);
}

// Computes intensity of specular reflection
double calculateSpecularIntensity(const Ray& ray, const Vec& reflectDirection, double shininess)
{
    Vec viewDirection = -1 * ray.calculateDirection();
    viewDirection.normalize();
    double dotProduct = viewDirection.dotProduct(reflectDirection);
    return (dotProduct > 0.0) ? pow(dotProduct, shininess) : 0.0;
}

// Calculates the specular component of shading
Color specularShade(Figure* obj, const Vec& normal, Light* light, const Vec& lightDirection, double dotProduct, const Ray& ray)
{
    Vec reflectDirection = calculateReflectDirection(normal, lightDirection, dotProduct);
    double specularity = calculateSpecularIntensity(ray, reflectDirection, obj->getShininess());
    return (specularity > 0.0) ? specularity * light->getShading().multiply(obj->getSpecularColor()) : Default_black;
}




// Finds the nearest intersection of a ray with any figure in the scene, considering transparency
pair<double, Figure*> nearestIntersection(const Ray& r, double minT, double maximumTimes, bool mayBeTransparent = true)
{
    double nearestT = maximumTimes + epsilon;
    Figure* nearestFigure = NULL;

    for (auto& figure : shapeList) {
        if (mayBeTransparent || !figure->hasTransmission()) {
            double currentT = figure->calculateIntersection(r, minT, maximumTimes);
            if (currentT < nearestT && currentT >= minT) {
                nearestT = currentT;
                nearestFigure = figure;
            }
        }
    }
    return { nearestT, nearestFigure };
}


// Check if the light source is visible from the intersection point
bool isVisible(const Vec& intersection, Light* light)
{
    Ray lightRay(intersection, light->getPosition());
    pair<double, Figure*> intersectionInfo = nearestIntersection(lightRay, epsilon, 1.0, false);
    return (intersectionInfo.second == NULL);
}

// Calculates the color contribution from a single light source
Color calculateLightContribution(Figure* obj, const Ray& ray, const Vec& intersection, const Vec& normal, Light* light)
{
    // Make a copy of the direction vector and normalize the copy
    Vec lightDirection = light->getPosition().subtract(intersection);
    lightDirection.normalize();  // Normalize the local copy

    double dotProduct = normal.dotProduct(lightDirection);
    Color color = diffuseShade(obj, light, dotProduct);
    color = color.add(specularShade(obj, normal, light, lightDirection, dotProduct, ray));

    return color;
}




// Calculate the combined color from all lights interacting with the intersection point
Color RT_lights(Figure* obj, const Ray& ray, const Vec& intersection, const Vec& normal)
{
    Color totalColor;
    for (Light* light : lightList)
    {
        if (isVisible(intersection, light)) {
            double attenuation = light->calculateAttenuation(intersection);
            totalColor = totalColor.add(attenuation * calculateLightContribution(obj, ray, intersection, normal, light));
        }
    }
    return totalColor;
}



// Computes the ambient component of the shading
Color calculateAmbientShading(Figure* obj)
{
    return ambient.multiply(obj->getAmbientColor());
}


// Calculates the final shading at a point taking into account ambient, diffuse, and specular contributions as well as reflections and transmissions
Color RT_shade(Figure* obj, const Ray& ray, const Vec& intersection, const Vec& normal, bool entering, double depth)
{
    Color shadedColor = calculateAmbientShading(obj);
    shadedColor = shadedColor.add(RT_lights(obj, ray, intersection, normal));

    if (depth < maxDepth) {
        shadedColor = shadedColor.add(RT_reflect(obj, ray, intersection, normal, depth));
        shadedColor = shadedColor.add(RT_transmit(obj, ray, intersection, normal, entering, depth));
    }

    return shadedColor;
}



// Trace a ray into the scene and determine the color at the intersection point
Color RT_trace(const Ray& r, double depth)
{
    pair<double, Figure*> intersection = nearestIntersection(r, epsilon, maximumTimes);
    if (!intersection.second) {
        return backgroundColor;  // No intersection, return the background color
    }

    // Calculate intersection point and normal
    Vec intersectionPoint = r.computePosition(intersection.first);
    pair<Vec, bool> normalData = intersection.second->calculateNormal(intersectionPoint, r);
    Vec normal = normalData.first;
    bool entering = normalData.second;

    // Compute shading for the intersection
    return RT_shade(intersection.second, r, intersectionPoint, normal, entering, depth);
}


// Calculate the center coordinates of a pixel in the image grid
pair<double, double> pixelCenter(int i, int j)
{
    double x = pixStartX + j * pixW;
    double y = pixStartY - i * pixH;
    return { x, y };  // Using modern C++ initializer list
}



// Main rendering loop that processes each pixel to generate the final image
void RayTraceAlgo()
{
    Vec cameraOrigin(camera1, camera2, camera3);  // Define the camera origin
    for (int i = 0; i < verticalRes; i++)
    {
        for (int j = 0; j < horizontalRes; j++)
        {
            pair<double, double> pc = pixelCenter(i, j);  // Get pixel center
            double x = pc.first;
            double y = pc.second;
            Vec pixelPosition(x, y, zCoor);   // Position in space for the pixel
            Ray ray(cameraOrigin, pixelPosition);  // Create a ray from camera to pixel
            image[i][j] = RT_trace(ray, 1);  // Trace the ray and store the result in the image
        }
    }
}



double getMaxColorDifference(const Color& pixelColor, const Color& neighborColor) {
    return std::max(0.0, pixelColor.difference(neighborColor));
}

double maxDifference(int i, int j) {
    double maximum = 0.0;
    const Color& pixelColor = image[i][j];

    for (int m = -1; m <= 1; m++) {
        int mi = i + m;
        if (mi < 0 || mi >= verticalRes) continue;

        for (int n = -1; n <= 1; n++) {
            int nj = j + n;
            if (nj < 0 || nj >= horizontalRes) continue;
            if (m != 0 || n != 0) {
                const Color& neighborColor = image[mi][nj];
                maximum = std::max(maximum, getMaxColorDifference(pixelColor, neighborColor));
            }
        }
    }
    return maximum;
}

Color traceSubPixel(const Vec& cameraOrigin, double x, double y, int depth) {
    Vec pixelPosition(x, y, zCoor);
    Ray ray(cameraOrigin, pixelPosition);
    return RT_trace(ray, depth);
}

void superSample(int i, int j) {
    Vec cameraOrigin(camera1, camera2, camera3);
    pair<double, double> pc = pixelCenter(i, j);
    double pcX = pc.first;
    double pcY = pc.second;

    int numSamples = 2 * ssRate + 1;
    double subPixelWidth = pixW / numSamples;
    double subPixelHeight = pixH / numSamples;
    Color accumulate;

    for (int ii = -ssRate; ii <= ssRate; ii++) {
        for (int jj = -ssRate; jj <= ssRate; jj++) {
            double x = pcX + jj * subPixelWidth;
            double y = pcY + ii * subPixelHeight;
            if (ii == 0 && jj == 0)
                accumulate = accumulate.add(image[i][j]); // Use the original pixel color for the center
            else
                accumulate = accumulate.add(traceSubPixel(cameraOrigin, x, y, 1)); // Trace and add the color from super sampled pixel
        }
    }

    Color average = (1.0 / (numSamples * numSamples))*accumulate; // Calculate average color
    image[i][j] = average;
}


vector<vector<bool>> determineSupersamplingNeeds() {
    vector<vector<bool>> needsSuperSample(verticalRes, vector<bool>(horizontalRes, false));
    for (int i = 0; i < verticalRes; i++) {
        for (int j = 0; j < horizontalRes; j++) {
            needsSuperSample[i][j] = (maxDifference(i, j) > aaThreshold);
        }
    }
    return needsSuperSample;
}

void performSuperSampling(const vector<vector<bool>>& needsSuperSample) {
    for (int i = 0; i < verticalRes; i++) {
        for (int j = 0; j < horizontalRes; j++) {
            if (needsSuperSample[i][j]) {
                superSample(i, j);
            }
        }
    }
}

void totalAntiAliasing() {
    vector<vector<bool>> needsSuperSample = determineSupersamplingNeeds();
    performSuperSampling(needsSuperSample);
}


void outputImage()
{
    ofstream ofs;
    openOutputFile("picture.ppm", ofs);
    for (int i = 0; i < verticalRes; i++)
        for (int j = 0; j < horizontalRes; j++)
            image[i][j].writeToStream(ofs);
    ofs.close();
}

int main(int, char* argv[])
{
    sceneParsing(argv[1]);

    beginImageGeneration();

    RayTraceAlgo();

    totalAntiAliasing();
    outputImage();
    return 0;
}

