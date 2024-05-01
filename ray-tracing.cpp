#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <iostream>

using namespace std;

#include "ray-tracing.h"


vector< vector<Color> > image;

list<Figure*> shapeList;
list<Light*> lightList;


const double epsilon = 0.00000000001;
const int colorScale = 255;

const double aaThreshold = 0.25;
int ssRate = 1;

double cameraX, cameraY, cameraZ;
int horizontalResolution, verticalResolution;
double zCoor;
double minX, maxX;
double minY, maxY;
double pixelWidth;
double pixelHeight;
double pixelStartX;
double pixelStartY;

int maxDepth;
Color backgroundColor;
Color ambient;

const Color BLACK_COLOR;
const double maxT = 1000000;




// Utility functions for clamping and getting maximum values



// Color Implementation

int clamp(int value, int minVal, int maxVal) { return (value < minVal) ? minVal : (value > maxVal ? maxVal : value); }
int getMax(int first, int second) { return (first > second) ? first : second; }

Color::Color() : redComponent(0.0), greenComponent(0.0), blueComponent(0.0) {}

Color::Color(std::ifstream& input) {
    input >> redComponent >> greenComponent >> blueComponent;
}

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

double Sphere::calculateIntersection(const Ray& ray, double minT, double maxT) const {
    Vec rayDirection = ray.getEndPoint().subtract(ray.getOrigin());
    Vec toRayOrigin = ray.getOrigin().subtract(center);

    double a = rayDirection.dotProduct(rayDirection);
    double b = 2.0 * rayDirection.dotProduct(toRayOrigin);
    double c = toRayOrigin.dotProduct(toRayOrigin) - (radius * radius);

    double discriminant = (b * b) - (4.0 * a * c);
    if (discriminant < 0.0) {
        return maxT + 1.0; // No intersection
    }

    double sqrtDiscriminant = sqrt(discriminant);
    double invDenom = 0.5 / a;  // Inverse of denominator 2a

    // Try the first root
    double t1 = (-b - sqrtDiscriminant) * invDenom;
    if (t1 >= minT && t1 < maxT) {
        return t1;
    }

    // Try the second root
    double t2 = (-b + sqrtDiscriminant) * invDenom;
    if (t2 >= minT && t2 < maxT) {
        return t2;
    }

    return maxT + 1.0; // No valid intersection in the range
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

double Plane::calculateIntersection(const Ray& ray, double minT, double maxT) const {
    Vec rayDirection = ray.calculateDirection();
    double numerator = distanceFromOrigin - ray.getOrigin().dotProduct(normalVector);
    double denominator = rayDirection.dotProduct(normalVector);

    if (std::abs(denominator) < epsilon) // Using std::abs for clarity
        return maxT + 1.0;

    double t = numerator / denominator;
    return (t >= minT && t <= maxT) ? t : maxT + 1.0;
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

double Quadric::calculateIntersection(const Ray& ray, double minT, double maxT) const {
    Vect p0Vect(ray.getOrigin());
    Vect p1Vect(ray.getEndPoint());
    Vect pDiffVect = p1Vect.subtract(p0Vect);

    double a = pDiffVect.dotProduct(geometry.multiplyByVector(pDiffVect));
    double b = 2.0 * pDiffVect.dotProduct(geometry.multiplyByVector(p0Vect));
    double c = p0Vect.dotProduct(geometry.multiplyByVector(p0Vect));

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < epsilon) return maxT + 1.0; // No real roots

    double sqrtDiscriminant = sqrt(discriminant);
    double invDenom = 0.5 / a;  // Inverse of denominator 2a

    double t1 = (-b - sqrtDiscriminant) * invDenom;
    if (t1 >= minT && t1 < maxT && !clip(ray.computePosition(t1))) return t1;

    double t2 = (-b + sqrtDiscriminant) * invDenom;
    if (t2 >= minT && t2 < maxT && !clip(ray.computePosition(t2))) return t2;

    return maxT + 1.0;
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
    ifs >> cameraX >> cameraY >> cameraZ >> zCoor;
}

void parseDimensions(ifstream& ifs) {
    ifs >> minX >> maxX >> minY >> maxY;
}

Color readColor(ifstream& ifs) {
    double r, g, b;
    ifs >> r >> g >> b;
    return Color(r, g, b);
}

void calculatePixels() {
    pixelWidth = (maxX - minX) / horizontalResolution;
    pixelHeight = (maxY - minY) / verticalResolution;
    pixelStartX = (minX + pixelWidth / 2.0);
    pixelStartY = (maxY - pixelHeight / 2.0);
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


void parseSceneFile(char* sceneName) {
    assert(sceneName != nullptr);
    ifstream ifs;
    openFile(ifs, sceneName);
    parseCameraSettings(ifs);
    parseDimensions(ifs);
    backgroundColor = readColor(ifs);
    ambient = readColor(ifs);
    ifs >> maxDepth;
    ifs >> horizontalResolution >> verticalResolution;
    calculatePixels();
    parseObjects(ifs);
    ifs.close();
}



void createImageStorage() {
    image.reserve(verticalResolution);
    for (int i = 0; i < verticalResolution; ++i) {
        image.push_back(vector<Color>());
        image[i].reserve(horizontalResolution);
    }
}

void initializePixelValues() {
    for (int i = 0; i < verticalResolution; ++i) {
        for (int j = 0; j < horizontalResolution; ++j) {
            image[i].push_back(BLACK_COLOR);
        }
    }
}

void initializeImage() {
    createImageStorage();
    initializePixelValues();
}

void writeP3Header(ofstream& ofs, const char* outputFile) {
    ofs << "P3" << "\n"
        << horizontalResolution << " " << verticalResolution << "\n"
        << colorScale << "\n";
}



void openOutputFile(char* outputFile, ofstream& ofs) {
    assert(outputFile != nullptr);
    ofs.open(outputFile);
    if (!ofs.is_open()) {
        throw runtime_error("Failed to open output file.");
    }
    writeP3Header(ofs, outputFile);
}



Color RT_reflect(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal, double depth)
{
    if (depth <= maxDepth && obj->hasReflection())
    {
        Vec viewDirection = -1.0 * ray.calculateDirection();
        viewDirection.normalize();
        double dotProduct = viewDirection.dotProduct(normal);
        Vec reflectDirection = ((2.0 * dotProduct) * normal).subtract(viewDirection);
        Ray reflectedRay(i, i.add(reflectDirection));
        return obj->getReflectivity().multiply(RT_trace(reflectedRay, depth + 1));
    }
    else return BLACK_COLOR;
}

Color RT_transmit(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
    bool entering, double depth)
{
    if (depth <= maxDepth && obj->hasTransmission())
    {
        Vec viewDirection = -1.0 * ray.calculateDirection();
        viewDirection.normalize();
        double dotProduct = viewDirection.dotProduct(normal);
        double dotProduct2 = dotProduct * dotProduct;
        double indexRatio = entering ? 1.0 / obj->getRefractionIndex() : obj->getRefractionIndex();
        double indexRatio2 = indexRatio * indexRatio;
        double temp1 = indexRatio * dotProduct;
        double temp2 = 1.0 - indexRatio2 + indexRatio2 * dotProduct2;
        if (temp2 >= 0.0)
        {
            double temp3 = temp1 - sqrt(temp2);
            Vec transmitDirection = (temp3 * normal).subtract(indexRatio * viewDirection);
            Ray transmittedRay(i, i.add(transmitDirection));
            return obj->getTransmissivity().multiply(RT_trace(transmittedRay, depth + 1));
        }
        else return BLACK_COLOR;
    }
    else return BLACK_COLOR;

}

Color diffuseShade(Figure* obj, Light* light, double dotProduct)
{
    if (dotProduct > 0.0)
        return dotProduct * (light->getShading().multiply(obj->getDiffuseColor()));
    else return BLACK_COLOR;
}

Color specularShade(Figure* obj, const Vec& normal,
    Light* light, const Vec& lightDirection, double dotProduct,
    const Ray& ray)
{
    Vec reflectDirection = ((2.0 * dotProduct) * normal).subtract(lightDirection);
    Vec viewDirection = -1.0 * ray.calculateDirection();
    viewDirection.normalize();
    double dotProduct2 = viewDirection.dotProduct(reflectDirection);
    if (dotProduct2 > 0.0)
    {
        double shineFactor = pow(dotProduct2, obj->getShininess());
        return shineFactor * (light->getShading().multiply(obj->getSpecularColor()));
    }
    return BLACK_COLOR;
}


pair<double, Figure*> nearestIntersection(const Ray& r,
    double minT, double maxT,
    bool mayBeTransparent = true)
{
    double t = maxT + epsilon;
    Figure* f = NULL;
    for (list<Figure*>::iterator iter = shapeList.begin();
        iter != shapeList.end();
        ++iter)
    {
        if (mayBeTransparent || !((*iter)->hasTransmission()))
        {
            double newT = (*iter)->calculateIntersection(r, minT, maxT);
            if (newT < t && newT >= minT)
            {
                t = newT;
                f = *iter;
            }
        }
    }
    return pair<double, Figure*>(t, f);
}


Color RT_lights(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal)
{
    Color color;
    for (list<Light*>::iterator iter = lightList.begin();
        iter != lightList.end();
        ++iter)
    {
        Light* light = *iter;
        Ray L_Ray(i, light->getPosition());
        double attenuation = light->calculateAttenuation(i);
        pair<double, Figure*> nint = nearestIntersection(L_Ray, epsilon, 1.0, false);
        if (nint.second == NULL)
        {
            Vec lightDirection = L_Ray.calculateDirection();
            lightDirection.normalize();
            double dotProduct = lightDirection.dotProduct(normal);
            color = color.add(attenuation * (diffuseShade(obj, light, dotProduct)));
            color = color.add(attenuation * (specularShade(obj, normal, light,
                lightDirection, dotProduct,
                ray)));
        }
    }
    return color;
}


Color RT_shade(Figure* obj, const Ray& ray, const Vec& i, const Vec& normal,
    bool entering, double depth)
{
    Color newColor = ambient.multiply(obj->getAmbientColor());
    newColor = newColor.add(RT_lights(obj, ray, i, normal));
    if (depth < maxDepth)
    {
        newColor = newColor.add(RT_reflect(obj, ray, i, normal, depth));
        newColor = newColor.add(RT_transmit(obj, ray, i, normal, entering, depth));
    }
    return newColor;
}



Color RT_trace(const Ray& r, double depth)
{

    pair<double, Figure*> intersection = nearestIntersection(r, epsilon, maxT);
    if (intersection.second == NULL) return backgroundColor;
    else
    {
        Figure* figure = intersection.second;
        Vec i = r.computePosition(intersection.first);
        pair<Vec, bool> normalData = figure->calculateNormal(i, r);
        Vec normal = normalData.first;
        bool entering = normalData.second;
        return RT_shade(figure, r, i, normal, entering, depth);
    }
}

pair<double, double> pixelCenter(int i, int j)
{
    double x = pixelStartX + j * pixelWidth;
    double y = pixelStartY - i * pixelHeight;
    return pair<double, double>(x, y);
}

void RT_algorithm()
{
    Vec p0(cameraX, cameraY, cameraZ);
    for (int i = 0; i < verticalResolution; i++)
    {
        for (int j = 0; j < horizontalResolution; j++)
        {
            pair<double, double> pc = pixelCenter(i, j);
            double x = pc.first;
            double y = pc.second;
            Vec p1(x, y, zCoor);
            Ray R(p0, p1);
            image[i][j] = RT_trace(R, 1);
        }
      
    }
}

double maxDifference(int i, int j)
{
    double maximum = 0.0;
    Color& pixelColor = image[i][j];
    for (int m = -1; m < 2; m++)
    {
        if (i + m == -1 || i + m == verticalResolution) continue;
        for (int n = -1; n < 2; n++)
        {
            if (j + n == -1 || j + n == horizontalResolution) continue;
            if (m != 0 || n != 0)
            {
                Color& neighbor = image[i + m][j + n];
                double difference = pixelColor.difference(neighbor);
                maximum = getMax(maximum, difference);
            }
        }
    }
    return maximum;
}

void superSample(int i, int j)
{
    Color accumulate;
    Vec p0(cameraX, cameraY, cameraZ);
    pair<double, double> pc = pixelCenter(i, j);
    double pcX = pc.first;
    double pcY = pc.second;
    int num = 2 * ssRate + 1;
    double w = pixelWidth / num;
    double h = pixelHeight / num;
    for (int ii = -ssRate; ii <= ssRate; ii++)
        for (int jj = -ssRate; jj <= ssRate; jj++)
        {
            if (ii == 0 && jj == 0) accumulate = accumulate.add(image[i][j]);
            else {
                double x = pcX + jj * w;
                double y = pcY + ii * h;
                Vec p1(x, y, zCoor);
                Ray R(p0, p1);
                accumulate = accumulate.add(RT_trace(R, 1));
            }
        }
    Color average = (1.0 / (num * num)) * accumulate;
    image[i][j] = average;
}

void antiAliasing()
{
    vector< vector<bool> > needsSuperSample;
    for (int i = 0; i < verticalResolution; i++)
    {
        needsSuperSample.push_back(vector<bool>());
        needsSuperSample[i].reserve(horizontalResolution);
        for (int j = 0; j < horizontalResolution; j++)
            if (maxDifference(i, j) > aaThreshold)
                needsSuperSample[i].push_back(true);
            else needsSuperSample[i].push_back(false);
    }
    for (int i = 0; i < verticalResolution; i++)
    {
        for (int j = 0; j < horizontalResolution; j++)
            if (needsSuperSample[i][j])  superSample(i, j);
        
    }
}

void writeImageFile()
{
    ofstream ofs;
    openOutputFile("picture.ppm", ofs);
    for (int i = 0; i < verticalResolution; i++)
        for (int j = 0; j < horizontalResolution; j++)
            image[i][j].writeToStream(ofs);
    ofs.close();
}

int main(int, char* argv[])
{
    parseSceneFile(argv[1]);
    initializeImage();
    RT_algorithm();
    antiAliasing();
    writeImageFile();
    return 0;
}
