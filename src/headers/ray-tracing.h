#include <list>
using namespace std;

class Color;
class Vec;
class Vect;
class Matrix;
class Ray;
class Light;
class Figure;
class Plane;
class Sphere;
class Quadric;


// Color Header

class Color {
    friend Color operator*(double scalar, const Color& color);

protected:
    double redComponent, greenComponent, blueComponent;

public:
    Color();
    Color(std::ifstream& input);
    Color(double red, double green, double blue);
    void writeToStream(std::ostream& output);
    Color add(const Color& other) const;
    Color subtract(const Color& other) const;
    Color multiply(const Color& other) const;
    double magnitude() const;
    double difference(const Color& other) const;
};

// Vec Header


class Vec {
    friend Vec operator*(double scalar, const Vec& vector);

private:
    double components[3];

public:
    Vec();
    Vec(std::ifstream& input);
    Vec(double x, double y, double z);
    Vec(const Vec& other);

    void normalize();
    double magnitude() const;
    double dotProduct(const Vec& other) const;
    Vec add(const Vec& other) const;
    Vec subtract(const Vec& other) const;

    const double& operator[](int index) const;
    double& operator[](int index);

    explicit Vec(const Vect& v);
};


// Vect header


class Vect {
private:
    std::vector<double> data;

public:
    Vect();
    explicit Vect(const Vec& vec);

    // Basic Operations
    Vect subtract(const Vect& other) const;

    // Element Access
    double getElement(int index) const;
    void setElement(int index, double value);

    // Vector Operations
    double dotProduct(const Vect& other) const;
    double calculateNorm() const;

    // Utility
    void initializeFromVec(const Vec& vec);
};

// Matrix Header

class Matrix {
private:
    std::vector<std::vector<double>> data;

public:
    Matrix();
    explicit Matrix(std::ifstream& ifs);

    const std::vector<double>& getRow(int index) const;
    std::vector<double>& getRow(int index);

    Vect multiplyByVector(const Vect& vector) const;  // More descriptive function name
};


// Ray Header
class Ray {
private:
    Vec origin;
    Vec endPoint;

public:
    Ray();
    Ray(const Vec& origin, const Vec& endPoint);

    Vec computePosition(double t) const;  // Explicitly named for clarity
    Vec calculateDirection() const;       // Same as above
    const Vec& getOrigin() const;         // Renamed for clarity
    const Vec& getEndPoint() const;       // Renamed for clarity
};

/// <summary>
/// Light Header
/// </summary>


class Light {
public:
    explicit Light(std::ifstream& ifs);

    Vec getPosition() const;
    Color getShading() const;
    double calculateAttenuation(const Vec& intersection) const;

private:
    Vec position;
    Color shading;
    double constantAttenuation;
    double linearAttenuation;
    double quadraticAttenuation;

    void readAttenuationCoefficients(std::ifstream& ifs);
    double computeDistanceSquared(const Vec& point) const;
};


class Figure {
protected:
    Color ambient;
    Color diffuse;
    Color specular;
    Color reflectivity;
    Color transmissivity;
    double shininess;
    double indexOfRefraction;
    int reflectionFlag, transmissionFlag;

public:
    Figure();

    // Initialization
    void loadMaterialProperties(std::ifstream& ifs);

    // Property Accessors
    Color getAmbientColor() const;
    Color getDiffuseColor() const;
    Color getSpecularColor() const;
    Color getReflectivity() const;
    Color getTransmissivity() const;
    double getShininess() const;
    double getRefractionIndex() const;
    bool hasTransmission() const;
    bool hasReflection() const;

    // Pure virtual functions for derived classes
    virtual double calculateIntersection(const Ray& r, double minT, double maxT) const = 0;
    virtual std::pair<Vec, bool> calculateNormal(const Vec& v, const Ray& r) const = 0;
};


// Plane Header
class Plane : public Figure {
private:
    Vec direction1;
    Vec direction2;
    Vec normalVector; 
    double distanceFromOrigin; 
    

public:
    explicit Plane(std::ifstream& ifs);

    virtual double calculateIntersection(const Ray& ray, double minT, double maxT) const override;
    virtual std::pair<Vec, bool> calculateNormal(const Vec& intersectionPoint, const Ray& ray) const override;
};

// Sphere header

class Sphere : public Figure {
private:
    Vec center;
    double radius;

public:
    explicit Sphere(std::ifstream& ifs);

    virtual double calculateIntersection(const Ray& ray, double minT, double maxT) const override;
    virtual std::pair<Vec, bool> calculateNormal(const Vec& intersectionPoint, const Ray& ray) const override;
};

// Quadric Header

class Quadric : public Figure {
private:
    Matrix geometry;
    std::list<std::unique_ptr<Matrix>> clips;  // Using smart pointers for better memory management

public:
    explicit Quadric(std::ifstream& ifs);

    virtual double calculateIntersection(const Ray& ray, double minT, double maxT) const override;
    virtual std::pair<Vec, bool> calculateNormal(const Vec& intersectionPoint, const Ray& ray) const override;
    bool clip(const Vec& point) const;
};


Vec operator*(double num, const Vec& v);

Color operator*(double num, const Color& c);

Color RT_trace(const Ray& r, double depth);
