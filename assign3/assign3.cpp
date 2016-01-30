/*
 CSCI 480
 Assignment 3 rayDirectiontracer
 Name: Ryan Yuan Lu
 */

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/Dense"
#include <math.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

float kEpsilon = 1e-8;
float currentT;
float currentU;
float currentV;

using Eigen::MatrixXd;
using Eigen::Vector3d;
using namespace std;

char* filename = 0;

//different display modes
#define MODE_DISPLAY 10
#define MODE_JPEG 2
int mode = MODE_DISPLAY;

//you may want to make  these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
// int HEIGHT = 480;
// int WIDTH = 640;

double aspectRatio;

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex {
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

typedef struct _Triangle {
    struct Vertex v[3];
} Triangle;

typedef struct _Sphere {
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
} Sphere;

typedef struct _Light {
    double position[3];
    double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int currentIntersectedObj = 0; // 1: triangle, 2: Sphere
Triangle currentIntersectedTri;
Sphere currentIntersectedSphere;

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

double currentIntersectionPos[3];

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

#define EPSILON 0.000001

int triangle_intersection(const Vector3d V1, // Triangle vertices
    const Vector3d V2,
    const Vector3d V3,
    const Vector3d O, //Ray origin
    const Vector3d D, //Ray direction
    float* out)
{
    Vector3d e1, e2; //Edge1, Edge2
    Vector3d P, Q, T;
    float det, inv_det, u, v;
    float t;

    //Find vectors for two edges sharing V1
    e1 = V2 - V1;
    e2 = V3 - V1;
    // SUB(e1, V2, V1);
    // SUB(e2, V3, V1);
    //Begin calculating determinant - also used to calculate u parameter
    P = D.cross(e2);
    // CROSS(P, D, e2);
    //if determinant is near zero, ray lies in plane of triangle
    det = e1.dot(P);
    // det = DOT(e1, P);
    //NOT CULLING
    if (det > -EPSILON && det < EPSILON)
        return 0;
    inv_det = 1.f / det;

    //calculate distance from V1 to ray origin
    T = O - V1;
    // SUB(T, O, V1);

    //Calculate u parameter and test bound
    u = T.dot(P) * inv_det;
    // u = DOT(T, P) * inv_det;
    //The intersection lies outside of the triangle
    if (u < 0.f || u > 1.f)
        return 0;

    //Prepare to test v parameter
    Q = T.cross(e1);
    // CROSS(Q, T, e1);

    //Calculate V parameter and test bound
    v = D.dot(Q) * inv_det;
    // v = DOT(D, Q) * inv_det;
    //The intersection lies outside of the triangle
    if (v < 0.f || u + v > 1.f)
        return 0;
    t = e2.dot(Q) * inv_det;
    // t = DOT(e2, Q) * inv_det;

    if (t > EPSILON) { //ray intersection
        *out = t;
        return 1;
    }

    // No hit, no win
    return 0;
}

bool rayTriangleIntersect(
    const Vector3d& orig, const Vector3d& dir,
    const Vector3d& v0, const Vector3d& v1, const Vector3d& v2,
    float& t, float& u, float& v)
{
#ifdef MOLLER_TRUMBORE
    Vector3d v0v1 = v1 - v0;
    Vector3d v0v2 = v2 - v0;
    Vector3d pvec = dir.cross(v0v2);
    float det = v0v1.dot(pvec);
#ifdef CULLING
    // if the determinant is negative the triangle is backfacing
    // if the determinant is close to 0, the rayDirection misses the triangle
    if (det < kEpsilon)
        return false;
#else
    // rayDirection and triangle are parallel if det is close to 0
    if (fabs(det) < kEpsilon)
        return false;
#endif
    float invDet = 1 / det;

    Vector3d tvec = orig - v0;
    u = tvec.dot(pvec) * invDet;
    if (u < 0 || u > 1)
        return false;

    Vector3d qvec = tvec.cross(v0v1);
    v = dir.dot(qvec) * invDet;
    if (v < 0 || u + v > 1)
        return false;

    t = v0v2.dot(qvec) * invDet;

    return true;
#else
    // compute plane's normal
    Vector3d v0v1 = v1 - v0;
    Vector3d v0v2 = v2 - v0;
    // no need to normalize
    Vector3d N = v0v1.cross(v0v2); // N
    float denom = N.dot(N);

    // Step 1: finding P

    // check if rayDirection and plane are parallel ?
    float NdotrayDirection = N.dot(dir);
    if (fabs(NdotrayDirection) < kEpsilon) // almost 0
        return false; // they are parallel so they don't intersect !

    // compute d parameter using equation 2
    float d = N.dot(v0);

    // compute t (equation 3)
    t = (N.dot(orig) + d) / NdotrayDirection;
    // check if the triangle is in behind the rayDirection
    if (t < 0)
        return false; // the triangle is behind

    // compute the intersection point using equation 1
    Vector3d P = orig + t * dir;

    // Step 2: inside-outside test
    Vector3d C; // vector perpendicular to triangle's plane

    // edge 0
    Vector3d edge0 = v1 - v0;
    Vector3d vp0 = P - v0;
    C = edge0.cross(vp0);
    if (N.dot(C) < 0)
        return false; // P is on the right side

    // edge 1
    Vector3d edge1 = v2 - v1;
    Vector3d vp1 = P - v1;
    C = edge1.cross(vp1);
    if ((u = N.dot(C)) < 0)
        return false; // P is on the right side

    // edge 2
    Vector3d edge2 = v0 - v2;
    Vector3d vp2 = P - v2;
    C = edge2.cross(vp2);
    if ((v = N.dot(C)) < 0)
        return false; // P is on the right side;

    u /= denom;
    v /= denom;

    return true; // this rayDirection hits the triangle
#endif
}

bool checkIntersections(Vector3d rayDirection)
{
    bool intTri = false;
    bool intSphere = false;
    Triangle t;
    Sphere s;
    float triTVals[num_triangles];
    float sphereTVals[num_spheres];
    float triUvals[num_triangles];
    float triVvals[num_triangles];

    bool doesIntersectTriangle[num_spheres];
    for (int i = 0; i < num_triangles; ++i) {
        triTVals[i] = 5000;
        triVvals[i] = 100;
        triUvals[i] = 100;
        doesIntersectTriangle[i] = false;
    }
    bool doesHaveTriangleInt = false;

    for (int i = 0; i < num_triangles; ++i) {
        t = triangles[i];
        Vector3d v0 = Vector3d(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
        Vector3d v1 = Vector3d(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
        Vector3d v2 = Vector3d(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);

        float u;
        float v;
        doesIntersectTriangle[i] = rayTriangleIntersect(Vector3d(0, 0, 0), rayDirection, v0, v1, v2, triTVals[i], triUvals[i], triVvals[i]);
        if (triVvals[i] < kEpsilon || !doesIntersectTriangle[i]) {
            doesIntersectTriangle[i] = false;
            triTVals[i] = 5000;
            triVvals[i] = 100;
            triUvals[i] = 100;
        }
    }

    bool doesIntersectSphere[num_spheres];
    for (int i = 0; i < num_spheres; ++i) {
        sphereTVals[i] = 5000;
        doesIntersectSphere[i] = false;
    }
    bool doesHaveSphereInt = false;

    for (int i = 0; i < num_spheres; ++i) {
        Vector3d unitVec = Vector3d(0, 0, 0);
        s = spheres[i];
        double a = rayDirection.x() * rayDirection.x() + rayDirection.y() * rayDirection.y() + rayDirection.z() * rayDirection.z();
        double b = 2 * (rayDirection.x() * (unitVec.x() - s.position[0]) + rayDirection.y() * (unitVec.y() - s.position[1]) + rayDirection.z() * (unitVec.z() - s.position[2]));
        double c = (unitVec.x() - s.position[0]) * (unitVec.x() - s.position[0]) + (unitVec.y() - s.position[1]) * (unitVec.y() - s.position[1]) + (unitVec.z() - s.position[2]) * (unitVec.z() - s.position[2]) - s.radius * s.radius;
        if ((b * b - 4 * c) < 0) {

            doesIntersectSphere[i] = false;
            sphereTVals[i] = 5000;
        }
        else {

            double t0 = (-b + sqrt(b * b - 4 * c)) / 2;
            double t1 = (-b - sqrt(b * b - 4 * c)) / 2;
            if (t0 > kEpsilon && t1 > kEpsilon) {
                // printf(" here 3\n");
                if (t1 < t0) {
                    sphereTVals[i] = t1;
                }
                else if (t0 < t1) {
                    sphereTVals[i] = t0;
                    /* code */
                }

                doesIntersectSphere[i] = true;
            }
            else if (t0 > kEpsilon && t1 <= kEpsilon) {
                sphereTVals[i] = t0;
                doesIntersectSphere[i] = true;
                /* code */
            }
            else if (t1 > kEpsilon && t0 <= kEpsilon) {
                sphereTVals[i] = t1;
                doesIntersectSphere[i] = true;
                /* code */
            }
            else {
                doesIntersectSphere[i] = false;
            }
        }
    }
    for (int i = 0; i < num_triangles; ++i) {

        if (doesIntersectTriangle[i]) {
            doesHaveTriangleInt = true;
        }
    }

    for (int i = 0; i < num_spheres; ++i) {

        if (doesIntersectSphere[i]) {
            doesHaveSphereInt = true;
        }
    }
    double minTriTVal = 5000;
    int triCurIndex = -1;
    double minSphTval = 5000;
    int sphCurIndex = -1;

    if (doesHaveTriangleInt && doesHaveSphereInt) {

        for (int i = 0; i < num_triangles; ++i) {

            if (triTVals[i] < minTriTVal) {
                minTriTVal = triTVals[i];
                triCurIndex = i;
            }
        }

        for (int i = 0; i < num_spheres; ++i) {
            if (sphereTVals[i] < minSphTval) {
                minSphTval = sphereTVals[i];
                sphCurIndex = i;
            }
        }

        if (minTriTVal < minSphTval) {
            currentT = minTriTVal;
            currentIntersectedObj = 1;
            currentU = triUvals[triCurIndex];
            currentV = triVvals[triCurIndex];
            currentIntersectedTri = triangles[triCurIndex];
        }
        else if (minSphTval < minTriTVal) {
            currentT = minSphTval;

            currentIntersectedObj = 2;
            currentIntersectedSphere = spheres[sphCurIndex];
        }
        return true;
    }
    else if (doesHaveTriangleInt) {
        for (int i = 0; i < num_triangles; ++i) {

            if (triTVals[i] < minTriTVal && triTVals[i] > kEpsilon) {
                minTriTVal = triTVals[i];
                triCurIndex = i;
            }
        }

        currentT = minTriTVal;
        currentU = triUvals[triCurIndex];
        currentV = triVvals[triCurIndex];
        currentIntersectedObj = 1;
        currentIntersectedTri = triangles[triCurIndex];
        return true;
    }
    else if (doesHaveSphereInt) {
        for (int i = 0; i < num_spheres; ++i) {
            if (sphereTVals[i] < minSphTval) {
                minSphTval = sphereTVals[i];
                sphCurIndex = i;
            }
        }
        currentT = minSphTval;
        currentIntersectedObj = 2;
        currentIntersectedSphere = spheres[sphCurIndex];

        return true;
    }
    currentT = -1;
    currentIntersectedObj = 0;
    return false;
}

bool checkIfShadowrayDirectionHit(Vector3d rayDir, Vector3d intPoint)
{
    // return false;
    bool intTri = false;
    bool intSphere = false;
    Triangle t;
    Sphere s;
    Vector3d rayDirection = rayDir.normalized();

    bool doesIntersectTriangle[num_triangles];
    int doesIntersectTriangle2[num_triangles];
    float tValCheckers[num_triangles];

    for (int i = 0; i < num_spheres; ++i) {
        tValCheckers[i] = -1.0;
        doesIntersectTriangle[i] = false;
        doesIntersectTriangle2[i] = 0;
    }
    bool doesHaveTriangleInt = false;
    for (int i = 0; i < num_triangles; ++i) {
        t = triangles[i];

        Vector3d v0 = Vector3d(t.v[0].position[0], t.v[0].position[1], t.v[0].position[2]);
        Vector3d v1 = Vector3d(t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]);
        Vector3d v2 = Vector3d(t.v[2].position[0], t.v[2].position[1], t.v[2].position[2]);

        float v = 0;
        float u = 0;
        float t = 0;
        doesIntersectTriangle2[i] = triangle_intersection(v0, v1, v2, intPoint, rayDirection, &tValCheckers[i]);
        if (doesIntersectTriangle2[i] == 0) {
            tValCheckers[i] = -1;
            doesIntersectTriangle[i] = false;
        }
        else {
            doesIntersectTriangle[i] = true;
        }
    }

    bool doesIntersectSphere[num_spheres];
    for (int i = 0; i < num_spheres; ++i) {
        doesIntersectSphere[i] = false;
    }
    bool doesHaveSphereInt = false;
    for (int i = 0; i < num_spheres; ++i) {
        Vector3d unitVec = intPoint;
        s = spheres[i];

        // calculate a, b ,c for the sphere
        double a = rayDirection.x() * rayDirection.x() + rayDirection.y() * rayDirection.y() + rayDirection.z() * rayDirection.z();
        double b = 2 * (rayDirection.x() * (unitVec.x() - s.position[0]) + rayDirection.y() * (unitVec.y() - s.position[1]) + rayDirection.z() * (unitVec.z() - s.position[2]));
        double c = (unitVec.x() - s.position[0]) * (unitVec.x() - s.position[0]) + (unitVec.y() - s.position[1]) * (unitVec.y() - s.position[1]) + (unitVec.z() - s.position[2]) * (unitVec.z() - s.position[2]) - s.radius * s.radius;
        double m = b * b - 4 * c;

        if (m <= 0) {
            doesIntersectSphere[i] = false;
        }
        else {

            double t0 = -b + sqrt(b * b - 4 * c) / 2;
            double t1 = -b - sqrt(b * b - 4 * c) / 2;

            if (t0 > kEpsilon || t1 > kEpsilon) {

                doesIntersectSphere[i] = true;
            }
        }
    }

    for (int i = 0; i < num_triangles; ++i) {
        if (doesIntersectTriangle2[i] == 1) {

            doesHaveTriangleInt = true;
        }
    }

    for (int i = 0; i < num_spheres; ++i) {
        if (doesIntersectSphere[i]) {

            doesHaveSphereInt = true;
        }
    }
    if (doesHaveTriangleInt == 1 || doesHaveSphereInt) {

        return true;
    }
    return false;
}

void getLightColor(Vector3d intPoint, double& lightColorX, double& lightColorY, double& lightColorZ)
{
    double I[3];
    double lightCoef[3];
    double diffusedCoef[3];
    double specularCoef[3];
    Vector3d dirToLight;
    Vector3d normal;
    Vector3d dirToCam;
    Vector3d refDir;
    double shininess;

    dirToCam = Vector3d((0 - intPoint.x()), (0 - intPoint.y()), (0 - intPoint.z())).normalized();
    //initialize light color by ambient
    I[0] = ambient_light[0];
    I[1] = ambient_light[1];
    I[2] = ambient_light[2];

    if (currentIntersectedObj == 2) {
        Sphere obj = currentIntersectedSphere;
        normal = 1 / obj.radius * Vector3d((intPoint.x() - obj.position[0]), (intPoint.y() - obj.position[1]), (intPoint.z() - obj.position[2]));
        shininess = obj.shininess;

        for (int i = 0; i < num_lights; ++i) {
            Light l = lights[i];
            Vector3d lightPosition = Vector3d(l.position[0], l.position[1], l.position[2]);
            Vector3d dirToLight;
            dirToLight = (lightPosition - intPoint).normalized();
            refDir = (2 * (dirToLight.dot(normal)) * normal - dirToLight).normalized();

            if (checkIfShadowrayDirectionHit(dirToLight, intPoint)) {
                // if blocked
                I[0] += 0;
                I[1] += 0;
                I[2] += 0;
            }
            else {
                // if no blockage
                //computer light color 3 times for 3 vertices, then interpolate
                for (int k = 0; k < 3; ++k) {

                    lightCoef[k] = l.color[k];
                    diffusedCoef[k] = obj.color_diffuse[k];
                    specularCoef[k] = obj.color_specular[k];

                    if (dirToLight.dot(normal) < 0 && refDir.dot(dirToCam) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (0) + specularCoef[k] * (pow(0, shininess)));
                    }
                    else if (dirToLight.dot(normal) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (0) + specularCoef[k] * (pow(refDir.dot(dirToCam), shininess)));
                    }
                    else if (refDir.dot(dirToCam) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (dirToLight.dot(normal)) + specularCoef[k] * (pow(0, shininess)));
                    }
                    else {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (dirToLight.dot(normal)) + specularCoef[k] * (pow(refDir.dot(dirToCam), shininess)));
                    }
                    if (I[k] > 1.0) {
                        I[k] = 1.0;
                    }
                }
            }

            lightColorX = I[0];
            lightColorY = I[1];
            lightColorZ = I[2];
        }
    }
    else if (currentIntersectedObj == 1) {

        Triangle tri = currentIntersectedTri;

        Vector3d n1 = Vector3d(tri.v[0].normal[0], tri.v[0].normal[1], tri.v[0].normal[2]);
        Vector3d n2 = Vector3d(tri.v[1].normal[0], tri.v[1].normal[1], tri.v[1].normal[2]);
        Vector3d n3 = Vector3d(tri.v[2].normal[0], tri.v[2].normal[1], tri.v[2].normal[2]);

        normal = (1 - currentU - currentV) * n1 + currentU * n2 + currentV * n3;

        for (int i = 0; i < num_lights; ++i) {
            Light l = lights[i];
            Vector3d dirToLight;
            Vector3d lightPosition = Vector3d(l.position[0], l.position[1], l.position[2]);
            dirToLight = (lightPosition - intPoint).normalized();
            refDir = (2 * (dirToLight.dot(normal)) * normal - dirToLight).normalized();

            if (checkIfShadowrayDirectionHit(dirToLight, intPoint)) {
                //if blocked, stack black color
                I[0] += 0;
                I[1] += 0;
                I[2] += 0;
            }
            else {
                //if doesnt have blockage
                double shininess1 = 0;
                double shininess2 = 0;
                double shininess3 = 0;

                double diffusedCoef1[3];
                double diffusedCoef2[3];
                double diffusedCoef3[3];

                double specularCoef1[3];
                double specularCoef2[3];
                double specularCoef3[3];

                shininess1 = tri.v[0].shininess;
                shininess2 = tri.v[1].shininess;
                shininess3 = tri.v[2].shininess;

                //interpolate shininess
                shininess = (1 - currentU - currentV) * shininess1 + currentU * shininess2 + currentV * shininess3;

                for (int k = 0; k < 3; ++k) {

                    lightCoef[k] = l.color[k];

                    //diffuse colors for 3 vertices
                    diffusedCoef1[k] = tri.v[0].color_diffuse[k];
                    diffusedCoef2[k] = tri.v[1].color_diffuse[k];
                    diffusedCoef3[k] = tri.v[2].color_diffuse[k];

                    //specular colors for 3 vertices
                    specularCoef1[k] = tri.v[0].color_specular[k];
                    specularCoef2[k] = tri.v[1].color_specular[k];
                    specularCoef3[k] = tri.v[2].color_specular[k];

                    //interpolate the colors
                    diffusedCoef[k] = (1 - currentU - currentV) * diffusedCoef1[k] + currentU * diffusedCoef2[k] + currentV * diffusedCoef3[k];
                    specularCoef[k] = (1 - currentU - currentV) * specularCoef1[k] + currentU * specularCoef2[k] + currentV * specularCoef3[k];

                    //check if any of the dot products are negative, then computer light color
                    if (dirToLight.dot(normal) < 0 && refDir.dot(dirToCam) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (0) + specularCoef[k] * (pow(0, shininess)));
                    }
                    else if (dirToLight.dot(normal) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (0) + specularCoef[k] * (pow(refDir.dot(dirToCam), shininess)));
                    }
                    else if (refDir.dot(dirToCam) < 0) {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (dirToLight.dot(normal)) + specularCoef[k] * (pow(0, shininess)));
                    }
                    else {
                        I[k] += lightCoef[k] * (diffusedCoef[k] * (dirToLight.dot(normal)) + specularCoef[k] * (pow(refDir.dot(dirToCam), shininess)));
                    }
                    if (I[k] > 1.0) {
                        I[k] = 1.0;
                    }
                }
            }
            lightColorX = I[0];
            lightColorY = I[1];
            lightColorZ = I[2];
        }
    }
}

void draw_scene()
{
    unsigned int x, y;
    aspectRatio = (double)WIDTH / (double)HEIGHT;

    Vector3d rayDirection;
    Vector3d rayDir;
    Vector3d rayOrigin;
    Vector3d intPoint;
    rayOrigin = Vector3d(0, 0, 0);
    rayDir = Vector3d(0, 0, -1);
    double fovR = fov * 3.14159 / 180;

    //find edges of image plane
    double topLeftX = -aspectRatio * tan(fovR / 2);
    double topLeftY = tan(fovR / 2);
    double topLeftZ = -1;

    double topRightX = aspectRatio * tan(fovR / 2);
    double topRightY = tan(fovR / 2);
    double topRightZ = -1;

    double bottomLeftX = -(double)aspectRatio * tan(fovR / 2);
    double bottomLeftY = -tan(fovR / 2);
    double bottomLeftZ = -1;

    for (x = 0; x < WIDTH; x++) {
        //1. send out rayDirections
        //find x for ray
        rayDir.x() = topLeftX + topRightX * 2 / WIDTH * x;
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for (y = 0; y < HEIGHT; y++) {
            //find y for ray
            rayDir.y() = -(topLeftY + bottomLeftY * 2 / HEIGHT * y);
            //2. check intersection
            if (checkIntersections(rayDir.normalized())) {
                intPoint = rayOrigin + rayDir.normalized() * currentT;
                double lightColorX = 0;
                double lightColorY = 0;
                double lightColorZ = 0;
                //3. illumination & shadows
                getLightColor(intPoint, lightColorX, lightColorY, lightColorZ);
                //4. plot based on color found
                plot_pixel(x, y, lightColorX * 255, lightColorY * 255, lightColorZ * 255);
            }
            else {
                //if no intersection we return white background
                plot_pixel(x, y, 255, 255, 255);
            }
        }
        glEnd();
        glFlush();
    }
    printf("Done!\n");
    fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    glColor3f(((double)r) / 256.f, ((double)g) / 256.f, ((double)b) / 256.f);
    glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[HEIGHT - y - 1][x][0] = r;
    buffer[HEIGHT - y - 1][x][1] = g;
    buffer[HEIGHT - y - 1][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x, y, r, g, b);
    if (mode == MODE_JPEG)
        plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
    Pic* in = NULL;

    in = pic_alloc(640, 480, 3, NULL);
    printf("Saving JPEG file: %s\n", filename);

    memcpy(in->pix, buffer, 3 * WIDTH * HEIGHT);
    if (jpeg_write(filename, in))
        printf("File saved Successfully\n");
    else
        printf("Error in Saving\n");

    pic_free(in);
}

void parse_check(char* expected, char* found)
{
    if (strcasecmp(expected, found)) {
        char error[100];
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, char* check, double p[3])
{
    char str[100];
    fscanf(file, "%s", str);
    parse_check(check, str);
    fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
    printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE* file, double* r)
{
    char str[100];
    fscanf(file, "%s", str);
    parse_check("rad:", str);
    fscanf(file, "%lf", r);
    printf("rad: %f\n", *r);
}

void parse_shi(FILE* file, double* shi)
{
    char s[100];
    fscanf(file, "%s", s);
    parse_check("shi:", s);
    fscanf(file, "%lf", shi);
    printf("shi: %f\n", *shi);
}

int loadScene(char* argv)
{
    FILE* file = fopen(argv, "r");
    int number_of_objects;
    char type[50];
    int i;
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file, "%i", &number_of_objects);

    printf("number of objects: %i\n", number_of_objects);
    char str[200];

    parse_doubles(file, "amb:", ambient_light);

    for (i = 0; i < number_of_objects; i++) {
        fscanf(file, "%s\n", type);
        printf("%s\n", type);
        if (strcasecmp(type, "triangle") == 0) {

            printf("found triangle\n");
            int j;

            for (j = 0; j < 3; j++) {
                parse_doubles(file, "pos:", t.v[j].position);
                parse_doubles(file, "nor:", t.v[j].normal);
                parse_doubles(file, "dif:", t.v[j].color_diffuse);
                parse_doubles(file, "spe:", t.v[j].color_specular);
                parse_shi(file, &t.v[j].shininess);
            }

            if (num_triangles == MAX_TRIANGLES) {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if (strcasecmp(type, "sphere") == 0) {
            printf("found sphere\n");

            parse_doubles(file, "pos:", s.position);
            parse_rad(file, &s.radius);
            parse_doubles(file, "dif:", s.color_diffuse);
            parse_doubles(file, "spe:", s.color_specular);
            parse_shi(file, &s.shininess);

            if (num_spheres == MAX_SPHERES) {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if (strcasecmp(type, "light") == 0) {
            printf("found light\n");
            parse_doubles(file, "pos:", l.position);
            parse_doubles(file, "col:", l.color);

            if (num_lights == MAX_LIGHTS) {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else {
            printf("unknown type in scene description:\n%s\n", type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    // printf("HERE\n");
    static int once = 0;
    if (!once) {
        draw_scene();
        if (mode == MODE_JPEG)
            save_jpg();
    }
    once = 1;
}

void MyKeyboardFunc(unsigned char Key, int x, int y)
{
    switch (Key) {
    case '1': {
        printf("Pressed Z\n");
        num_lights = 0;
        num_spheres = 0;
        num_triangles = 0;
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
        loadScene("test2.scene");
        // init();

        draw_scene();

        break;
    }
    case '2': {
        printf("Pressed Z\n");
        num_lights = 0;
        num_spheres = 0;
        num_triangles = 0;
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
        loadScene("spheres.scene");
        // init();

        draw_scene();

        break;
    }
    case '3': {
        printf("Pressed Z\n");
        num_lights = 0;
        num_spheres = 0;
        num_triangles = 0;
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
        loadScene("table.scene");
        // init();

        draw_scene();

        break;
    }
    case '4': {
        printf("Pressed Z\n");
        num_lights = 0;
        num_spheres = 0;
        num_triangles = 0;
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
        loadScene("test2-o.scene");
        // init();

        draw_scene();

        break;
    }
    case '5': {
        printf("Pressed Z\n");
        num_lights = 0;
        num_spheres = 0;
        num_triangles = 0;
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);
        loadScene("test2-ng.scene");
        // init();

        draw_scene();

        break;
    }
    default:
        break;
    };
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4) {
        printf("usage: %s <scenefile> [jpegname] (type'D' to turn on debug mode)\n", argv[0]);
        exit(0);
    }
    printf("ARC:: %d\n", argc);
    if (argc == 3) {
        if (strcmp(argv[2], "D") == 0) {
            printf("not goes here\n");
            mode = MODE_DISPLAY;
            //turn on debug mode
        }
        else {
            printf("goes here\n");
            mode = MODE_JPEG;
            filename = argv[2];
        }
    }
    else if (argc == 2)
        mode = MODE_DISPLAY;
    else if (argc == 4) {
        //turn on debug mode
    }

    glutInit(&argc, argv);
    loadScene(argv[1]);

    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);

    glutInitWindowSize(WIDTH, HEIGHT);
    int window = glutCreateWindow("rayDirection Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutKeyboardFunc(MyKeyboardFunc);
    glutMainLoop();
}
