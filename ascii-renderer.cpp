#include <float.h>
#include <math.h>

#include <iostream>

using namespace std;

const double EXTERNAL_RADIUS = 2;
const double INTERNAL_RADIUS = 1.2;

const int CANVAS_WIDTH = 60;
const int CANVAS_HEIGHT = 60;

const double FOV_WIDTH = 8;
const double FOV_HEIGHT = 8;

const double CAMERA_Z = 10;
const double FOV_Z = 3;

const double LIGHT_X = -10;
const double LIGHT_Y = -10;
const double LIGHT_Z = 10;

const double U_INTERVAL_START = 0;
const double U_INTERVAL_END = 2 * M_PI;
const double U_INTERVAL_STEP = 0.05;

const double V_INTERVAL_START = 0;
const double V_INTERVAL_END = 2 * M_PI;
const double V_INTERVAL_STEP = 0.05;

/*
  Vector3D represents a 3D vector of doubles:
  
  [x, y, z]
*/

struct Vector3D {
  double x, y, z;
};

/*
  Vector2D represents a 2D vector of integers:
  
  [x, y]
*/

struct Vector2D {
  int x, y;
};

/*
  Matrix3X3 represents a 3X3 matrix of doubles:
  
  | x1 y1 z1 |
  | x2 y2 z2 |
  | x3 y3 z3 |
*/

struct Matrix3X3 {
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
};

/* Vectorial operations */

double magnitude(Vector3D vector) {
  return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

Vector3D computeVector(Vector3D pointA, Vector3D pointB) {
  struct Vector3D vector;

  vector.x = pointB.x - pointA.x;
  vector.y = pointB.y - pointA.y;
  vector.z = pointB.z - pointA.z;

  return vector;
}

Vector3D computeVersor(Vector3D vector) {
  struct Vector3D versor;

  double length = magnitude(vector);

  versor.x = vector.x / length;
  versor.y = vector.y / length;
  versor.z = vector.z / length;

  return versor;
}

double scalarProduct(Vector3D vectorA, Vector3D vectorB) {
  return vectorA.x * vectorB.x + vectorA.y * vectorB.y + vectorA.z * vectorB.z;
}

double angleBetweenVectors(Vector3D vectorA, Vector3D vectorB) {
  double vectorAMagnitude = magnitude(vectorA);
  double vectorBMagnitude = magnitude(vectorB);

  if (vectorAMagnitude == 0 || vectorBMagnitude == 0) {
    return M_PI;
  }

  return acos(scalarProduct(vectorA, vectorB) / (vectorAMagnitude * vectorBMagnitude));
}

Vector3D projectOnPlane(Vector3D sourcePoint, Vector3D vector, double planeZ) {
  if (vector.z == 0) {
    struct Vector3D projection;

    projection.x = DBL_MAX;
    projection.y = DBL_MAX;
    projection.z = planeZ;

    return projection;
  }

  double distance = abs(sourcePoint.z - planeZ);

  struct Vector3D projection;

  projection.x = (vector.x / abs(vector.z)) * distance;
  projection.y = (vector.y / abs(vector.z)) * distance;
  projection.z = planeZ;

  return projection;
}

Vector3D reflectAcrossNormal(Vector3D vector, Vector3D normal) {
  struct Vector3D reflection;

  if (angleBetweenVectors(vector, normal) > M_PI / 2) {
    reflection.x = 0;
    reflection.y = 0;
    reflection.z = 0;

    return reflection;
  }

  double vectorAndNormalScalarProduct = scalarProduct(vector, normal);
  double squaredNormalMagnitude = pow(magnitude(normal), 2);
  double coefficient = 2 * vectorAndNormalScalarProduct / squaredNormalMagnitude;

  reflection.x = coefficient * normal.x - vector.x;
  reflection.y = coefficient * normal.y - vector.y;
  reflection.z = coefficient * normal.z - vector.z;

  return reflection;
}

Vector3D matrixVector3DProduct(Matrix3X3 matrix, Vector3D vector) {
  struct Vector3D result;

  result.x = matrix.x1 * vector.x + matrix.y1 * vector.y + matrix.z1 * vector.z;
  result.y = matrix.x2 * vector.x + matrix.y2 * vector.y + matrix.z2 * vector.z;
  result.z = matrix.x3 * vector.x + matrix.y3 * vector.y + matrix.z3 * vector.z;

  return result;
}

Matrix3X3 computeRotationMatrix(double xAngle, double yAngle, double zAngle) {
  struct Matrix3X3 result;

  double cosXAngle = cos(xAngle);
  double sinXAngle = sin(xAngle);
  double cosYAngle = cos(yAngle);
  double sinYAngle = sin(yAngle);
  double cosZAngle = cos(zAngle);
  double sinZAngle = sin(zAngle);

  result.x1 = cosYAngle * cosZAngle;
  result.x2 = -cosYAngle * sinZAngle;
  result.x3 = sinYAngle;
  result.y1 = cosZAngle * sinXAngle * sinYAngle + cosXAngle * sinZAngle;
  result.y2 = cosZAngle * cosXAngle - sinXAngle * sinYAngle * sinZAngle;
  result.y3 = -cosYAngle * sinXAngle;
  result.z1 = -cosXAngle * cosZAngle * sinYAngle + sinXAngle * sinZAngle;
  result.z2 = cosZAngle * sinXAngle + cosXAngle * sinYAngle * sinZAngle;
  result.z3 = cosXAngle * cosYAngle;

  return result;
}

/* Toroidal surface  */

double toroidalSurfaceX(double R, double r, double u, double v) {
  return (R + r * cos(u)) * cos(v);
}

double toroidalSurfaceY(double R, double r, double u, double v) {
  return (R + r * cos(u)) * sin(v);
}

double toroidalSurfaceZ(double R, double r, double u, double v) {
  return r * sin(u);
}

Vector3D computeToroidalSurfacePoint(double R, double r, double u, double v) {
  struct Vector3D point;

  point.x = toroidalSurfaceX(R, r, u, v);
  point.y = toroidalSurfaceY(R, r, u, v);
  point.z = toroidalSurfaceZ(R, r, u, v);

  return point;
}

/* Toroidal normal vectors */

double toroidalSurfaceNormalX(double R, double r, double u, double v) {
  return toroidalSurfaceX(R, r, u, v) * r * cos(u);
}

double toroidalSurfaceNormalY(double R, double r, double u, double v) {
  return toroidalSurfaceY(R, r, u, v) * r * cos(u);
}

double toroidalSurfaceNormalZ(double R, double r, double u, double v) {
  return toroidalSurfaceZ(R, r, u, v) * (R + r * cos(u));
}

Vector3D computeToroidalSurfaceNormalVersor(double R, double r, double u, double v) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = toroidalSurfaceNormalX(R, r, u, v);
  vector.y = toroidalSurfaceNormalY(R, r, u, v);
  vector.z = toroidalSurfaceNormalZ(R, r, u, v);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

/* Möbius Strip surface  */

double mobiusSurfaceX(double u, double v) {
  return (1 + (v / 2) * cos(u / 2)) * cos(u);
}

double mobiusSurfaceY(double u, double v) {
  return (1 + (v / 2) * cos(u / 2)) * sin(u);
}

double mobiusSurfaceZ(double u, double v) {
  return (v / 2) * sin(u / 2);
}

Vector3D computeMobiusSurfacePoint(double u, double v) {
  struct Vector3D point;

  point.x = mobiusSurfaceX(u, v);
  point.y = mobiusSurfaceY(u, v);
  point.z = mobiusSurfaceZ(u, v);

  return point;
}

/* Möbius Strip normal vectors */

double mobiusSurfaceNormalX(double R, double r, double u, double v) {
  return (1 / 2) * cos(u) * sin(u / 2) - (v / 4) * cos(u / 2) * sin(u / 2) * sin(u) - (v / 8) * sin(u);
}

double mobiusSurfaceNormalY(double R, double r, double u, double v) {
  return (1 / 2) * sin(u) * sin(u / 2) + (v / 4) * cos(u / 2) * sin(u / 2) * sin(u) + (v / 8) * cos(u);
}

double mobiusSurfaceNormalZ(double R, double r, double u, double v) {
  return -(1 / 2) * cos(u / 2) * (1 + (v / 2) * cos(u / 2));
}

Vector3D computeMobiusSurfaceNormalVersor(double R, double r, double u, double v) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = mobiusSurfaceNormalX(R, r, u, v);
  vector.y = mobiusSurfaceNormalY(R, r, u, v);
  vector.z = mobiusSurfaceNormalZ(R, r, u, v);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

/* Helpers */

Vector2D projectionToCanvasMatrix(
    Vector3D projection,
    double fovWidth,
    double fovHeight,
    int canvasWidth,
    int canvasHeight) {
  struct Vector2D matrixCoordinates;

  double horizontalUnit = fovWidth / canvasWidth;
  double verticalUnit = fovHeight / canvasHeight;

  double horizontalCoordinate = projection.x / horizontalUnit;
  double verticalCoordinate = projection.y / verticalUnit;

  matrixCoordinates.x = round(horizontalCoordinate + (canvasWidth / 2));
  matrixCoordinates.y = round(verticalCoordinate + (canvasHeight / 2));

  return matrixCoordinates;
}

char computeIntensityCharacter(double intensity, double intensityMaxValue = 1) {
  int charactersCount = 12;
  char characters[] = {'.', ',', '-', '~', ':', ';', '=', '!', '*', '#', '$', '@'};

  int characterIndex = round(intensity / intensityMaxValue * (charactersCount - 1));

  return characters[characterIndex];
}

void resetMatrices(
    double depthBufferMatrix[][CANVAS_WIDTH],
    double absoluteIntensityMatrix[][CANVAS_WIDTH]) {
  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      depthBufferMatrix[y][x] = -DBL_MAX;
      absoluteIntensityMatrix[y][x] = 0;
    }
  }
}

void render(
    Vector3D cameraPoint,
    Vector3D lightPoint,
    double depthBufferMatrix[][CANVAS_WIDTH],
    double absoluteIntensityMatrix[][CANVAS_WIDTH],
    double xAngle = 0,
    double yAngle = 0,
    double zAngle = 0) {
  struct Matrix3X3 rotationMatrix = computeRotationMatrix(xAngle, yAngle, zAngle);

  for (double u = U_INTERVAL_START; u < U_INTERVAL_END; u = u + U_INTERVAL_STEP) {
    for (double v = V_INTERVAL_START; v < V_INTERVAL_END; v = v + V_INTERVAL_STEP) {
      struct Vector3D surfacePoint = computeToroidalSurfacePoint(EXTERNAL_RADIUS, INTERNAL_RADIUS, u, v);
      struct Vector3D surfaceNormal = computeToroidalSurfaceNormalVersor(EXTERNAL_RADIUS, INTERNAL_RADIUS, u, v);

      struct Vector3D rotatedSurfacePoint = matrixVector3DProduct(rotationMatrix, surfacePoint);
      struct Vector3D rotatedSurfaceNormal = matrixVector3DProduct(rotationMatrix, surfaceNormal);

      struct Vector3D surfaceToCameraVector = computeVector(rotatedSurfacePoint, cameraPoint);
      struct Vector3D surfaceToCameraVersor = computeVersor(surfaceToCameraVector);
      struct Vector3D lightVector = computeVector(rotatedSurfacePoint, lightPoint);
      struct Vector3D reflectionVector = reflectAcrossNormal(lightVector, rotatedSurfaceNormal);
      struct Vector3D projection = projectOnPlane(cameraPoint, surfaceToCameraVersor, FOV_Z);

      double reflectionAngle = angleBetweenVectors(surfaceToCameraVector, reflectionVector);
      double lightIntensity = 1 - (reflectionAngle / M_PI);

      struct Vector2D matrixCoordinates = projectionToCanvasMatrix(
          projection,
          FOV_WIDTH,
          FOV_HEIGHT,
          CANVAS_HEIGHT,
          CANVAS_WIDTH);

      if (rotatedSurfacePoint.z > depthBufferMatrix[matrixCoordinates.y][matrixCoordinates.x]) {
        absoluteIntensityMatrix[matrixCoordinates.y][matrixCoordinates.x] = lightIntensity;
        depthBufferMatrix[matrixCoordinates.y][matrixCoordinates.x] = rotatedSurfacePoint.z;
      }
    }
  }
}

void print(
    double depthBufferMatrix[][CANVAS_WIDTH],
    double absoluteIntensityMatrix[][CANVAS_WIDTH]) {
  double maxLightIntensity = 0;

  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      if (absoluteIntensityMatrix[y][x] > maxLightIntensity) {
        maxLightIntensity = absoluteIntensityMatrix[y][x];
      }
    }
  }

  cout << "\x1b[H";

  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      if (depthBufferMatrix[y][x] == -DBL_MAX) {
        cout << ' ';
      } else {
        cout << computeIntensityCharacter(absoluteIntensityMatrix[y][x], maxLightIntensity);
      }
    }

    cout << endl;
  }
}

/* Main */

int main() {
  double depthBufferMatrix[CANVAS_HEIGHT][CANVAS_WIDTH];
  double absoluteIntensityMatrix[CANVAS_HEIGHT][CANVAS_WIDTH];

  for (double i = 0;; i += 1) {
    /* Light vector */

    struct Vector3D lightPoint;

    lightPoint.x = LIGHT_X;
    lightPoint.y = LIGHT_Y;
    lightPoint.z = LIGHT_Z;

    /* Camera point */

    struct Vector3D cameraPoint;

    cameraPoint.x = 0;
    cameraPoint.y = 0;
    cameraPoint.z = CAMERA_Z;

    /* Rotation angles */

    double xAngle = M_PI * i / 300;
    double yAngle = 0;
    double zAngle = M_PI * i / 240;

    /* Resetting matrices */

    resetMatrices(depthBufferMatrix, absoluteIntensityMatrix);

    /* Rendering */

    render(
        cameraPoint,
        lightPoint,
        depthBufferMatrix,
        absoluteIntensityMatrix,
        xAngle,
        yAngle,
        zAngle);

    /* Printing */

    print(depthBufferMatrix, absoluteIntensityMatrix);
  }

  /* Exiting */

  return 0;
}
