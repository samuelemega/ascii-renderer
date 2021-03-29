/**
 * 
 * Author: Samuele Mega <samuelemega@gmail.com>
 * Date: 2021-03-28
 *
 * Inspired by Andy Sloane's donut (https://www.a1k0n.net/2006/09/15/obfuscated-c-donut.html).
 *
 * Intended for didactic purposes.
 * 
 * */

#include <float.h>
#include <math.h>

#include <iostream>

using namespace std;

const int CANVAS_WIDTH = 60;
const int CANVAS_HEIGHT = 60;

const double LIGHT_X = 10;
const double LIGHT_Y = -10;
const double LIGHT_Z = 10;

const double CAMERA_Z = 10;

const double TOROID_EXTERNAL_RADIUS = 2;
const double TOROID_INTERNAL_RADIUS = 1.2;

const double TOROID_FOV_WIDTH = 6;
const double TOROID_FOV_HEIGHT = 6;
const double TOROID_FOV_Z = 3;

const double TOROID_U_INTERVAL_START = 0;
const double TOROID_U_INTERVAL_END = 2 * M_PI;
const double TOROID_U_INTERVAL_STEP = 0.05;

const double TOROID_V_INTERVAL_START = 0;
const double TOROID_V_INTERVAL_END = 2 * M_PI;
const double TOROID_V_INTERVAL_STEP = 0.02;

const int TOROID_U_STEPS_COUNT = (TOROID_U_INTERVAL_END - TOROID_U_INTERVAL_START) / TOROID_U_INTERVAL_STEP;
const int TOROID_V_STEPS_COUNT = (TOROID_V_INTERVAL_END - TOROID_V_INTERVAL_START) / TOROID_V_INTERVAL_STEP;
const int TOROID_PRECOMPUTED_DATA_COUNT = 4;

const double TOROID_X_ANGLE_PERIOD = 1000;
const double TOROID_Y_ANGLE_PERIOD = 800;
const double TOROID_Z_ANGLE_PERIOD = 0;

const double MOBIUS_FOV_WIDTH = 2.5;
const double MOBIUS_FOV_HEIGHT = 2.5;
const double MOBIUS_FOV_Z = 3;

const double MOBIUS_U_INTERVAL_START = 0;
const double MOBIUS_U_INTERVAL_END = 2 * M_PI;
const double MOBIUS_U_INTERVAL_STEP = 0.05;

const double MOBIUS_V_INTERVAL_START = -1;
const double MOBIUS_V_INTERVAL_END = 1;
const double MOBIUS_V_INTERVAL_STEP = 0.02;

const int MOBIUS_U_STEPS_COUNT = (MOBIUS_U_INTERVAL_END - MOBIUS_U_INTERVAL_START) / MOBIUS_U_INTERVAL_STEP;
const int MOBIUS_V_STEPS_COUNT = (MOBIUS_V_INTERVAL_END - MOBIUS_V_INTERVAL_START) / MOBIUS_V_INTERVAL_STEP;
const int MOBIUS_PRECOMPUTED_DATA_COUNT = 4;

const double MOBIUS_X_ANGLE_PERIOD = 2000;
const double MOBIUS_Y_ANGLE_PERIOD = 1800;
const double MOBIUS_Z_ANGLE_PERIOD = 0;

const double WHIRLIGIG_FOV_WIDTH = 5;
const double WHIRLIGIG_FOV_HEIGHT = 5;
const double WHIRLIGIG_FOV_Z = 3;

const double WHIRLIGIG_U_INTERVAL_START = 0;
const double WHIRLIGIG_U_INTERVAL_END = 2 * M_PI;
const double WHIRLIGIG_U_INTERVAL_STEP = 0.05;

const double WHIRLIGIG_V_INTERVAL_START = 0;
const double WHIRLIGIG_V_INTERVAL_END = M_PI;
const double WHIRLIGIG_V_INTERVAL_STEP = 0.05;

const int WHIRLIGIG_U_STEPS_COUNT = (WHIRLIGIG_U_INTERVAL_END - WHIRLIGIG_U_INTERVAL_START) / WHIRLIGIG_U_INTERVAL_STEP;
const int WHIRLIGIG_V_STEPS_COUNT = (WHIRLIGIG_V_INTERVAL_END - WHIRLIGIG_V_INTERVAL_START) / WHIRLIGIG_V_INTERVAL_STEP;
const int WHIRLIGIG_PRECOMPUTED_DATA_COUNT = 6;

const double WHIRLIGIG_X_ANGLE_PERIOD = 1500;
const double WHIRLIGIG_Y_ANGLE_PERIOD = 1200;
const double WHIRLIGIG_Z_ANGLE_PERIOD = 0;

enum SurfaceType { Toroid,
                   Mobius,
                   Whirligig };

const SurfaceType SURFACE = Toroid;

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

  if (length == 0) {
    versor.x = vector.x / length;
    versor.y = vector.y / length;
    versor.z = vector.z / length;

    return versor;
  }

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

  double squaredNormalMagnitude = pow(magnitude(normal), 2);

  if (angleBetweenVectors(vector, normal) > M_PI / 2 || squaredNormalMagnitude == 0) {
    reflection.x = 0;
    reflection.y = 0;
    reflection.z = 0;

    return reflection;
  }

  double vectorAndNormalScalarProduct = scalarProduct(vector, normal);
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

double toroidalSurfaceX(double R, double r, double cosU, double cosV) {
  return (R + r * cosU) * cosV;
}

double toroidalSurfaceY(double R, double r, double cosU, double sinV) {
  return (R + r * cosU) * sinV;
}

double toroidalSurfaceZ(double R, double r, double sinU) {
  return r * sinU;
}

Vector3D computeToroidalSurfacePoint(double R, double r, double cosU, double cosV, double sinU, double sinV) {
  struct Vector3D point;

  point.x = toroidalSurfaceX(R, r, cosU, cosV);
  point.y = toroidalSurfaceY(R, r, cosU, sinV);
  point.z = toroidalSurfaceZ(R, r, sinU);

  return point;
}

/* Toroidal normal vectors */

double toroidalSurfaceNormalX(double R, double r, double cosU, double cosV) {
  return (R + r * cosU) * cosV * r * cosU;
}

double toroidalSurfaceNormalY(double R, double r, double cosU, double sinV) {
  return (R + r * cosU) * sinV * r * cosU;
}

double toroidalSurfaceNormalZ(double R, double r, double cosU, double sinU) {
  return r * sinU * (R + r * cosU);
}

Vector3D computeToroidalSurfaceNormalVersor(double R, double r, double cosU, double cosV, double sinU, double sinV) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = toroidalSurfaceNormalX(R, r, cosU, cosV);
  vector.y = toroidalSurfaceNormalY(R, r, cosU, sinV);
  vector.z = toroidalSurfaceNormalZ(R, r, cosU, sinU);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

void precomputeToroidalData(double*** precomputedDataMatrix) {
  for (int uIndex = 0; uIndex < TOROID_U_STEPS_COUNT; uIndex += 1) {
    double u = TOROID_U_INTERVAL_START + (double)uIndex * TOROID_U_INTERVAL_STEP;

    for (int vIndex = 0; vIndex < TOROID_V_STEPS_COUNT; vIndex += 1) {
      double v = TOROID_V_INTERVAL_START + (double)vIndex * TOROID_V_INTERVAL_STEP;

      precomputedDataMatrix[uIndex][vIndex][0] = cos(u);
      precomputedDataMatrix[uIndex][vIndex][1] = cos(v);
      precomputedDataMatrix[uIndex][vIndex][2] = sin(u);
      precomputedDataMatrix[uIndex][vIndex][3] = sin(v);
    }
  }
}

/* Möbius Strip surface  */

double mobiusSurfaceX(double v, double cosU, double cosUHalf) {
  return (1 + (v / 2) * cosUHalf) * cosU;
}

double mobiusSurfaceY(double v, double sinU, double cosUHalf) {
  return (1 + (v / 2) * cosUHalf) * sinU;
}

double mobiusSurfaceZ(double v, double sinUHalf) {
  return (v / 2) * sinUHalf;
}

Vector3D computeMobiusSurfacePoint(double v, double cosU, double sinU, double cosUHalf, double sinUHalf) {
  struct Vector3D point;

  point.x = mobiusSurfaceX(v, cosU, cosUHalf);
  point.y = mobiusSurfaceY(v, sinU, cosUHalf);
  point.z = mobiusSurfaceZ(v, sinUHalf);

  return point;
}

/* Möbius Strip normal vectors */

double mobiusSurfaceNormalX(double v, double cosU, double sinU, double cosUHalf, double sinUHalf) {
  return (1 / 2) * cosU * sinUHalf - (v / 4) * cosUHalf * sinUHalf * sinU - (v / 8) * sinU;
}

double mobiusSurfaceNormalY(double v, double cosU, double sinU, double cosUHalf, double sinUHalf) {
  return (1 / 2) * sinU * sinUHalf + (v / 4) * cosUHalf * sinUHalf * sinU + (v / 8) * cosU;
}

double mobiusSurfaceNormalZ(double v, double cosUHalf) {
  return -(1 / 2) * cosUHalf * (1 + (v / 2) * cosUHalf);
}

Vector3D computeMobiusSurfaceNormalVersor(double v, double cosU, double sinU, double cosUHalf, double sinUHalf) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = mobiusSurfaceNormalX(v, cosU, sinU, cosUHalf, sinUHalf);
  vector.y = mobiusSurfaceNormalY(v, cosU, sinU, cosUHalf, sinUHalf);
  vector.z = mobiusSurfaceNormalZ(v, cosUHalf);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

void precomputeMobiusData(double*** precomputedDataMatrix) {
  for (int uIndex = 0; uIndex < MOBIUS_U_STEPS_COUNT; uIndex += 1) {
    for (int vIndex = 0; vIndex < MOBIUS_V_STEPS_COUNT; vIndex += 1) {
      double u = MOBIUS_U_INTERVAL_START + (double)uIndex * MOBIUS_U_INTERVAL_STEP;
      double v = MOBIUS_V_INTERVAL_START + (double)vIndex * MOBIUS_V_INTERVAL_STEP;

      precomputedDataMatrix[uIndex][vIndex][0] = cos(u);
      precomputedDataMatrix[uIndex][vIndex][1] = sin(u);
      precomputedDataMatrix[uIndex][vIndex][2] = cos(u / 2);
      precomputedDataMatrix[uIndex][vIndex][3] = sin(u / 2);
    }
  }
}

/* Whirligig surface (from https://www.mathworks.com/help/matlab/ref/fsurf.html) */

double whirligigSurfaceRadius(double sin7U5V) {
  return 2 + sin7U5V;
}

double whirligigSurfaceX(double cosU, double sinV, double sin7U5V) {
  return whirligigSurfaceRadius(sin7U5V) * sinV * cosU;
}

double whirligigSurfaceY(double sinU, double sinV, double sin7U5V) {
  return whirligigSurfaceRadius(sin7U5V) * sinV * sinU;
}

double whirligigSurfaceZ(double cosV, double sin7U5V) {
  return whirligigSurfaceRadius(sin7U5V) * cosV;
}

Vector3D computeWhirligigSurfacePoint(double cosU, double sinU, double cosV, double sinV, double sin7U5V) {
  struct Vector3D point;

  point.x = whirligigSurfaceX(cosU, sinV, sin7U5V);
  point.y = whirligigSurfaceY(sinU, sinV, sin7U5V);
  point.z = whirligigSurfaceZ(cosV, sin7U5V);

  return point;
}

/* Whirligig normal vectors */

double whirligigSurfaceNormalX(double cosU, double sinU, double cosV, double sinV, double cos7U5V, double sin7U5V) {
  return sinV * (7 * sinU * cos7U5V + cosU * (sin7U5V + 2)) * (5 * cosV * cos7U5V - sinV * (sin7U5V + 2)) - 7 * sinU * cosV * cos7U5V * (5 * sinV * cos7U5V + cosV * (sin7U5V + 2));
}

double whirligigSurfaceNormalY(double cosU, double sinU, double cosV, double sinV, double cos7U5V, double sin7U5V) {
  return 7 * cosU * cosV * cos7U5V * (5 * sinV * cos7U5V + cosV * (sin7U5V + 2)) - sinV * (7 * cosU * cos7U5V - sinU * (sin7U5V + 2)) * (5 * cosV * cos7U5V - sinV * (sin7U5V + 2));
}

double whirligigSurfaceNormalZ(double cosU, double sinU, double cosV, double sinV, double cos7U5V, double sin7U5V) {
  return sinU * sinV * (5 * sinV * cos7U5V + cosV * (sin7U5V + 2)) * (7 * cosU * cos7U5V - sinU * (sin7U5V + 2)) - cosU * sinV * (7 * sinU * cos7U5V + cosU * (sin7U5V + 2)) * (5 * sinV * cos7U5V + cosV * (sin7U5V + 2));
}

Vector3D computeWhirligigSurfaceNormalVersor(double cosU, double sinU, double cosV, double sinV, double cos7U5V, double sin7U5V) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = whirligigSurfaceNormalX(cosU, sinU, cosV, sinV, cos7U5V, sin7U5V);
  vector.y = whirligigSurfaceNormalY(cosU, sinU, cosV, sinV, cos7U5V, sin7U5V);
  vector.z = whirligigSurfaceNormalZ(cosU, sinU, cosV, sinV, cos7U5V, sin7U5V);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

void precomputeWhirligigData(double*** precomputedDataMatrix) {
  for (int uIndex = 0; uIndex < WHIRLIGIG_U_STEPS_COUNT; uIndex += 1) {
    for (int vIndex = 0; vIndex < WHIRLIGIG_V_STEPS_COUNT; vIndex += 1) {
      double u = WHIRLIGIG_U_INTERVAL_START + (double)uIndex * WHIRLIGIG_U_INTERVAL_STEP;
      double v = WHIRLIGIG_V_INTERVAL_START + (double)vIndex * WHIRLIGIG_V_INTERVAL_STEP;

      precomputedDataMatrix[uIndex][vIndex][0] = cos(u);
      precomputedDataMatrix[uIndex][vIndex][1] = sin(u);
      precomputedDataMatrix[uIndex][vIndex][2] = cos(v);
      precomputedDataMatrix[uIndex][vIndex][3] = sin(v);
      precomputedDataMatrix[uIndex][vIndex][4] = cos(7 * u + 5 * v);
      precomputedDataMatrix[uIndex][vIndex][5] = sin(7 * u + 5 * v);
    }
  }
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
    double lightIntensityMatrix[][CANVAS_WIDTH]) {
  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      depthBufferMatrix[y][x] = -DBL_MAX;
      lightIntensityMatrix[y][x] = 0;
    }
  }
}

void render(
    SurfaceType surfaceType,
    Vector3D cameraPoint,
    Vector3D lightPoint,
    double depthBufferMatrix[][CANVAS_WIDTH],
    double lightIntensityMatrix[][CANVAS_WIDTH],
    double*** precomputedDataMatrix,
    Matrix3X3 rotationMatrix) {
  double fovWidth,
      fovHeight,
      fovZ,
      uIntervalStart,
      uIntervalEnd,
      uIntervalStep,
      uIntervalStepsCount,
      vIntervalStart,
      vIntervalEnd,
      vIntervalStep,
      vIntervalStepsCount;

  switch (surfaceType) {
    case Toroid:
      fovWidth = TOROID_FOV_WIDTH;
      fovHeight = TOROID_FOV_HEIGHT;
      fovZ = TOROID_FOV_Z;
      uIntervalStart = TOROID_U_INTERVAL_START;
      uIntervalEnd = TOROID_U_INTERVAL_END;
      uIntervalStep = TOROID_U_INTERVAL_STEP;
      uIntervalStepsCount = TOROID_U_STEPS_COUNT;
      vIntervalStart = TOROID_V_INTERVAL_START;
      vIntervalEnd = TOROID_V_INTERVAL_END;
      vIntervalStep = TOROID_V_INTERVAL_STEP;
      vIntervalStepsCount = TOROID_V_STEPS_COUNT;

      break;
    case Mobius:
      fovWidth = MOBIUS_FOV_WIDTH;
      fovHeight = MOBIUS_FOV_HEIGHT;
      fovZ = MOBIUS_FOV_Z;
      uIntervalStart = MOBIUS_U_INTERVAL_START;
      uIntervalEnd = MOBIUS_U_INTERVAL_END;
      uIntervalStep = MOBIUS_U_INTERVAL_STEP;
      uIntervalStepsCount = MOBIUS_U_STEPS_COUNT;
      vIntervalStart = MOBIUS_V_INTERVAL_START;
      vIntervalEnd = MOBIUS_V_INTERVAL_END;
      vIntervalStep = MOBIUS_V_INTERVAL_STEP;
      vIntervalStepsCount = MOBIUS_V_STEPS_COUNT;

      break;
    case Whirligig:
      fovWidth = WHIRLIGIG_FOV_WIDTH;
      fovHeight = WHIRLIGIG_FOV_HEIGHT;
      fovZ = WHIRLIGIG_FOV_Z;
      uIntervalStart = WHIRLIGIG_U_INTERVAL_START;
      uIntervalEnd = WHIRLIGIG_U_INTERVAL_END;
      uIntervalStep = WHIRLIGIG_U_INTERVAL_STEP;
      uIntervalStepsCount = WHIRLIGIG_U_STEPS_COUNT;
      vIntervalStart = WHIRLIGIG_V_INTERVAL_START;
      vIntervalEnd = WHIRLIGIG_V_INTERVAL_END;
      vIntervalStep = WHIRLIGIG_V_INTERVAL_STEP;
      vIntervalStepsCount = WHIRLIGIG_V_STEPS_COUNT;

      break;
  }

  for (int uIndex = 0; uIndex < uIntervalStepsCount; uIndex++) {
    double u = uIntervalStart + (double)uIndex * uIntervalStep;

    for (int vIndex = 0; vIndex < vIntervalStepsCount; vIndex++) {
      double v = vIntervalStart + (double)vIndex * vIntervalStep;

      /* Computing current surface point and normal vector */

      struct Vector3D surfacePoint, surfaceNormal;

      switch (surfaceType) {
        case Toroid:
          surfacePoint = computeToroidalSurfacePoint(TOROID_EXTERNAL_RADIUS,
                                                     TOROID_INTERNAL_RADIUS,
                                                     precomputedDataMatrix[uIndex][vIndex][0],
                                                     precomputedDataMatrix[uIndex][vIndex][1],
                                                     precomputedDataMatrix[uIndex][vIndex][2],
                                                     precomputedDataMatrix[uIndex][vIndex][3]);

          surfaceNormal = computeToroidalSurfaceNormalVersor(TOROID_EXTERNAL_RADIUS,
                                                             TOROID_INTERNAL_RADIUS,
                                                             precomputedDataMatrix[uIndex][vIndex][0],
                                                             precomputedDataMatrix[uIndex][vIndex][1],
                                                             precomputedDataMatrix[uIndex][vIndex][2],
                                                             precomputedDataMatrix[uIndex][vIndex][3]);

          break;
        case Mobius:
          surfacePoint = computeMobiusSurfacePoint(v,
                                                   precomputedDataMatrix[uIndex][vIndex][0],
                                                   precomputedDataMatrix[uIndex][vIndex][1],
                                                   precomputedDataMatrix[uIndex][vIndex][2],
                                                   precomputedDataMatrix[uIndex][vIndex][3]);

          surfaceNormal = computeMobiusSurfaceNormalVersor(v,
                                                           precomputedDataMatrix[uIndex][vIndex][0],
                                                           precomputedDataMatrix[uIndex][vIndex][1],
                                                           precomputedDataMatrix[uIndex][vIndex][2],
                                                           precomputedDataMatrix[uIndex][vIndex][3]);

          break;
        case Whirligig:
          surfacePoint = computeWhirligigSurfacePoint(precomputedDataMatrix[uIndex][vIndex][0],
                                                      precomputedDataMatrix[uIndex][vIndex][1],
                                                      precomputedDataMatrix[uIndex][vIndex][2],
                                                      precomputedDataMatrix[uIndex][vIndex][3],
                                                      precomputedDataMatrix[uIndex][vIndex][5]);

          surfaceNormal = computeWhirligigSurfaceNormalVersor(precomputedDataMatrix[uIndex][vIndex][0],
                                                              precomputedDataMatrix[uIndex][vIndex][1],
                                                              precomputedDataMatrix[uIndex][vIndex][2],
                                                              precomputedDataMatrix[uIndex][vIndex][3],
                                                              precomputedDataMatrix[uIndex][vIndex][4],
                                                              precomputedDataMatrix[uIndex][vIndex][5]);

          break;
      }

      /* Rotating surface point and normal vector */

      struct Vector3D rotatedSurfacePoint = matrixVector3DProduct(rotationMatrix, surfacePoint);
      struct Vector3D rotatedSurfaceNormal = matrixVector3DProduct(rotationMatrix, surfaceNormal);

      /* Computing point projection on the canvas */

      struct Vector3D surfaceToCameraVector = computeVector(rotatedSurfacePoint, cameraPoint);
      struct Vector3D surfaceToCameraVersor = computeVersor(surfaceToCameraVector);
      struct Vector3D projection = projectOnPlane(cameraPoint, surfaceToCameraVersor, fovZ);

      struct Vector2D matrixCoordinates = projectionToCanvasMatrix(
          projection,
          fovWidth,
          fovHeight,
          CANVAS_HEIGHT,
          CANVAS_WIDTH);

      /* Computing reflecting light intensity */

      if (rotatedSurfacePoint.z > depthBufferMatrix[matrixCoordinates.y][matrixCoordinates.x]) {
        struct Vector3D lightVector = computeVector(rotatedSurfacePoint, lightPoint);
        struct Vector3D reflectionVector = reflectAcrossNormal(lightVector, rotatedSurfaceNormal);

        double reflectionAngle = angleBetweenVectors(surfaceToCameraVector, reflectionVector);
        double lightIntensity = 1 - (reflectionAngle / M_PI);

        lightIntensityMatrix[matrixCoordinates.y][matrixCoordinates.x] = lightIntensity;
        depthBufferMatrix[matrixCoordinates.y][matrixCoordinates.x] = rotatedSurfacePoint.z;
      }
    }
  }
}

void print(
    double depthBufferMatrix[][CANVAS_WIDTH],
    double lightIntensityMatrix[][CANVAS_WIDTH]) {
  double maxLightIntensity = 0;

  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      if (lightIntensityMatrix[y][x] > maxLightIntensity) {
        maxLightIntensity = lightIntensityMatrix[y][x];
      }
    }
  }

  cout << "\x1b[H";

  for (int y = 0; y < CANVAS_HEIGHT; y++) {
    for (int x = 0; x < CANVAS_WIDTH; x++) {
      if (depthBufferMatrix[y][x] == -DBL_MAX) {
        cout << ' ';
      } else {
        cout << computeIntensityCharacter(lightIntensityMatrix[y][x], maxLightIntensity);
      }
    }

    cout << endl;
  }
}

/* Main */

int main() {
  double depthBufferMatrix[CANVAS_HEIGHT][CANVAS_WIDTH];
  double lightIntensityMatrix[CANVAS_HEIGHT][CANVAS_WIDTH];
  double*** precomputedDataMatrix;

  /* Camera point */

  struct Vector3D cameraPoint;

  cameraPoint.x = 0;
  cameraPoint.y = 0;
  cameraPoint.z = CAMERA_Z;

  /* Light point */

  struct Vector3D lightPoint;

  lightPoint.x = LIGHT_X;
  lightPoint.y = LIGHT_Y;
  lightPoint.z = LIGHT_Z;

  /* Rotation periods */

  double xAnglePeriod, yAnglePeriod, zAnglePeriod;

  switch (SURFACE) {
    case Toroid:
      xAnglePeriod = TOROID_X_ANGLE_PERIOD;
      yAnglePeriod = TOROID_Y_ANGLE_PERIOD;
      zAnglePeriod = TOROID_Z_ANGLE_PERIOD;

      break;
    case Mobius:
      xAnglePeriod = MOBIUS_X_ANGLE_PERIOD;
      yAnglePeriod = MOBIUS_Y_ANGLE_PERIOD;
      zAnglePeriod = MOBIUS_Z_ANGLE_PERIOD;

      break;

    case Whirligig:
      xAnglePeriod = WHIRLIGIG_X_ANGLE_PERIOD;
      yAnglePeriod = WHIRLIGIG_Y_ANGLE_PERIOD;
      zAnglePeriod = WHIRLIGIG_Z_ANGLE_PERIOD;

      break;
  }

  /* Precomputing data */

  int uStepsCount, vStepsCount, precomputedDataCount;

  switch (SURFACE) {
    case Toroid:
      uStepsCount = TOROID_U_STEPS_COUNT;
      vStepsCount = TOROID_V_STEPS_COUNT;
      precomputedDataCount = TOROID_PRECOMPUTED_DATA_COUNT;

      break;
    case Mobius:
      uStepsCount = MOBIUS_U_STEPS_COUNT;
      vStepsCount = MOBIUS_V_STEPS_COUNT;
      precomputedDataCount = MOBIUS_PRECOMPUTED_DATA_COUNT;

      break;
    case Whirligig:
      uStepsCount = WHIRLIGIG_U_STEPS_COUNT;
      vStepsCount = WHIRLIGIG_V_STEPS_COUNT;
      precomputedDataCount = WHIRLIGIG_PRECOMPUTED_DATA_COUNT;

      break;
  }

  precomputedDataMatrix = new double**[uStepsCount];

  for (int uIndex = 0; uIndex < uStepsCount; uIndex++) {
    precomputedDataMatrix[uIndex] = new double*[vStepsCount];

    for (int vIndex = 0; vIndex < vStepsCount; vIndex++) {
      precomputedDataMatrix[uIndex][vIndex] = new double[precomputedDataCount];
    }
  }

  switch (SURFACE) {
    case Toroid:
      precomputeToroidalData(precomputedDataMatrix);

      break;
    case Mobius:
      precomputeMobiusData(precomputedDataMatrix);

      break;
    case Whirligig:
      precomputeWhirligigData(precomputedDataMatrix);

      break;
  }

  /* Rendering cycle */

  for (double i = 0;; i += 1) {
    // /* Rotation angles */

    double xAngle = xAnglePeriod != 0 ? i * 2 * M_PI / xAnglePeriod : 0;
    double yAngle = yAnglePeriod != 0 ? i * 2 * M_PI / yAnglePeriod : 0;
    double zAngle = zAnglePeriod != 0 ? i * 2 * M_PI / zAnglePeriod : 0;

    struct Matrix3X3 rotationMatrix = computeRotationMatrix(xAngle, yAngle, zAngle);

    // /* Resetting matrices */

    resetMatrices(depthBufferMatrix, lightIntensityMatrix);

    /* Rendering */

    render(
        SURFACE,
        cameraPoint,
        lightPoint,
        depthBufferMatrix,
        lightIntensityMatrix,
        precomputedDataMatrix,
        rotationMatrix);

    /* Printing */

    print(depthBufferMatrix, lightIntensityMatrix);
  }

  /* Exiting */

  return 0;
}
