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

const double MOBIUS_X_ANGLE_PERIOD = 2000;
const double MOBIUS_Y_ANGLE_PERIOD = 1800;
const double MOBIUS_Z_ANGLE_PERIOD = 0;

const double WHIRLIGIG_FOV_WIDTH = 6;
const double WHIRLIGIG_FOV_HEIGHT = 6;
const double WHIRLIGIG_FOV_Z = 3;

const double WHIRLIGIG_U_INTERVAL_START = 0;
const double WHIRLIGIG_U_INTERVAL_END = 2 * M_PI;
const double WHIRLIGIG_U_INTERVAL_STEP = 0.05;

const double WHIRLIGIG_V_INTERVAL_START = 0;
const double WHIRLIGIG_V_INTERVAL_END = M_PI;
const double WHIRLIGIG_V_INTERVAL_STEP = 0.05;

const double WHIRLIGIG_X_ANGLE_PERIOD = 1000;
const double WHIRLIGIG_Y_ANGLE_PERIOD = 800;
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

double mobiusSurfaceNormalX(double u, double v) {
  return (1 / 2) * cos(u) * sin(u / 2) - (v / 4) * cos(u / 2) * sin(u / 2) * sin(u) - (v / 8) * sin(u);
}

double mobiusSurfaceNormalY(double u, double v) {
  return (1 / 2) * sin(u) * sin(u / 2) + (v / 4) * cos(u / 2) * sin(u / 2) * sin(u) + (v / 8) * cos(u);
}

double mobiusSurfaceNormalZ(double u, double v) {
  return -(1 / 2) * cos(u / 2) * (1 + (v / 2) * cos(u / 2));
}

Vector3D computeMobiusSurfaceNormalVersor(double u, double v) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = mobiusSurfaceNormalX(u, v);
  vector.y = mobiusSurfaceNormalY(u, v);
  vector.z = mobiusSurfaceNormalZ(u, v);

  double vectorMagnitude = magnitude(vector);

  versor.x = vector.x / vectorMagnitude;
  versor.y = vector.y / vectorMagnitude;
  versor.z = vector.z / vectorMagnitude;

  return versor;
}

/* Whirligig surface (from https://www.mathworks.com/help/matlab/ref/fsurf.html) */

double whirligigSurfaceRadius(double u, double v) {
  return 2 + sin(7 * u + 5 * v);
}

double whirligigSurfaceX(double u, double v) {
  return whirligigSurfaceRadius(u, v) * sin(v) * cos(u);
}

double whirligigSurfaceY(double u, double v) {
  return whirligigSurfaceRadius(u, v) * sin(v) * sin(u);
}

double whirligigSurfaceZ(double u, double v) {
  return whirligigSurfaceRadius(u, v) * cos(v);
}

Vector3D computeWhirligigSurfacePoint(double u, double v) {
  struct Vector3D point;

  point.x = whirligigSurfaceX(u, v);
  point.y = whirligigSurfaceY(u, v);
  point.z = whirligigSurfaceZ(u, v);

  return point;
}

/* Whirligig normal vectors */

double whirligigSurfaceNormalX(double u, double v) {
  return sin(v) * (7 * sin(u) * cos(7 * u + 5 * v) + cos(u) * (sin(7 * u + 5 * v) + 2)) * (5 * cos(v) * cos(7 * u + 5 * v) - sin(v) * (sin(7 * u + 5 * v) + 2)) - 7 * sin(u) * cos(v) * cos(7 * u + 5 * v) * (5 * sin(v) * cos(7 * u + 5 * v) + cos(v) * (sin(7 * u + 5 * v) + 2));
}

double whirligigSurfaceNormalY(double u, double v) {
  return 7 * cos(u) * cos(v) * cos(7 * u + 5 * v) * (5 * sin(v) * cos(7 * u + 5 * v) + cos(v) * (sin(7 * u + 5 * v) + 2)) - sin(v) * (7 * cos(u) * cos(7 * u + 5 * v) - sin(u) * (sin(7 * u + 5 * v) + 2)) * (5 * cos(v) * cos(7 * u + 5 * v) - sin(v) * (sin(7 * u + 5 * v) + 2));
}

double whirligigSurfaceNormalZ(double u, double v) {
  return sin(u) * sin(v) * (5 * sin(v) * cos(7 * u + 5 * v) + cos(v) * (sin(7 * u + 5 * v) + 2)) * (7 * cos(u) * cos(7 * u + 5 * v) - sin(u) * (sin(7 * u + 5 * v) + 2)) - cos(u) * sin(v) * (7 * sin(u) * cos(7 * u + 5 * v) + cos(u) * (sin(7 * u + 5 * v) + 2)) * (5 * sin(v) * cos(7 * u + 5 * v) + cos(v) * (sin(7 * u + 5 * v) + 2));
}

Vector3D computeWhirligigSurfaceNormalVersor(double u, double v) {
  struct Vector3D vector;
  struct Vector3D versor;

  vector.x = whirligigSurfaceNormalX(u, v);
  vector.y = whirligigSurfaceNormalY(u, v);
  vector.z = whirligigSurfaceNormalZ(u, v);

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
    SurfaceType surfaceType,
    Vector3D cameraPoint,
    Vector3D lightPoint,
    double depthBufferMatrix[][CANVAS_WIDTH],
    double absoluteIntensityMatrix[][CANVAS_WIDTH],
    double xAngle = 0,
    double yAngle = 0,
    double zAngle = 0) {
  struct Matrix3X3 rotationMatrix = computeRotationMatrix(xAngle, yAngle, zAngle);

  double fovWidth,
      fovHeight,
      fovZ,
      uIntervalStart,
      uIntervalEnd,
      uIntervalStep,
      vIntervalStart,
      vIntervalEnd,
      vIntervalStep;

  switch (surfaceType) {
    case Toroid:
      fovWidth = TOROID_FOV_WIDTH;
      fovHeight = TOROID_FOV_HEIGHT;
      fovZ = TOROID_FOV_Z;
      uIntervalStart = TOROID_U_INTERVAL_START;
      uIntervalEnd = TOROID_U_INTERVAL_END;
      uIntervalStep = TOROID_U_INTERVAL_STEP;
      vIntervalStart = TOROID_V_INTERVAL_START;
      vIntervalEnd = TOROID_V_INTERVAL_END;
      vIntervalStep = TOROID_V_INTERVAL_STEP;

      break;
    case Mobius:
      fovWidth = MOBIUS_FOV_WIDTH;
      fovHeight = MOBIUS_FOV_HEIGHT;
      fovZ = MOBIUS_FOV_Z;
      uIntervalStart = MOBIUS_U_INTERVAL_START;
      uIntervalEnd = MOBIUS_U_INTERVAL_END;
      uIntervalStep = MOBIUS_U_INTERVAL_STEP;
      vIntervalStart = MOBIUS_V_INTERVAL_START;
      vIntervalEnd = MOBIUS_V_INTERVAL_END;
      vIntervalStep = MOBIUS_V_INTERVAL_STEP;

      break;
    case Whirligig:
      fovWidth = WHIRLIGIG_FOV_WIDTH;
      fovHeight = WHIRLIGIG_FOV_HEIGHT;
      fovZ = WHIRLIGIG_FOV_Z;
      uIntervalStart = WHIRLIGIG_U_INTERVAL_START;
      uIntervalEnd = WHIRLIGIG_U_INTERVAL_END;
      uIntervalStep = WHIRLIGIG_U_INTERVAL_STEP;
      vIntervalStart = WHIRLIGIG_V_INTERVAL_START;
      vIntervalEnd = WHIRLIGIG_V_INTERVAL_END;
      vIntervalStep = WHIRLIGIG_V_INTERVAL_STEP;

      break;
  }

  for (double u = uIntervalStart; u < uIntervalEnd; u = u + uIntervalStep) {
    for (double v = vIntervalStart; v < vIntervalEnd; v = v + vIntervalStep) {
      struct Vector3D surfacePoint, surfaceNormal;

      switch (surfaceType) {
        case Toroid:
          surfacePoint = computeToroidalSurfacePoint(TOROID_EXTERNAL_RADIUS, TOROID_INTERNAL_RADIUS, u, v);
          surfaceNormal = computeToroidalSurfaceNormalVersor(TOROID_EXTERNAL_RADIUS, TOROID_INTERNAL_RADIUS, u, v);

          break;
        case Mobius:
          surfacePoint = computeMobiusSurfacePoint(u, v);
          surfaceNormal = computeMobiusSurfaceNormalVersor(u, v);

          break;
        case Whirligig:
          surfacePoint = computeWhirligigSurfacePoint(u, v);
          surfaceNormal = computeWhirligigSurfaceNormalVersor(u, v);

          break;
      }

      struct Vector3D rotatedSurfacePoint = matrixVector3DProduct(rotationMatrix, surfacePoint);
      struct Vector3D rotatedSurfaceNormal = matrixVector3DProduct(rotationMatrix, surfaceNormal);

      struct Vector3D surfaceToCameraVector = computeVector(rotatedSurfacePoint, cameraPoint);
      struct Vector3D surfaceToCameraVersor = computeVersor(surfaceToCameraVector);
      struct Vector3D lightVector = computeVector(rotatedSurfacePoint, lightPoint);
      struct Vector3D reflectionVector = reflectAcrossNormal(lightVector, rotatedSurfaceNormal);
      struct Vector3D projection = projectOnPlane(cameraPoint, surfaceToCameraVersor, fovZ);

      double reflectionAngle = angleBetweenVectors(surfaceToCameraVector, reflectionVector);
      double lightIntensity = 1 - (reflectionAngle / M_PI);

      struct Vector2D matrixCoordinates = projectionToCanvasMatrix(
          projection,
          fovWidth,
          fovHeight,
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

  /* Camera point */

  struct Vector3D cameraPoint;

  cameraPoint.x = 0;
  cameraPoint.y = 0;
  cameraPoint.z = CAMERA_Z;

  for (double i = 0;; i += 1) {
    /* Light vector */

    struct Vector3D lightPoint;

    lightPoint.x = LIGHT_X;
    lightPoint.y = LIGHT_Y;
    lightPoint.z = LIGHT_Z;

    // /* Rotation angles */

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

    double xAngle = xAnglePeriod != 0 ? i * 2 * M_PI / xAnglePeriod : 0;
    double yAngle = yAnglePeriod != 0 ? i * 2 * M_PI / yAnglePeriod : 0;
    double zAngle = zAnglePeriod != 0 ? i * 2 * M_PI / zAnglePeriod : 0;

    // /* Resetting matrices */

    resetMatrices(depthBufferMatrix, absoluteIntensityMatrix);

    /* Rendering */

    render(
        SURFACE,
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
