#pragma once

#include <math.h>
#include <ostream>
#include <iostream>

class Math
{
public:
  static constexpr double Clamp(double, double = 0.0, double = 1.0);

  // Minimum and maximum
  static constexpr double Min(double, double);
  static constexpr double Max(double, double);
  static constexpr double Min(double, double, double);
  static constexpr double Max(double, double, double);

  static constexpr double DegreeToRadian(double);
  static constexpr double RadianToDegree(double);
};

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline constexpr double Math::Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Min(double a, double b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Max(double a, double b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Max(double a, double b, double c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Min(double a, double b, double c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline constexpr double Math::DegreeToRadian(double a)
{
  return a * 3.14159265358979323846 / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline constexpr double Math::RadianToDegree(double a)
{
  return a * 180.0 / 3.14159265358979323846;
}

// Class
class Vector
{
protected:
  double c[3]; //!< Components.
public:
  //! Empty
  Vector() {}

  explicit Vector(double);
  explicit Vector(double, double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vector operator+ () const;
  Vector operator- () const;

  // Assignment operators
  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (const Vector&);
  Vector& operator/= (const Vector&);
  Vector& operator*= (double);
  Vector& operator/= (double);

  // Binary operators
  friend int operator> (const Vector&, const Vector&);
  friend int operator< (const Vector&, const Vector&);

  friend int operator>= (const Vector&, const Vector&);
  friend int operator<= (const Vector&, const Vector&);

  // Binary operators
  friend Vector operator+ (const Vector&, const Vector&);
  friend Vector operator- (const Vector&, const Vector&);

  friend constexpr double operator* (const Vector&, const Vector&);

  friend Vector operator* (const Vector&, double);
  friend Vector operator* (double, const Vector&);
  friend Vector operator/ (const Vector&, double);

  friend Vector operator/ (const Vector&, const Vector&);

  // Boolean functions
  friend int operator==(const Vector&, const Vector&);
  friend int operator!=(const Vector&, const Vector&);

  // Norm
  friend double Norm(const Vector&);
  friend double SquaredNorm(const Vector&);

  friend void Normalize(Vector&);
  friend Vector Normalized(const Vector&);

  // Compare functions
  static Vector Min(const Vector&, const Vector&);
  static Vector Max(const Vector&, const Vector&);

  // Abs
  friend Vector Abs(const Vector&);

  // Orthogonal and orthonormal vectors
  Vector Orthogonal() const;
  void Orthonormal(Vector&, Vector&) const;

  friend Vector Lerp(const Vector&, const Vector&, double);
  static Vector Bilinear(const Vector&, const Vector&, const Vector&, const Vector&, double, double);

  // Scale
  Vector Scaled(const Vector&) const;
  Vector Inverse() const;

  friend std::ostream& operator<<(std::ostream&, const Vector&);




public:
  static const Vector Null; //!< Null vector.
  static const Vector X; //!< Vector(1,0,0).
  static const Vector Y; //!< Vector(0,1,0).
  static const Vector Z; //!< Vector(0,0,1).
};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vector::Vector(double a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vector::Vector(double a, double b, double c)
{
  Vector::c[0] = a;
  Vector::c[1] = b;
  Vector::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline double& Vector::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline double Vector::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vector Vector::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vector Vector::operator- () const
{
  return Vector(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vector& Vector::operator+= (const Vector& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vector& Vector::operator-= (const Vector& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vector& Vector::operator*= (double a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Scale a vector.
\param a Scaling vector.
*/
inline Vector Vector::Scaled(const Vector& a) const
{
  return Vector(c[0] * a[0], c[1] * a[1], c[2] * a[2]);
}

/*!
\brief Inverse of a vector.

This function inverses the components of the vector. This is the same as:
\code
Vector v=Vector(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vector Vector::Inverse() const
{
  return Vector(1.0 / c[0], 1.0 / c[1], 1.0 / c[2]);
}

//! Destructive division by a scalar.
inline Vector& Vector::operator/= (double a)
{
  c[0] /= a; c[1] /= a; c[2] /= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.

This is the same as Scale:
\code
Vector u(2.0,-1.0,1.0);
u=u.Scaled(Vector(3.0,1.0,2.0)); // u*=Vector(3.0,1.0,2.0);
\endcode
*/
inline Vector& Vector::operator*= (const Vector& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Destructively divide the components of a vector by another vector.
inline Vector& Vector::operator/= (const Vector& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1]; c[2] /= u.c[2];
  return *this;
}

//! Compare two vectors.
inline int operator> (const Vector& u, const Vector& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]) && (u.c[2] > v.c[2]));
}

//! Compare two vectors.
inline int operator< (const Vector& u, const Vector& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]) && (u.c[2] < v.c[2]));
}

//! Overloaded
inline int operator>= (const Vector& u, const Vector& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]) && (u.c[2] >= v.c[2]));
}

//! Overloaded
inline int operator<= (const Vector& u, const Vector& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]) && (u.c[2] <= v.c[2]));
}

//! Adds up two vectors.
inline Vector operator+ (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vector operator- (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

//! Scalar product.
inline constexpr double operator* (const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}

//! Right multiply by a scalar.
inline Vector operator* (const Vector& u, double a)
{
  return Vector(u.c[0] * a, u.c[1] * a, u.c[2] * a);
}

//! Left multiply by a scalar.
inline Vector operator* (double a, const Vector& v)
{
  return v * a;
}

//! Cross product.
inline Vector operator/ (const Vector& u, const Vector& v)
{
  return Vector(u.c[1] * v.c[2] - u.c[2] * v.c[1], u.c[2] * v.c[0] - u.c[0] * v.c[2], u.c[0] * v.c[1] - u.c[1] * v.c[0]);
}

//! Left multiply by a scalar
inline Vector operator/ (const Vector& u, double a)
{
  return Vector(u.c[0] / a, u.c[1] / a, u.c[2] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vector& u, const Vector& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline int operator!= (const Vector& u, const Vector& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a vector.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a vector instead.
\param u %Vector.
\sa SquaredNorm
*/
inline double Norm(const Vector& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vector.
\sa Norm
*/
inline double SquaredNorm(const Vector& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Return a normalized vector.

Compute the inverse of its norm and scale the components.

This function does not check if the vector is null.
\param u %Vector.
*/
inline Vector Normalized(const Vector& u)
{
  return u * (1.0 / Norm(u));
}

/*!
\brief Computes the absolute value of a vector.
\param u %Vector.
*/
inline Vector Abs(const Vector& u)
{
  return Vector(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

/*!
\brief Return a vector with coordinates set to the minimum coordinates
of the two argument vectors.
*/
inline Vector Vector::Min(const Vector& a, const Vector& b)
{
  return Vector(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1], a[2] < b[2] ? a[2] : b[2]);
}

/*!
\brief Return a vector with coordinates set to the maximum coordinates
of the two argument vectors.
*/
inline Vector Vector::Max(const Vector& a, const Vector& b)
{
  return Vector(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1], a[2] > b[2] ? a[2] : b[2]);
}

/*!
\brief Linear interpolation between two vectors.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vector Lerp(const Vector& a, const Vector& b, double t)
{
  return a + t * (b - a);
}

/*!
\brief Bi-linear interpolation between four vectors.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated vectors.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vector Vector::Bilinear(const Vector& a00, const Vector& a10, const Vector& a11, const Vector& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

class Matrix {
private:
    double tabMatrice[9];
public:
    Matrix();
    Matrix(const Matrix &m);
    Matrix(double a, double b, double c, double d, double e, double f, double g, double h, double i);
    ~Matrix();

    double& operator[](int i);
    double operator[](int i) const;

    Matrix operator+ (const Matrix &m) const;
    Matrix operator- (const Matrix &m) const;
    Matrix operator* (const Matrix &m) const;
    Matrix operator* (double x) const;
    Matrix operator/ (double x) const;

    Matrix& operator=(const Matrix &m);
    Matrix& operator+= (const Matrix &m);
    Matrix& operator-= (const Matrix &m);
    Matrix& operator*= (const Matrix &m);
    Matrix& operator*= (double x);
    Matrix& operator/= (double x);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);

    Matrix transpose();
    Matrix inverse();

    static const Matrix rot_x(const double &m);
    static const Matrix rot_y(const double &m);
    static const Matrix rot_z(const double &m);

    static const Matrix omth(const double &x, const double &y, const double &z);

    Vector operator* (const Vector &vec) const;
};


/*!
 * \brief Constructeur d'une matrice
 */
inline Matrix::Matrix() {
    for(int i=0; i<9; i++) {
        tabMatrice[i] = 0;
    }
}

/*!
 * \brief Constructeur par copie d'une matrice
 * \param m Matrix
 */
inline Matrix::Matrix(const Matrix &m) {
    for(int i=0; i<8; i++) {
        tabMatrice[i] = m[i];
    }
}

/*!
 * \brief Constructeur avec double en paramètre
 * \param a double
 * \param b double
 * \param c double
 * \param d double
 * \param e double
 * \param f double
 * \param g double
 * \param h double
 * \param i double
 */
inline Matrix::Matrix(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    tabMatrice[0] = a;
    tabMatrice[1] = b;
    tabMatrice[2] = c;
    tabMatrice[3] = d;
    tabMatrice[4] = e;
    tabMatrice[5] = f;
    tabMatrice[6] = g;
    tabMatrice[7] = h;
    tabMatrice[8] = i;
}

/*!
 * \brief Destructeur d'une matrice
 */
inline Matrix::~Matrix() {

}

/*!
 * \brief Surcharge de l'opérateur []
 * \param i entier pour savoir quel élément supprimer
 * \return entier valeur de la matrice
 */
inline double& Matrix::operator[](int i) {
  return tabMatrice[i];
}

/*!
 * \brief Surcharge de l'opérateur +
 * \param m Matrix
 * \return Matrix LATRIX additioné avec la matrice m
 */
inline Matrix Matrix::operator+(const Matrix &m) const {
  Matrix LATRIX;

  for(int i=0; i<9; i++) {
    LATRIX[i] = tabMatrice[i] + m[i];
  }
  return LATRIX;
}

/*!
 * \brief Surcharge de l'opérateur -
 * \param m Matrix
 * \return Matrix LATRIX soutrait avec la matrice m
 */
inline Matrix Matrix::operator-(const Matrix &m) const {
  Matrix LATRIX;

  for(int i=0; i<9; i++) {
    LATRIX[i] = tabMatrice[i] - m[i];
  }
  return LATRIX;
}

/*!
 * \brief Surcharge de l'opérateur *
 * \param m Matrix
 * \return Matrix LATRIX multiplié avec la matrice m
 */
inline Matrix Matrix::operator*(const Matrix &m) const {
  Matrix LATRIX;

  LATRIX[0]  = tabMatrice[0] * m[0] + tabMatrice[1] * m[3] + tabMatrice[2] * m[6];
  LATRIX[1]  = tabMatrice[0] * m[1] + tabMatrice[1] * m[4] + tabMatrice[2] * m[7];
  LATRIX[2]  = tabMatrice[0] * m[2] + tabMatrice[1] * m[5] + tabMatrice[2] * m[8];

  LATRIX[3]  = tabMatrice[3] * m[0] + tabMatrice[4] * m[3] + tabMatrice[5] * m[6];
  LATRIX[4]  = tabMatrice[3] * m[1] + tabMatrice[4] * m[4] + tabMatrice[5] * m[7];
  LATRIX[5]  = tabMatrice[3] * m[2] + tabMatrice[4] * m[5] + tabMatrice[5] * m[8];

  LATRIX[6]  = tabMatrice[6] * m[0] + tabMatrice[7] * m[3] + tabMatrice[8] * m[6];
  LATRIX[7]  = tabMatrice[6] * m[1] + tabMatrice[7] * m[4] + tabMatrice[8] * m[7];
  LATRIX[8]  = tabMatrice[6] * m[2] + tabMatrice[7] * m[5] + tabMatrice[8] * m[8];

  return LATRIX;
}

/*!
 * \brief Surcharge de l'opérateur *
 * \param x double
 * \return Matrix LATRIX multiplié avec une valeur
 */
inline Matrix Matrix::operator*(double x) const {
  Matrix LATRIX;

  for(int i=0; i<9; i++) {
    LATRIX[i] = tabMatrice[i] * x;
  }
  return LATRIX;
}

/*!
 * \brief Surcharge de l'opérateur /
 * \param x double
 * \return Matrix LATRIX divisé avec la matrice m
 */
inline Matrix Matrix::operator/(double x) const {
  Matrix LATRIX;

  for(int i=0; i<9; i++) {
    LATRIX[i] = tabMatrice[i] / x;
  }
  return LATRIX;
}


/*!
 * \brief Surcharge de l'opérateur =
 * \param m Matrix
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator=(const Matrix &m) {
    for(int i=0; i<9; i++) {
        tabMatrice[i] = m[i];
    }
    return *this;
}

/*!
 * \brief Surcharge de l'opérateur +=
 * \param m Matrix
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator+=(const Matrix &m) {
  for(int i=0; i<9; i++) {
    tabMatrice[i] = tabMatrice[i] + m[i];
  }
  return *this;
}


/*!
 * \brief Surcharge de l'opérateur -=
 * \param m Matrix
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator-=(const Matrix &m) {
  for(int i=0; i<9; i++) {
    tabMatrice[i] = tabMatrice[i] - m[i];
  }
  return *this;
}

/*!
 * \brief Surcharge de l'opérateur *=
 * \param m Matrix
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator*=(const Matrix &m) {
  Matrix LATRIX;

  LATRIX[0]  = tabMatrice[0] * m[0] + tabMatrice[1] * m[3] + tabMatrice[2] * m[6];
  LATRIX[1]  = tabMatrice[0] * m[1] + tabMatrice[1] * m[4] + tabMatrice[2] * m[7];
  LATRIX[2]  = tabMatrice[0] * m[2] + tabMatrice[1] * m[5] + tabMatrice[2] * m[8];

  LATRIX[3]  = tabMatrice[3] * m[0] + tabMatrice[4] * m[3] + tabMatrice[5] * m[6];
  LATRIX[4]  = tabMatrice[3] * m[1] + tabMatrice[4] * m[4] + tabMatrice[5] * m[7];
  LATRIX[5]  = tabMatrice[3] * m[2] + tabMatrice[4] * m[5] + tabMatrice[5] * m[8];

  LATRIX[6]  = tabMatrice[6] * m[0] + tabMatrice[7] * m[3] + tabMatrice[8] * m[6];
  LATRIX[7]  = tabMatrice[6] * m[1] + tabMatrice[7] * m[4] + tabMatrice[8] * m[7];
  LATRIX[8]  = tabMatrice[6] * m[2] + tabMatrice[7] * m[5] + tabMatrice[8] * m[8];

  *this = LATRIX;
  return *this;
}

/*!
 * \brief Surcharge de l'opérateur *=
 * \param x double
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator*=(double x) {
    for(int i=0; i<9; i++) {
        tabMatrice[i] = tabMatrice[i] * x;
    }
    return *this;
}

/*!
 * \brief Surcharge de l'opérateur /=
 * \param x double
 * \return Nouvelle matrice
 */
inline Matrix& Matrix::operator/=(double x) {
    for(int i=0; i<9; i++) {
        tabMatrice[i] = tabMatrice[i] / x;
    }
    return *this;
}

inline std::ostream& operator<<(std::ostream &os, const Matrix &m) {
    int v=0;

    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            os << m[v] << " ";
            v++;
        }
        os<<std::endl;
    }
    return os;
}

/*!
 * \brief Transpose une matrice
 * \return Matrix
 */
inline Matrix Matrix::transpose() {
    return Matrix(tabMatrice[0], tabMatrice[3], tabMatrice[6], tabMatrice[1], tabMatrice[4], tabMatrice[7], tabMatrice[2], tabMatrice[5], tabMatrice[8]);
}

/*!
 * \brief Inverse une matrice
 * \return Matrix
 */
inline Matrix Matrix::inverse() {
    double m11 = tabMatrice[4] * tabMatrice[8] - tabMatrice[5] * tabMatrice[7];
    double m12 = tabMatrice[3] * tabMatrice[8] - tabMatrice[5] * tabMatrice[6];
    double m13 = tabMatrice[3] * tabMatrice[7] - tabMatrice[4] * tabMatrice[6];
    double m21 = tabMatrice[1] * tabMatrice[8] - tabMatrice[2] * tabMatrice[7];
    double m22 = tabMatrice[0] * tabMatrice[8] - tabMatrice[2] * tabMatrice[6];
    double m23 = tabMatrice[0] * tabMatrice[7] - tabMatrice[1] * tabMatrice[6];
    double m31 = tabMatrice[1] * tabMatrice[5] - tabMatrice[2] * tabMatrice[4];
    double m32 = tabMatrice[0] * tabMatrice[5] - tabMatrice[2] * tabMatrice[3];
    double m33 = tabMatrice[0] * tabMatrice[4] - tabMatrice[1] * tabMatrice[3];

    Matrix com(m11, -m12, m13, -m21, m22, -m23, m31, -m32, m33);
    com = com.transpose();
    double det = tabMatrice[0] * m11 - tabMatrice[1] * m12 + tabMatrice[2] * m13;

    if(det != 0) {
        return com/det;
    }
    else{
        std::cout<<"L'inverse de cette matrice n'existe pas";
        return *this;
    }
}

/*!
 * \brief Applique une rotation en x sur une matrice
 * \return Matrix
 */
inline const Matrix Matrix::rot_x(const double &m) {
    return Matrix(1, 0, 0, 0, cos(Math::DegreeToRadian(m)), -sin(Math::DegreeToRadian(m)), 0, sin(Math::DegreeToRadian(m)), cos(Math::DegreeToRadian(m)));
}

/*!
 * \brief Applique une rotation en y sur une matrice
 * \return Matrix
 */
inline const Matrix Matrix::rot_y(const double &m) {
    return Matrix(cos(m), 0, sin(Math::DegreeToRadian(m)), 0, 1, 0, -sin(Math::DegreeToRadian(m)), 0, cos(Math::DegreeToRadian(m)));
}

/*!
 * \brief Applique une rotation en z sur une matrice
 * \param m double
 * \return Matrix
 */
inline const Matrix Matrix::rot_z(const double &m) {
    return Matrix(cos(Math::DegreeToRadian(m)), -sin(Math::DegreeToRadian(m)), 0, sin(Math::DegreeToRadian(m)), cos(Math::DegreeToRadian(m)), 0, 0, 0, 1);
}

/*!
 * \brief Applique une omotethie sur une matrice
 * \param x double
 * \param y double
 * \param z double
 * \return
 */
inline const Matrix Matrix::omth(const double &x, const double &y, const double &z) {
    return Matrix(x, 0, 0, 0, y, 0, 0, 0, z);
}

/*!
 * \brief Surcharge de l'operateur *
 * \param vec Vector
 * \return Vecteur multiplié avec une matrice
 */
inline Vector Matrix::operator* (const Vector &vec) const {
  return Vector(vec[0] * tabMatrice[0] + vec[1] * tabMatrice[1] + vec[2] * tabMatrice[2],
          vec[0] * tabMatrice[3] + vec[1] * tabMatrice[4] + vec[2] * tabMatrice[5],
          vec[0] * tabMatrice[6] + vec[1] * tabMatrice[7] + vec[2] * tabMatrice[8]);
}
