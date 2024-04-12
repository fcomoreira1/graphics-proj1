#ifndef VECTOR_H
#define VECTOR_H
class Vector {
  public:
    Vector();
    explicit Vector(double x, double y, double z);
    double norm2() const;
    double norm() const;
    void normalize();
    double operator[](int i) const;
    double &operator[](int i);
    void print();
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b);
Vector operator-(const Vector &a, const Vector &b);
Vector operator*(const double a, const Vector &b);
Vector operator*(const Vector &a, const double b);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector &a, const double b);
double dot(const Vector &a, const Vector &b);
Vector cross(const Vector &a, const Vector &b);
#endif
