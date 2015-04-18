#ifndef CS488_PRIMITIVE_HPP
#define CS488_PRIMITIVE_HPP

#include "algebra.hpp"
#include "polyroots.hpp"
#include <cmath>
#include <vector>

extern bool bounding;
extern bool fast;

enum class PrimitiveType {
  Sphere,
  NonhierSphere,
  Cube,
  NonhierBox,
  Mesh,
  None
};

//---------------------------------------------------------------------------------

class Primitive {
public:
  Primitive();
  virtual ~Primitive();
  PrimitiveType get_type() { return m_type; }
  double determinant(Matrix3x3 M);
  bool intersectTriangle( Point3D A, Point3D B, Point3D C,
                          Point3D ray_origin, Vector3D ray_dir, double & tnear );

  virtual bool intersects(  Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                            double & tnear, Vector3D & normal );
  virtual Vector3D get_normal( Point3D point );
  virtual Point3D getMidpoint(Point3D p1, Point3D p2);

  virtual bool intersects_sphere(  Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                                  double & tnear, Vector3D & normal  );

  virtual void set_bounding_sphere();
  virtual void scale_sphere( Matrix4x4 trans );
  virtual void translate_sphere( Matrix4x4 trans );
  virtual double distance3D( Point3D point1, Point3D point2 );

protected:
  PrimitiveType m_type;
  std::vector<Point3D> m_data;
  Point3D m_pos;
  Point3D max_coord;
  Point3D min_coord;
  double m_radius;
};

//---------------------------------------------------------------------------------

class Sphere : public Primitive {
public:
  Sphere();
  virtual ~Sphere();

private:
  void triangleRecursion( Point3D A, Point3D B, Point3D C, 
                          int current_depth, int max_depth, 
                          int start, int end, double radius);
  Point3D normalize(Point3D point, Point3D origin, double length);

};

//---------------------------------------------------------------------------------

class Cube : public Primitive {
public:
  Cube();
  virtual ~Cube();  
};

//---------------------------------------------------------------------------------

class NonhierSphere : public Primitive {
public:
  NonhierSphere(const Point3D& pos, double radius)
  {
    m_pos = pos;
    m_radius = radius;
    max_coord = Point3D(radius, radius, radius);
    min_coord = Point3D(-radius, -radius, -radius);
  }
  virtual ~NonhierSphere();

  Point3D get_center() { return m_pos; }
  double get_radius() { return m_radius; }

  void scale_sphere( Matrix4x4 trans );
  void translate_sphere( Matrix4x4 trans );

  bool intersects(  Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                    double & tnear, Vector3D & normal );


};

//---------------------------------------------------------------------------------

class NonhierBox : public Primitive {
public:
  NonhierBox(const Point3D& pos, double size)
    : m_pos(pos), m_size(size)
  {
  }
  
  virtual ~NonhierBox();

private:
  Point3D m_pos;
  double m_size;
};

//---------------------------------------------------------------------------------


#endif
