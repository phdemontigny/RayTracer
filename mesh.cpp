//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include "mesh.hpp"
#include <iostream>

Mesh::Mesh(const std::vector<Point3D>& verts,
           const std::vector< std::vector<int> >& faces)
  : m_verts(verts),
    m_faces(faces)
{
  set_bounding_sphere();
}

//---------------------------------------------------------------------------------

void Mesh::set_bounding_sphere() {
  
  min_coord = Point3D(INFINITY,INFINITY,INFINITY);
  max_coord = Point3D(-INFINITY,-INFINITY,-INFINITY);

  for (Point3D vertex : m_verts)
  {
    for (int i=0; i<3; ++i) 
    {
      if (vertex[i] < min_coord[i]) {
        min_coord[i] = vertex[i];
      }
      else if (vertex[i] > max_coord[i]) {
        max_coord[i] = vertex[i];
      }
    }
  }

  double x_radius = (max_coord[0] - min_coord[0])/2.0;
  double y_radius = (max_coord[1] - min_coord[1])/2.0;
  double z_radius = (max_coord[2] - min_coord[2])/2.0;

  m_pos = getMidpoint(min_coord, max_coord);
  m_radius = distance3D(m_pos, max_coord);
}

//---------------------------------------------------------------------------------

bool Mesh::intersects(  Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                        double & tnear, Vector3D & normal ) 
{

  // Test intersection with bounding sphere
  if (bounding) {
    return intersects_sphere( ray_origin, ray_dir, trans, tnear, normal );
  }

  if ( fast && !intersects_sphere( ray_origin, ray_dir, trans, tnear, normal ) ) {
    return false;
  }

  double epsilon = .001;
  double tmin = INFINITY;
  // Optimize later by including a bounding box
  for( std::vector<int> face : m_faces) 
  {
    // get vertices
    Point3D A = trans*m_verts[ face[0] ];
    Point3D B = trans*m_verts[ face[1] ];
    Point3D C = trans*m_verts[ face[2] ];
    double t = INFINITY;
    // std::cerr << A << ", " << B << ", " << C << std::endl;
    if ( intersectTriangle(A,B,C,ray_origin,ray_dir,t) ) {
      if (t > epsilon && t < tmin) {
        tmin = t;
        Vector3D vector1 = A - B;
        Vector3D vector2 = A - C; 
        normal = vector1.cross(vector2);
      }
    }
  }
  if (tmin == INFINITY) {
    return false;
  }
  tnear = tmin;
  return true;
}


//---------------------------------------------------------------------------------


std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  std::cerr << "mesh({";
  for (std::vector<Point3D>::const_iterator I = mesh.m_verts.begin(); I != mesh.m_verts.end(); ++I) {
    if (I != mesh.m_verts.begin()) std::cerr << ",\n      ";
    std::cerr << *I;
  }
  std::cerr << "},\n\n     {";
  
  for (std::vector<Mesh::Face>::const_iterator I = mesh.m_faces.begin(); I != mesh.m_faces.end(); ++I) {
    if (I != mesh.m_faces.begin()) std::cerr << ",\n      ";
    std::cerr << "[";
    for (Mesh::Face::const_iterator J = I->begin(); J != I->end(); ++J) {
      if (J != I->begin()) std::cerr << ", ";
      std::cerr << *J;
    }
    std::cerr << "]";
  }
  std::cerr << "});" << std::endl;
  return out;
}
