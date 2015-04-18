//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include "primitive.hpp"

Primitive::Primitive()
{
}

//---------------------------------------------------------------------------------

Vector3D Primitive::get_normal( Point3D point ) {
  Vector3D n = point - m_pos;
  n.normalize();
  return n;
}

//---------------------------------------------------------------------------------

Primitive::~Primitive()
{
}

//---------------------------------------------------------------------------------

void Primitive::set_bounding_sphere() {
  
  min_coord = Point3D(INFINITY,INFINITY,INFINITY);
  max_coord = Point3D(-INFINITY,-INFINITY,-INFINITY);

  for (Point3D vertex : m_data)
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

  // std::cerr << m_radius << ", " << min_coord << ", " << max_coord << std::endl;
}

//---------------------------------------------------------------------------------

void Primitive::scale_sphere(Matrix4x4 trans) {
  max_coord = trans * max_coord;
  min_coord = trans * min_coord;
  m_radius = distance3D(m_pos, max_coord);
}

//---------------------------------------------------------------------------------

void Primitive::translate_sphere(Matrix4x4 trans) {
  // std::cerr << trans << std::endl;
  m_pos = trans * m_pos;
}

//---------------------------------------------------------------------------------

bool Primitive::intersects_sphere(  Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                                    double & tnear, Vector3D & normal  ) 
{

  Point3D center = m_pos;
  double radius = m_radius;
  double epsilon = .001;

  // std::cerr << m_radius << std::endl;

  double roots[2] = {INFINITY,INFINITY};

  // find roots of at^2 + bt + c
  Vector3D oc = ray_origin-center;
  double a = ray_dir.dot(ray_dir);    // r^2
  double b = 2*oc.dot(ray_dir);     // 2(origin-center)d
  double c = oc.dot(oc) - pow(radius,2);  // (origin-center)^2 - radius^2 

  int num_roots = quadraticRoots(a,b,c,roots);

  if ( num_roots > 0 ) {
    if (roots[0] < epsilon && roots[1] < epsilon ) { 
      return false;
    }
    else {
      if (roots[0] < epsilon) { tnear = roots[1]; }
      else if (roots[1] < epsilon) { tnear = roots[0]; }
      else { tnear = std::min(roots[0],roots[1]); }
      Point3D hit_point = ray_origin + (tnear*ray_dir);
      normal = get_normal( hit_point );
      return true;
    }
  }
  return false;
}

//---------------------------------------------------------------------------------

double Primitive::distance3D( Point3D point1, Point3D point2 )
{
  Vector3D segment = point1 - point2;
  double dist = segment.dot(segment);
  return sqrt(dist);
}

//---------------------------------------------------------------------------------


double Primitive::determinant(Matrix3x3 M)
{
  double a1 = M[0][0] * M[1][1] * M[2][2];
  double a2 = M[1][0] * M[2][1] * M[0][2];
  double a3 = M[2][0] * M[0][1] * M[1][2];

  double s1 = M[2][0] * M[1][1] * M[0][2];
  double s2 = M[2][1] * M[1][2] * M[0][0];
  double s3 = M[2][2] * M[1][0] * M[0][1];

  return a1 + a2 + a3 - s1 - s2 - s3;
}

//---------------------------------------------------------------------------------

Point3D Primitive::getMidpoint(Point3D p1, Point3D p2)
{
    double x = (p1[0] + p2[0]) / 2;
    double y = (p1[1] + p2[1]) / 2;
    double z = (p1[2] + p2[2]) / 2;

    return Point3D(x,y,z);
}

//---------------------------------------------------------------------------------


// pre: Assume A,B,C in correct right-hand order
// updates t to be the point of intersection
bool Primitive::intersectTriangle( Point3D A, Point3D B, Point3D C,
                              Point3D ray_origin, Vector3D ray_dir, double & tnear )
{

  // Source: http://people.cs.umass.edu/~elm/Teaching/591B_S06/assign_5_help.cpp

  double a, b, c, d, e, f, g, h, i, j, k, l;
  a = A[0] - B[0];
  b = A[1] - B[1];
  c = A[2] - B[2];

  d = A[0] - C[0];
  e = A[1] - C[1];
  f = A[2] - C[2];

  g = ray_dir[0];
  h = ray_dir[1];
  i = ray_dir[2];

  j = A[0] - ray_origin[0];
  k = A[1] - ray_origin[1];
  l = A[2] - ray_origin[2];

  // now apply Cramer's rule to solve for t, beta, gamma
  // We can short-circuit, meaning that if the bary coords are outside the triangle,
  // or t isn't in the given interval, we can return false without evaluating the rest
  
  double M      =   a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);
  
  tnear         = -(f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c)) / M;
  
  double gamma  =  (i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c)) / M;
  if (gamma < 0 || gamma > 1)
    return false;

  double beta   =  (j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g)) / M;
  if (beta < 0 || beta > (1 - gamma))
    return false;

  return true;
  /*
  My bad code

  Matrix3x3 R = Matrix3x3(  Vector3D( B[0]-A[0], C[0]-A[0], ray_dir[0] ),
                            Vector3D( B[1]-A[1], C[1]-A[1], ray_dir[1] ),
                            Vector3D( B[2]-A[2], C[2]-A[2], ray_dir[2] ) );

  Matrix3x3 R1 = Matrix3x3( Vector3D( ray_origin[0]-A[0], C[0]-A[0], ray_dir[0] ),
                            Vector3D( ray_origin[1]-A[1], C[1]-A[1], ray_dir[1] ),
                            Vector3D( ray_origin[2]-A[2], C[2]-A[2], ray_dir[2] ) );

  Matrix3x3 R2 = Matrix3x3( Vector3D( B[0]-A[0], ray_origin[0]-A[0], ray_dir[0] ),
                            Vector3D( B[1]-A[1], ray_origin[1]-A[1], ray_dir[1] ),
                            Vector3D( B[2]-A[2], ray_origin[2]-A[2], ray_dir[2] ) );

  Matrix3x3 R3 = Matrix3x3( Vector3D( B[0]-A[0], C[0]-A[0], ray_origin[0]-A[0] ),
                            Vector3D( B[1]-A[1], C[1]-A[1], ray_origin[1]-A[1] ),
                            Vector3D( B[2]-A[2], C[2]-A[2], ray_origin[2]-A[2] ) );

  double D = determinant( R );
  double epsilon = 1e-7;

  if (abs(D) == 0) {
    return false;
  }

  double D1 = determinant( R1 );
  double D2 = determinant( R2 );
  double D3 = determinant( R3 );

  double beta = D1/D;
  double gamma = D2/D;

  std::cerr << beta << ", " << gamma << std::endl;

  if (beta >= 0 && gamma >= 0 && (beta + gamma) <= 1) {
    tnear = D3/D;
    return true;
  }
  return false;

  */
}

//---------------------------------------------------------------------------------


bool Primitive::intersects( Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                            double & tnear, Vector3D & normal ) {

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
  for( unsigned i=0; i<m_data.size(); i+=3) 
  {
    // get vertices
    Point3D A = trans*m_data[i];
    Point3D B = trans*m_data[i+1];
    Point3D C = trans*m_data[i+2];
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

Sphere::Sphere()
{

  // Create a Mesh for the Sphere
	int max_depth = 4;
	// 4^(md-1) sub-triangles * 3 vertices per triangle
	int vert_per_triangle = pow(4,max_depth-1)*3;
	// 8 triangles in the diamond
	int total = vert_per_triangle*8;
  m_data.resize(total);

	// std::vector<Point3D> m_data (total);
  double radius = 1;

    // vertices of the original diamond
    Point3D v[] = {
        Point3D(0, 	0, 		1*radius),	 // A
        Point3D(0.5*radius,	0.5*radius,	0),	// B
        Point3D(0.5*radius,	-0.5*radius,	0),	// C
        Point3D(-0.5*radius,	-0.5*radius,	0),	// D
        Point3D(-0.5*radius,	0.5*radius,	0),	// E
        Point3D(0,		0,		-1*radius)	// F
    };

    Point3D diamondData[] = {
        v[0],v[2],v[1], // ACB
        v[0],v[3],v[2], // ADC
        v[0],v[4],v[3], // AED
        v[0],v[1],v[4], // ABE
        v[5],v[1],v[2], // FBC
        v[5],v[2],v[3], // FCD
        v[5],v[3],v[4], // FDE
        v[5],v[4],v[1], // FEB
    };

    for (int i=0; i<8; i++) {

        Point3D A = diamondData[3*i];
        Point3D B = diamondData[3*i+1];
        Point3D C = diamondData[3*i+2];

        triangleRecursion(A, B, C, 1, max_depth, 
            				      i*vert_per_triangle, (i+1)*vert_per_triangle, radius );
    }

    set_bounding_sphere();

}

//---------------------------------------------------------------------------------


Point3D Sphere::normalize(Point3D point, Point3D origin, double length) {

    double dx = point[0] - origin[0];
    double dy = point[1] - origin[1];
    double dz = point[2] - origin[2];

    double distance = pow(point[0],2) + pow(point[1],2) + pow(point[2],2);

    distance = sqrt(distance);

    dx = dx * (length / distance);
    dy = dy * (length / distance);
    dz = dz * (length / distance);

    return Point3D(origin[0] + dx, origin[1] + dy, origin[2] + dz);

}

//---------------------------------------------------------------------------------


void Sphere::triangleRecursion( Point3D A, Point3D B, Point3D C, 
                                int current_depth, int max_depth, 
                                int start, int end, double radius) {

    // std::cerr << "start: " << start << ", end: " << end << std::endl;
    Point3D origin = Point3D(0,0,0);
    A = normalize(A, origin, radius);
    B = normalize(B, origin, radius);
    C = normalize(C, origin, radius);
    if (current_depth < max_depth) {
        int i = (end-start)/4;
        Point3D D = getMidpoint(A,B);
        Point3D E = getMidpoint(B,C);
        Point3D F = getMidpoint(C,A);
        triangleRecursion(A,D,F,current_depth+1,max_depth,start,start+i,radius);
        triangleRecursion(D,B,E,current_depth+1,max_depth,start+i,start+2*i,radius);
        triangleRecursion(D,E,F,current_depth+1,max_depth,start+2*i,start+3*i,radius);
        triangleRecursion(F,E,C,current_depth+1,max_depth,start+3*i,end,radius);
    }
    else {
        m_data[start] = Point3D(A[0], A[1], A[2]);
        m_data[start+1] = Point3D(B[0], B[1], B[2]);
        m_data[start+2] = Point3D(C[0], C[1], C[2]);

        // code for finding normals
        // Vector3D vector1 = A - B;
        // Vector3D vector2 = A - C; 
        // Vector3D normal = Vector3D::normal(vector1,vector2);
    }
}

//---------------------------------------------------------------------------------

Sphere::~Sphere()
{
}

//---------------------------------------------------------------------------------

Cube::Cube()
{

  m_data.resize(36);

  // ORGANIZATION OF VERTICES:
  // Bottom: 0   1   Top: 4   5  -> y
  //         3   2        7   6
  //         FRONT
  // Origin: 5, Extends into positive XYZ region
  // NOTE: "cd" = "cubeData"
  Point3D cd[8] = {
            // X    Y    Z
      Point3D(-0.5, -0.5, -0.5),
      Point3D( 0.5, -0.5, -0.5),
      Point3D( 0.5, -0.5,  0.5),
      Point3D(-0.5, -0.5,  0.5),
      Point3D(-0.5,  0.5, -0.5),
      Point3D( 0.5,  0.5, -0.5),
      Point3D( 0.5,  0.5,  0.5),
      Point3D(-0.5,  0.5,  0.5)
  };

  // ORGANIZATION OF TRIANGLES
  // BOT - FRONT - RIGHT - BACK - LEFT - TOP
  // TOPLEFT - TOPRIGHT
  // Use the right hand rule for defining normal vectors
  m_data = {
      cd[0], cd[1], cd[3],
      cd[2], cd[3], cd[1],
      cd[7], cd[3], cd[6],
      cd[2], cd[6], cd[3],
      cd[6], cd[2], cd[5],
      cd[1], cd[5], cd[2],
      cd[5], cd[1], cd[4],
      cd[0], cd[4], cd[1],
      cd[4], cd[0], cd[7],
      cd[3], cd[7], cd[0],
      cd[4], cd[7], cd[5],
      cd[6], cd[5], cd[7]
  };

  set_bounding_sphere();
}

Cube::~Cube()
{
}

//---------------------------------------------------------------------------------

void NonhierSphere::scale_sphere(Matrix4x4 trans) {
}

//---------------------------------------------------------------------------------

void NonhierSphere::translate_sphere(Matrix4x4 trans) {
}

//---------------------------------------------------------------------------------

bool NonhierSphere::intersects( Point3D ray_origin, Vector3D ray_dir, Matrix4x4 trans,
                                double & tnear, Vector3D & normal  ) 
{

	return intersects_sphere( ray_origin, ray_dir, trans, tnear, normal );

}

//---------------------------------------------------------------------------------

NonhierSphere::~NonhierSphere()
{
}

NonhierBox::~NonhierBox()
{
}
