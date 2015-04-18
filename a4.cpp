//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include "a4.hpp"
#include "image.hpp"
#include <cassert>

//---------------------------------------------------------------------------------


// Returns a list of all nodes to be rendered
// TODO: update hierarchical nodes
void extract_nodes( SceneNode* root, 
                    std::vector<SceneNode *> & node_list) 
{   
    for(SceneNode* node : root->get_children())
    {
        Matrix4x4 current_scale = root->get_scale() * node->get_scale();
        Matrix4x4 current_rotation = root->get_rotation() * node->get_rotation();
        Matrix4x4 current_translation = root->get_translation() * node->get_translation();
        node->set_scale( current_scale );
        node->set_rotation( current_rotation );
        node->set_translation( current_translation );
        node->set_transform( current_translation * current_rotation * current_scale );
        if ( node->get_primitive() != NULL ) {
            node_list.push_back(node);
        }
        if ( !node->is_leaf() ) {
            extract_nodes(node, node_list);
        }
    }
}

//---------------------------------------------------------------------------------

Vector3D ggReflection( Vector3D V, Vector3D N ) {
  double VdN = V.dot(N);
  return V - 2*VdN*N;
}

//---------------------------------------------------------------------------------

Colour directLight( Point3D origin,
                    Light * source,
                    std::vector<SceneNode *> node_list,
                    SceneNode * previous_node )
{
  Vector3D shadow_ray = source->position - origin;
  double t = INFINITY;

  for(SceneNode * node : node_list)
  {
    Primitive *prim = node->get_primitive();
    Matrix4x4 trans = node->get_transform();
    Vector3D normal = Vector3D(0,0,0);
    if ( prim->intersects(origin, shadow_ray, trans, t, normal) ) {
      // If in shadow, return nothing
      // std::cerr << origin << std::endl;
      return Colour(0,0,0);
    }
  }
  // If no objects in the way
  return source->colour;
}

//---------------------------------------------------------------------------------

Colour generateBackground(int y, int height)
{
  float proportion = 1 - float(y) / float(height);
  return Colour(0,0,proportion);
}

//---------------------------------------------------------------------------------


Colour trace(	Point3D eye,
              Point3D origin,
						  Vector3D ray,
						  std::vector<SceneNode *> node_list,
              const Colour& ambient,
              const std::list<Light*>& lights, 
              int depth, int y, int height)
{

  // t value for which the ray first intersects the sphere
  double tnear            = INFINITY;
  SceneNode *closest_node = NULL;
  Colour background       = Colour(0,0,0);
  Vector3D N;

  for(SceneNode * node : node_list)
  {
    Primitive *prim = node->get_primitive();
    Matrix4x4 trans = node->get_transform();
    double t = INFINITY;
    Vector3D normal = Vector3D(0,0,0);
    if ( prim->intersects(origin, ray, trans, t, normal) ) {
      if (t < tnear) { 
        N             = normal;
        tnear         = t;
        closest_node  = node;
      }
    }
  }

  // if no intersection, return background color
  if (closest_node == NULL) {
    return generateBackground(y, height);
  }
  // return Colour(0,0,1);

  // std::cerr << N << std::endl;

  Material* mat       = closest_node->get_material();
  Primitive *prim     = closest_node->get_primitive();
  Point3D hit_point   = origin + tnear*ray;

  Vector3D V          = eye - hit_point;
  V.normalize();
  N.normalize();

  Colour kd           = mat->get_kd();
  Colour ks           = mat->get_ks();
  double alpha        = mat->get_shine();
  Colour pixel_color  = kd*ambient;

  for (Light * source : lights) {
    Vector3D L        = source->position - hit_point;
    L.normalize();
    Vector3D R        = 2 * (L.dot(N)) * N - L;
    R.normalize();

    Colour light      = directLight(hit_point, source, node_list, closest_node);
    double LdN        = L.dot(N);
    Colour diffuse    = LdN * light;
    double RdV        = pow(R.dot(V),alpha);
    Colour specular   = RdV * light;
    pixel_color       = pixel_color + (kd * diffuse) + (ks * specular);

  }
  // Specular Lighting
  if ( depth < 5 ) {
    Vector3D reflected  = ggReflection( ray, N );
    reflected.normalize();

    pixel_color = pixel_color + 
                  ks * trace( eye,hit_point,reflected,node_list,
                              ambient,lights,depth+1,y,height);
  }
  return pixel_color;
}


//---------------------------------------------------------------------------------


void a4_render(// What to render
               SceneNode* root,
               // Where to output the image
               const std::string& filename,
               // Image size
               int width, int height,
               // Viewing parameters
               const Point3D& eye, const Vector3D& view,
               const Vector3D& up, double fov,
               // Lighting parameters
               const Colour& ambient,
               const std::list<Light*>& lights
               )
{

  // Fill in raytracing code here.
  // img(x,y,0) = amount of RGB Red
  // img(x,y,1) = amount of RGB Green
  // img(x,y,2) = amount of RGB Blue
  Image img(width, height, 3);

  std::vector<SceneNode *> node_list; 
  extract_nodes(root, node_list);

  float aspect = float(width) / float(height);
  float angle = tan((M_PI * 0.5 * fov)/180.); 

  // change notation to match notes
  float nx = float(width);
  float ny = float(height);

  // d = distance to near plane
  float d = abs(float(eye[2])) + 1;
  // w,h = width, height of the view plane
  float h = 2*d*tan(angle);

  // timer for rendering
  float timer = 0;
  float max_time = width*height;
  int ticker = 0;

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
    	// NOTE: Assume view frame = height x width = pixels
    	Point3D Pk = Point3D(x,y,d);

      // Step 1: Translate pixel point onto the near plane
      Matrix4x4 T1 = Matrix4x4( Vector4D(1,0,0,-nx/2),
      													Vector4D(0,1,0,-ny/2),
      													Vector4D(0,0,1,d),
      													Vector4D(0,0,0,1) ); 

      // Step 2: Scale to preserve aspect ratio and correct sign
      Matrix4x4 S2 = Matrix4x4( Vector4D(-h/ny,	0,		0,0),
      													Vector4D(0,			h/ny,	0,0),
      													Vector4D(0,			0,		1,0),
      													Vector4D(0,			0,		0,1) ); 

      // Step 3: Rotate to superimpose WCS onto VCS
      Vector3D w = view;
      w.normalize();
      Vector3D u = up.cross(w);
      u.normalize();
      Vector3D v = w.cross(u);

      Matrix4x4 R3 = Matrix4x4( Vector4D(u[0], 	v[0], w[0],	0),
      													Vector4D(u[1], 	v[1], w[1],	0),
      													Vector4D(u[2], 	v[2], w[2],	0),
      													Vector4D(0, 		0, 		0,		1) );

      // Step 4: Translate to World Coordinates
      Matrix4x4 T4 = Matrix4x4( Vector4D(1,0,0,eye[0]),
      													Vector4D(0,1,0,eye[1]),
      													Vector4D(0,0,1,eye[2]),
      													Vector4D(0,0,0,1) );

      float invWidth  = 1 / float(width);
      float invHeight = 1 / float(height); 
      float xx = (2 * (x * invWidth) - 1) * angle * aspect;
      float yy = (1 - 2 * (y * invHeight)) * angle; 


      Point3D Pw = T4 * R3 * S2 * T1 * Pk;
      Vector3D ray = Vector3D(xx,yy,-1);
      ray.normalize();

      Colour pixel = trace(eye, eye, ray, node_list, ambient, lights, 1, y, height);
      img(x, y, 0) = pixel.R();
      img(x, y, 1) = pixel.G();
      img(x, y, 2) = pixel.B();

      timer++;
      if (timer/max_time >= .1) {
        ticker += 10;
        std::cerr << "---" << ticker << "% COMPLETE---" << std::endl;
        timer = 0.0;
      }
    }
  }

  std::cerr << "---100% COMPLETE---" << std::endl;

  for (std::list<Light*>::const_iterator I = lights.begin(); I != lights.end(); ++I) {
    if (I != lights.begin()) std::cerr << ", ";
    std::cerr << **I;
  }
  std::cerr << "});" << std::endl;
  
  // For now, just make a sample image.
  /*
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < height; x++) {
      // Red: increasing from top to bottom
      img(x, y, 0) = (double)y / height;
      // Green: increasing from left to right
      img(x, y, 1) = (double)x / width;
      // Blue: in lower-left and upper-right corners
      img(x, y, 2) = ((y < height/2 && x < height/2)
                      || (y >= height/2 && x >= height/2)) ? 1.0 : 0.0;
    }
  }
  */
  img.savePng(filename);
  
}
