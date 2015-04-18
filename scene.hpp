#ifndef SCENE_HPP
#define SCENE_HPP

#include <list>
#include "algebra.hpp"
#include "primitive.hpp"
#include "material.hpp"

class SceneNode {
public:
  SceneNode();
  SceneNode(const std::string& name);
  virtual ~SceneNode();
  
  // Hierarchy
  typedef std::list<SceneNode*> ChildList;

  const Matrix4x4& get_transform() const { return m_trans; }
  const Matrix4x4& get_translation() const { return m_translation; }
  const Matrix4x4& get_rotation() const { return m_rotation; }
  const Matrix4x4& get_scale() const { return m_scale; }

  const Matrix4x4& get_inverse() const { return m_invtrans; }
  const ChildList get_children() const { return m_children; }
  virtual Primitive* get_primitive() { return NULL; }
  virtual Material* get_material() { return NULL; }

  // Returns true if there are no child nodes
  virtual bool is_leaf();
  
  void set_transform(const Matrix4x4& m)
  {
    m_trans = m;
    // m_invtrans = m.invert();
  }

  virtual void set_scale(const Matrix4x4& m)
  {
    m_scale = m;
  }

  virtual void set_rotation(const Matrix4x4& m)
  {
    m_rotation = m;
  }

  virtual void set_translation(const Matrix4x4& m)
  {
    m_translation = m;
  }

  void set_transform(const Matrix4x4& m, const Matrix4x4& i)
  {
    m_trans = m;
    m_invtrans = i;
  }

  void add_child(SceneNode* child)
  {
    m_children.push_back(child);
  }

  void remove_child(SceneNode* child)
  {
    m_children.remove(child);
  }

  // Callbacks to be implemented.
  // These will be called from Lua.
  void rotate(char axis, double angle);
  void scale(const Vector3D& amount);
  void translate(const Vector3D& amount);

  // Returns true if and only if this node is a JointNode
  virtual bool is_joint() const;
  
protected:
  
  // Useful for picking
  int m_id;
  std::string m_name;

  // Transformations
  Matrix4x4 m_trans;
  Matrix4x4 m_invtrans;

  Matrix4x4 m_translation;
  Matrix4x4 m_rotation;
  Matrix4x4 m_scale;

  ChildList m_children;

};

class JointNode : public SceneNode {
public:
  JointNode(const std::string& name);
  virtual ~JointNode();

  virtual bool is_joint() const;

  void set_joint_x(double min, double init, double max);
  void set_joint_y(double min, double init, double max);

  struct JointRange {
    double min, init, max;
  };

  
protected:

  JointRange m_joint_x, m_joint_y;
};

class GeometryNode : public SceneNode {
public:
  GeometryNode();
  GeometryNode(const std::string& name,
               Primitive* primitive);
  virtual ~GeometryNode();

  Material* get_material() { return m_material; }
  // Material* get_material();

  Primitive* get_primitive() { return m_primitive; }
  // Material* get_primitive();

  void set_material(Material* material)
  {
    m_material = material;
  }

  void set_scale(const Matrix4x4& m)
  {
    m_scale = m;
    m_primitive->scale_sphere(m);
  }

  void set_translation(const Matrix4x4& m)
  {
    m_translation = m;
    m_primitive->translate_sphere(m);
  }

protected:
  Material* m_material;
  Primitive* m_primitive;
};

#endif
