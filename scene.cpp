//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include "scene.hpp"
#include <iostream>

SceneNode::SceneNode()
{

}

SceneNode::SceneNode(const std::string& name)
  : m_name(name)
{
}

SceneNode::~SceneNode()
{
}

void SceneNode::rotate(char axis, double angle)
{

  double q = ((M_PI * angle)/180.);

  if (axis == 'x') {
    Matrix4x4 transformation_matrix = Matrix4x4(  Vector4D(1, 0,      0,      0),
                                        Vector4D(0, cos(q), sin(q), 0),
                                        Vector4D(0, -sin(q),cos(q), 0),
                                        Vector4D(0, 0,      0,      1));
    m_rotation = transformation_matrix * m_rotation;
  }
  else if (axis == 'y') {
    Matrix4x4 transformation_matrix = Matrix4x4(  Vector4D(cos(q), 0, -sin(q), 0),
                                        Vector4D(0,      1, 0,       0),
                                        Vector4D(sin(q), 0, cos(q),  0),
                                        Vector4D(0,      0, 0,       1));
    m_rotation = transformation_matrix * m_rotation;
  }
  else if (axis == 'z') {
    Matrix4x4 transformation_matrix = Matrix4x4(  Vector4D(cos(q),  sin(q), 0,  0),
                                        Vector4D(-sin(q), cos(q), 0,  0),
                                        Vector4D(0,       0,      1,  0),
                                        Vector4D(0,       0,      0,  1));
    m_rotation = transformation_matrix * m_rotation;
  }
  else {
    std::cerr << "ERROR: Invalid axis -- " << axis << std::endl;
  }
}

void SceneNode::scale(const Vector3D& amount)
{
  Matrix4x4 transformation_matrix = Matrix4x4(  Vector4D(amount[0],  0,  0,  0),
                                                Vector4D(0,  amount[1],  0,  0),
                                                Vector4D(0,  0,  amount[2],  0),
                                                Vector4D(0,  0,  0,  1));

  m_scale = transformation_matrix * m_scale;
}

void SceneNode::translate(const Vector3D& amount)
{
  Matrix4x4 transformation_matrix = Matrix4x4(  Vector4D(1,  0,  0,  amount[0]),
                                                Vector4D(0,  1,  0,  amount[1]),
                                                Vector4D(0,  0,  1,  amount[2]),
                                                Vector4D(0,  0,  0,  1));

  m_translation = transformation_matrix * m_translation;
}

bool SceneNode::is_joint() const
{
  return false;
}

bool SceneNode::is_leaf() 
{
  return (m_children.size() == 0);
}

JointNode::JointNode(const std::string& name)
  : SceneNode(name)
{
}

JointNode::~JointNode()
{
}

bool JointNode::is_joint() const
{
  return true;
}

void JointNode::set_joint_x(double min, double init, double max)
{
  m_joint_x.min = min;
  m_joint_x.init = init;
  m_joint_x.max = max;
}

void JointNode::set_joint_y(double min, double init, double max)
{
  m_joint_y.min = min;
  m_joint_y.init = init;
  m_joint_y.max = max;
}

GeometryNode::GeometryNode(const std::string& name, Primitive* primitive)
  : SceneNode(name),
    m_primitive(primitive)
{
}

GeometryNode::GeometryNode() : m_primitive(NULL), m_material(NULL)
{
}

GeometryNode::~GeometryNode()
{
}
 
