//---------------------------------------------------------------------------
//
// Name: Philippe Demontigny
//
// Student Number: 20557658
// User-id: pdemonti
// Assignment: 4
//
//---------------------------------------------------------------------------


#include "material.hpp"


Material::Material(const Colour& kd, const Colour& ks, double shininess)
  : m_kd(kd), m_ks(ks), m_shininess(shininess)
{
}

Material::~Material()
{
}

void Material::apply_gl() const
{
  // Perform OpenGL calls necessary to set up this material.
}
