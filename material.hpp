#ifndef CS488_MATERIAL_HPP
#define CS488_MATERIAL_HPP

#include "algebra.hpp"

class Material {
public:
  Material(const Colour& kd, const Colour& ks, double shininess);
  virtual ~Material();

  Colour get_kd() { return m_kd; }
  Colour get_ks() { return m_ks; }
  double get_shine() { return m_shininess; }

  virtual void apply_gl() const;

private:
  Colour m_kd;
  Colour m_ks;

  double m_shininess;
};


#endif
