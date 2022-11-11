#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"
#include <QtCore/qmath.h>

class Sphere {
private:
    int rayon;
    Vector centre;
public:
    Sphere(Vector centre, int rayon);
    ~Sphere();

    int getRayon() const ;
    Vector getCentre () const;
};
