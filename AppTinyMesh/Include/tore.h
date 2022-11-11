#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"
#include <QtCore/qmath.h>

class Tore {
private:
    int rayon;
    Vector centre;
public:
    Tore(Vector centre, int rayon);
    ~Tore();

    int getRayon() const ;
    Vector getCentre () const;
};
