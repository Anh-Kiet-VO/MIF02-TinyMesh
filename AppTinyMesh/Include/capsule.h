#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"
#include <QtCore/qmath.h>

class Capsule {
private:
    int rayon;
    Vector centre;
    int hauteur;
public:
    Capsule(Vector centre, int rayon, int hauteur);
    ~Capsule();

    int getRayon() const ;
    int getHauteur() const ;
    Vector getCentre () const;
};
