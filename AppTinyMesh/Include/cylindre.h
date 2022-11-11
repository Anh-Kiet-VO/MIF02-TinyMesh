#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"
#include <QtCore/qmath.h>

class Cylindre {
private:
    float rayon;
    Vector centre;
    int hauteur;

public:
    Cylindre(Vector centre, int rayon, int hauteur);
    ~Cylindre();

    int getHauteur() const;
    float getRayon() const;
    Vector getCentre() const;
};
