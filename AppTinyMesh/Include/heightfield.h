#pragma once

#include <vector>
#include <QImage>
#include "mathematics.h"

class HeightField {
private :
    Vector centre;
    QImage img;
    int echelle;

public :
    HeightField();
    HeightField(QImage, int);

    Vector getCentre() const {
        return centre;
    }

    QImage getImg() const {
        return img;
    }

    int getEchelle() const {
        return echelle;
    }
};
