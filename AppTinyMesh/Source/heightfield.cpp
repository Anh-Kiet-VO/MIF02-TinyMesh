// HeightField

// Self include
#include "heightfield.h"

/*!
 * \brief Génère un terrain à partir d'une image, on peut aussi modifier l'échelle du terrain
 */
HeightField::HeightField() {
    centre = Vector(0, 0, 0);
    img = QImage("..\\TinyMesh-master\\Data\\terrain1.png","png");
    echelle = 10;
}

/*!
 * \brief Génère un terrain à partir d'une image et d'une échelle passé en paramètre
 * \param i QImage image du terrain
 * \param e int Echelle
 */
HeightField::HeightField(QImage i, int e) {
    img = i;
    centre = Vector(0, 0, 0);
    echelle = e;
}
