// Sphere

// Self include
#include "sphere.h"

/*!
 * \brief Créer une sphere à partir d'un vecteur centre et d'un rayon donné en paramètre
 * \param c Vector centre
 * \param r int rayon
 */
Sphere::Sphere(Vector c, int r) {
    rayon = r;
    centre = c;
}

/*!
 * \brief Destructeur d'une sphere
 */
Sphere::~Sphere() {

}

/*!
 * \brief Getter de rayon
 * \return rayon
 */
int Sphere::getRayon() const {
    return rayon;
}


/*!
 * \brief Getter de centre
 * \return centre
 */
Vector Sphere::getCentre() const {
    return centre;
}
