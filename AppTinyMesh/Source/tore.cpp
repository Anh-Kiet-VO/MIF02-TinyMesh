// Tore

// Self include
#include "tore.h"

/*!
 * \brief Créer un tore à partir d'un vecteur centre et d'un rayon donné en paramètre
 * \param c
 * \param r
 */
Tore::Tore(Vector c, int r) {
    rayon = r;
    centre = c;
}

/*!
 * \brief Destructeur d'une sphere
 */
Tore::~Tore() {

}

/*!
 * \brief Getter de rayon
 * \return rayon
 */
int Tore::getRayon() const {
    return rayon;
}

/*!
 * \brief Getter de centre
 * \return centre
 */
Vector Tore::getCentre() const {
    return centre;
}
