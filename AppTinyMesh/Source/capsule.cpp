// Capsule

// Self include
#include "capsule.h"

/*!
 * \brief Créer une capsule à partir d'un vecteur centre, d'un rayon et d'une hauteur donné en paramètre
 * \param c Vector centre
 * \param r int rayon
 * \param h int hateur
 */
Capsule::Capsule(Vector c, int r, int h) {
    rayon = r;
    centre = c;
    hauteur = h;
}

/*!
 * \brief Destructeur d'une capsule
 */
Capsule::~Capsule() {

}

/*!
 * \brief Getter de rayon
 * \return rayon
 */
int Capsule::getRayon() const {
    return rayon;
}

/*!
 * \brief Getter de centre
 * \return centre
 */
Vector Capsule::getCentre() const {
    return centre;
}

/*!
 * \brief Getter de hauteur
 * \return hauteur
 */
int Capsule::getHauteur() const {
    return hauteur;
}
