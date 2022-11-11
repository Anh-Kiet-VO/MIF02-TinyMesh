// Cylindre

// Self include
#include "cylindre.h"

/*!
 * \brief Créer un cylindre à partir d'un vecteur centre, d'un rayon et d'une hauteur donné en paramètre
 * \param c Vector centre
 * \param r int rayon
 * \param h int hauteur
 */
Cylindre::Cylindre(Vector c, int r, int h) {
    centre = c;
    rayon = r;
    hauteur = h;
}

/*!
 * \brief Destructeur d'un cylindre
 */
Cylindre::~Cylindre() {

}

/*!
 * \brief Getter de hauteur
 * \return hauteur
 */
int Cylindre::getHauteur() const {
    return hauteur;
}

/*!
 * \brief Getter de rayon
 * \return rayon
 */
float Cylindre::getRayon() const {
    return rayon;
}


/*!
 * \brief Getter de centre
 * \return centre
 */
Vector Cylindre::getCentre() const {
    return centre;
}
