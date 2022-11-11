# Problème rencontré
J'utilise l'ancienne version de TinyMesh.
Donc je n'ai pas les boutons d'utilisable Pour voir les formes il faut décommenter à la main pour avoir les formes souhaités.
C'est pour cela qu'il y a 3 BoxMeshExample.
Le 1er avec toutes les formes, le 2e avec l'objet complexe et le 3e avec le terrain.

Lors de la compilation il ce peut qu'il y ait des erreurs dans le fichier moc_qte.cpp je ne sais pas pourquoi mais il faut supprimer ses lignes :
case 6: _t->cylindreExample(); break;
case 7: _t->CylindreMesh(); break;
case 8: _t->SphereMesh(); break;
case 9: _t->CapsuleMesh(); break;

 
