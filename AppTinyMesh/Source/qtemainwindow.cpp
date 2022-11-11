#include "qte.h"
#include "implicits.h"


#include "mathematics.h"

MainWindow::MainWindow()
{
    // Chargement de l'interface
    uiw.setupUi(this);

    // Chargement du GLWidget
    meshWidget = new MeshWidget;
    QGridLayout* GLlayout = new QGridLayout;
    GLlayout->addWidget(meshWidget, 0, 0);
    GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw.widget_GL->setLayout(GLlayout);

    // Creation des connect
    CreateActions();

    meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
    delete meshWidget;
}

void MainWindow::CreateActions()
{
    // Buttons
    connect(uiw.boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw.sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw.resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw.wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw.radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw.radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

    // Widget edition
    connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
    connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

/*
    /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\

  J'utilise l'ancienne version de TinyMesh. Donc je n'ai pas les boutons d'utilisable
  Pour voir les formes il faut décommenter à la main pour avoir les formes souhaités
  Dans le second BoxMeshExample il y a l'objet complexe

    /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
*/

void MainWindow::BoxMeshExample()
{
    // Cylindre :
    //Mesh MESH = Mesh(Cylindre(Vector(0,0,0),2, 2));

    // Capsule :
    //Mesh MESH = Mesh(Capsule(Vector(0,0,0),1,1));

    // Tore :
    //Mesh MESH = Mesh(Tore(Vector(0,0,0),1));

    // Sphere :
    Mesh MESH = Mesh(Sphere(Vector(0,0,0),1));

    // Avec le Sphere Warp d'appliqué :
    MESH.SphereWarp(Sphere(Vector(0,0,0),1),Vector(2,0,2));

    std::vector<Color> cols;
    cols.resize(MESH.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(MESH, cols, MESH.VertexIndexes());
    UpdateGeometry();
}

// Objet complexe
//void MainWindow::BoxMeshExample()
//{
//    Mesh Roue1 = Mesh(Tore(Vector(0,0,0),5));
//    Roue1.edit(Matrix::rot_x(90));
//    Mesh Roue2 = Mesh(Tore(Vector(30,0,0),5));
//    Roue2.edit(Matrix::rot_x(90));
//    Roue1.merge(Roue2);

//    Mesh Cylindre1 = Mesh(Cylindre(Vector(0,-22,0),1,8));
//    Cylindre1.edit((Matrix::rot_z(90)));
//    Roue1.merge(Cylindre1);

//    Mesh RondArriere = Mesh(Sphere(Vector(30,0,0),2));
//    Roue1.merge(RondArriere);

//    Mesh Cylindre2 = Mesh(Cylindre(Vector(0,10,0),1,10));
//    Cylindre2.edit((Matrix::rot_x(90)));
//    Roue1.merge(Cylindre2);

//    Mesh RoueDevant = Mesh(Sphere(Vector(0,0,0),2));
//    Roue1.merge(RoueDevant);

//    Mesh Cylindre3 = Mesh(Cylindre(Vector(0,3,11),1,11));
//    Cylindre3.edit((Matrix::rot_z(90)));
//    Cylindre3.edit((Matrix::rot_y(45)));
//    Roue1.merge(Cylindre3);

//    Mesh BasePedale = Mesh(Sphere(Vector(13,0,1),2));
//    BasePedale.edit((Matrix::omth(1,1,2)));
//    Roue1.merge(BasePedale);

//    Mesh Pedales = Mesh(Capsule(Vector(13,0,2),1,5));
//    Roue1.merge(Pedales);

//    Mesh Cylindre4 = Mesh(Cylindre(Vector(0,-45,20.7),1,12));
//    Cylindre4.edit((Matrix::rot_z(90)));
//    Cylindre4.edit((Matrix::rot_y(20)));
//    Cylindre4.edit((Matrix::omth(1,1,2)));
//    Roue1.merge(Cylindre4);

//    Mesh Cylindre5 = Mesh(Cylindre(Vector(0,-8,9),1,8));
//    Cylindre5.edit((Matrix::rot_z(90)));
//    Cylindre5.edit((Matrix::omth(1.2,1,2)));
//    Roue1.merge(Cylindre5);

//    Mesh Cylindre6 = Mesh(Cylindre(Vector(0,-15,-10),1,8));
//    Cylindre6.edit((Matrix::rot_x(-61)));
//    Cylindre6.edit((Matrix::rot_z(90)));
//    Roue1.merge(Cylindre6);

//    Mesh Siege = Mesh(Sphere(Vector(20,0,9),2));
//    Siege.SphereWarp(Sphere(Vector(22,0,9),1),Vector(5,0,0));
//    Siege.edit((Matrix::omth(1,2,2)));
//    Roue1.merge(Siege);

//    Mesh Volant = Mesh(Cylindre(Vector(-10,0,-5),1,10));
//    Volant.edit((Matrix::rot_y(90)));
//    Volant.edit((Matrix::omth(1,1,2)));
//    Roue1.merge(Volant);

//    std::vector<Color> cols;
//    cols.resize(Roue1.Vertexes());

//    for (int i = 0; i < cols.size(); i++) {
//        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);
//    }

//    meshColor = MeshColor(Roue1, cols, Roue1.VertexIndexes());

//    UpdateGeometry();
//}

// Terrain
//void MainWindow::BoxMeshExample() {
//    QImage img = QImage("..\\TinyMesh-master\\Data\\terrain1.png","png");
//    Mesh terrain = Mesh(HeightField(img, 10));

//    std::vector<Color> cols;
//    cols.resize(terrain.Vertexes());

//    for (int i = 0; i < cols.size(); i++) {
//        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);
//    }

//    meshColor = MeshColor(terrain, cols, terrain.VertexIndexes());
//    UpdateGeometry();
//}

void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (int i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();

}

void MainWindow::UpdateGeometry()
{
    meshWidget->ClearAll();
    meshWidget->AddMesh("BoxMesh", meshColor);

    uiw.lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw.lineEdit_2->setText(QString::number(meshColor.Triangles()));

    UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw.wireframe->isChecked());

    if (uiw.radioShadingButton_1->isChecked())
        meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
    else
        meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
    meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
