#include "mesh.h"
#include "cylindre.h"
#include "sphere.h"
#include "tore.h"

#include "mathematics.h"

#include <QVector3D>
/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize 
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}

/*!
 * \brief Créer un cylindre
 * \param cylindre Cylindre
 */
Mesh::Mesh(const Cylindre& cylindre) {
    int i;
    const int div = (25 * 2);
    float alpha;
    float step = 2.0 * M_PI / (div -2);

    vertices.resize(div+2);
    vertices[div+1] = cylindre.getCentre() + Vector(0, cylindre.getHauteur(), 0);
    vertices[div] = cylindre.getCentre() + Vector(0, -(cylindre.getHauteur()), 0);

    for(i=0; i<div-1; i+=2) {
        alpha = i*step;

        normals.push_back(Vector(cos(alpha), 0, sin(alpha)));
        vertices[i] = cylindre.getCentre() + Vector(cos(alpha)*cylindre.getRayon(), cylindre.getHauteur(), sin(alpha)*cylindre.getRayon());

        normals.push_back(Vector(cos(alpha), 0, sin(alpha)));
        vertices[i+1] = cylindre.getCentre() + Vector(cos(alpha)*cylindre.getRayon(), -(cylindre.getHauteur()), sin(alpha)*cylindre.getRayon());
    }

    normals.push_back(Vector(0,1,0));
    normals.push_back(Vector(0,-1,0));

    varray.reserve((div+2) * 3);
    narray.reserve((div+2) * 3);

    for(i=0; i<(div-2); ++i) {
        AddSmoothTriangle(i, i, (i+1)%(div-2), (i+1)%(div-2), (i+2)%(div-2), (i+2)%(div-2));
    }

    for(i=0; i<div; i+=2) {
        AddSmoothTriangle(div, div, i+1, i+1, (i+3)%(div), (i+3)%(div));
        AddSmoothTriangle(div+1, div+1, i, i, (i+2)%(div), (i+2)%(div));
    }
}

/*!
 * \brief Créer une sphere
 * \param Sphere sphere
 */
Mesh::Mesh(const Sphere& sphere) {
    int divBeta = 50;
    int divAlpha = divBeta/2;
    double beta, alpha, alpha2;
    int sm = 0;

    vertices.resize(divBeta * divAlpha);

    for (int i=0; i < divAlpha; i++) {
        alpha = -0.5f * M_PI + float(i) * M_PI / divAlpha;
        alpha2 = -0.5f * M_PI + float(i+1) * M_PI / divAlpha;
        for (int j = 0; j < divBeta; j+=2) {
            beta = float(j) * 2.f * M_PI / (divBeta);

            normals.push_back(Vector(cos(alpha) * cos(beta), sin(alpha), cos(alpha) * sin(beta)));
            vertices[sm] = sphere.getCentre() + sphere.getRayon() * Vector(cos(alpha) * cos(beta), sin(alpha), cos(alpha) * sin(beta));

            normals.push_back(Vector(cos(alpha2) * cos(beta), sin(alpha2), cos(alpha2) * sin(beta)));
            vertices[sm+1] = sphere.getCentre() + sphere.getRayon() * Vector(cos(alpha2) * cos(beta), sin(alpha2), cos(alpha2) * sin(beta));

            sm += 2;
        }
    }

        varray.reserve(vertices.size() * 3);
        narray.reserve(vertices.size() * 3);

        for (int i = 0; i <= vertices.size()-3; i++) {
            AddSmoothTriangle(i, i, (i+1) % vertices.size(), (i+1) % vertices.size(), (i+2) % vertices.size(), (i+2) % vertices.size());
        }
}

/*!
 * \brief Créer une capsule
 * \param capsule Capsule
 */
Mesh::Mesh(const Capsule& capsule) {
    int i;
    const int div = (50 * 2);
    float alpha;
    float step = 2.0 * M_PI / (div -2);

    vertices.resize(div);

    for(i=0; i<=div-4; i+=2) {
        alpha = i*step;

        normals.push_back(Vector(cos(alpha), 0, sin(alpha)));
        vertices[i] = capsule.getCentre() + Vector(cos(alpha)*capsule.getRayon(), capsule.getHauteur(), sin(alpha)*capsule.getRayon());

        normals.push_back(Vector(cos(alpha), 0, sin(alpha)));
        vertices[i+1] = capsule.getCentre() + Vector(cos(alpha)*capsule.getRayon(), -(capsule.getHauteur()), sin(alpha)*capsule.getRayon());
    }

    varray.reserve((div) * 3);
    narray.reserve((div) * 3);

    for(i=0; i<(div-2); ++i) {
        AddSmoothTriangle(i, i, (i+1)%(div-2), (i+1)%(div-2), (i+2)%(div-2), (i+2)%(div-2));
    }

    Mesh sphereTop (Sphere(capsule.getCentre() + Vector(0,capsule.getHauteur(),0),capsule.getRayon()));
    merge(sphereTop);

    Mesh sphereBottom (Sphere(capsule.getCentre() + Vector(0,-(capsule.getHauteur()),0),capsule.getRayon()));
    merge(sphereBottom);
}

/*!
 * \brief Créer un tore
 * \param tore Tore
 */
Mesh::Mesh(const Tore& tore) {
    int st = 15;
    int sl = 15;
    int innerR = 1;
    int outerR = 5;

    float phi = 0.0;
    float dp = (2 * M_PI) / sl;

    float theta = 0.0;
    float dt = (2 * M_PI) / st;

    for (int i = 0; i < st; i++) {
        theta = dt * i;
        for (int j = 0; j < sl; j++) {
            phi = dp * j;

            Vector v(cos(theta) * (outerR + cos(phi) * innerR), sin(theta) * (outerR + cos(phi) * innerR), sin(phi) * innerR);
            vertices.push_back(tore.getCentre() + v );
            normals.push_back(v - tore.getCentre());
        }
    }

    varray.reserve(vertices.size() * 3);
    narray.reserve(vertices.size() * 3);

    for (int i = 0; i < vertices.size(); ++i) {
        AddSmoothTriangle(i, i, (i + 1) % vertices.size(), (i + 1) % vertices.size(), (i + sl) % vertices.size(), (i + sl) % vertices.size());
        AddSmoothTriangle((i + sl) % vertices.size(), (i + sl) % vertices.size(), (i + 1 + sl) % vertices.size(), (i + 1 + sl) % vertices.size(), (i + 1) % vertices.size(), (i + 1) % vertices.size());
    }
}

/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}



#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}

/*!
 * \brief Merge 2 mesh
 * \param m Mesh
 */
void Mesh::merge(const Mesh &m) {
    for(int k=0; k<m.Triangles(); k++) {
        this->AddTriangle(m.VertexIndex(k, 0) + Vertexes(), m.VertexIndex(k, 1) + Vertexes(), m.VertexIndex(k, 2) + Vertexes(),  m.NormalIndex(k, 0) + Normals());
    }

    for(int i=0; i<m.vertices.size(); i++) {
        this->vertices.push_back(m.vertices[i]);
    }

    for(int j=0; j<m.normals.size(); j++) {
        this->normals.push_back(m.normals[j]);
    }
}

/*!
 * \brief Déforme une sphère dans la direction d'un vector donné en paramètre
 * \param s Sphere
 * \param t Vector
 */
void Mesh::SphereWarp(const Sphere &s, const Vector &t) {
    Vector v;
    for(int i=0;i<vertices.size();i++)
    {
        double d= SquaredNorm(vertices[i]-s.getCentre());

        if(d<(s.getRayon()*s.getRayon()))
        {
            v = t * (1-(d/(s.getRayon()*s.getRayon()))) ;
            vertices[i]= vertices[i] + v;
        }
    }
}

/*!
 * \brief Modifie un mesh avec une matrice passé en paramètre
 * \param m Matrix
 */
void Mesh::edit(const Matrix& m) {
    for(int i=0; i<vertices.size(); i++) {
        vertices[i] = (m) * vertices[i];
    }

    for(int i=0; i<normals.size(); i++) {
        normals[i]= m * normals[i];
    }
}

/*!
 * \brief Créer un terrain
 * \param heightField HeightField
 */
Mesh::Mesh(const HeightField& heightField) {
    vertices.resize(heightField.getImg().width() * heightField.getImg().width() + heightField.getImg().height());
    varray.reserve(heightField.getImg().width() * heightField.getImg().width() + heightField.getImg().height());
    narray.reserve(heightField.getImg().width() * heightField.getImg().width() + heightField.getImg().height());

    for(int i=0; i<heightField.getImg().width()-1; i++) {
        for(int j=0; j<heightField.getImg().height()-1; j++) {
            int altitude = heightField.getImg().pixelColor(i,j).black();

            vertices[i*heightField.getImg().width()+j] = Vector(j, i, -altitude / heightField.getEchelle()) + heightField.getCentre();
            normals.push_back(Vector(j, i, 1));

            if((i != heightField.getImg().width()-2) ) {
                AddTriangle(i * heightField.getImg().width() + j, (i+1) * heightField.getImg().width() + j + 1,
                            i * heightField.getImg().width() + j + 1, i * (heightField.getImg().width() - 1) + j);
            }

            if( j != heightField.getImg().height()-2) {
                AddTriangle(i * heightField.getImg().width() + j, (i + 1) * heightField.getImg().width() + j + 1,
                            (i + 1) * heightField.getImg().width() + j, i * (heightField.getImg().width() - 1) + j);
            }
        }
    }
}

