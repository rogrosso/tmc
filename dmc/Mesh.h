#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>

// project
#include "Vector.h"

namespace p_dmc {
	class Mesh {
	public:
		typedef typename Vector Vertex;
		typedef typename Vector Normal;
		typedef typename std::array<int, 3> Triangle;
		typedef typename std::array<int, 4> Quadrilateral;
		typedef typename std::array<int, 2> Line;
		typedef typename std::array<double, 4> Attribute;
		typedef typename std::array<Vertex, 8> BBox;
	public:
		void resizeVertices(const int size_)
		{
			m_vertices.resize(size_);
			m_normals.resize(size_);
			m_attributes.resize(size_);
		}
		void addVertex(const int pos, Vertex& v, Normal& n)
		{
			m_vertices[pos] = v;
			m_normals[pos] = n;
			/*const double r = 250.0 / 255.0;
			const double g = 235.0 / 255.0;
			const double b = 215.0 / 255.0;
			const double a = 1.0;*/
            /*const double r = 175.0 / 255.0;
            const double g = 170.0 / 255.0;
            const double b = 170.0 / 255.0;
            const double a = 1.0;*/
            const double r = 250.0 / 255.0;
            const double g = 250.0 / 255.0;
            const double b = 250.0 / 255.0;
            const double a = 1.0;

			m_attributes[pos] = { r,g,b,a };
		}
		void addVertex(Vertex& v, Normal& n)
		{
			m_vertices.push_back(v);
			m_normals.push_back(n);
			//const double r = 250.0 / 255.0;
			//const double g = 235.0 / 255.0;
			//const double b = 215.0 / 255.0;
			//const double a = 1.0;
            /*const double r = 235.0 / 255.0;
            const double g = 235.0 / 255.0;
            const double b = 235.0 / 255.0;
            const double a = 1.0;*/
            /*const double r = 175.0 / 255.0;
            const double g = 170.0 / 255.0;
            const double b = 170.0 / 255.0;
            const double a = 1.0;*/
            const double r = 250.0 / 255.0;
            const double g = 250.0 / 255.0;
            const double b = 250.0 / 255.0;
            const double a = 1.0;
			m_attributes.push_back({ r,g,b,a });
		}
		void addVertex(Vertex& v, Normal& n, Attribute& a)
		{
			m_vertices.push_back(v);
			m_normals.push_back(n);
			m_attributes.push_back(a);
		}

		void resizeTriangles(const int size_)
		{
			m_triangles.resize(size_);
		}
		void addTriangle(const int pos, Triangle& t)
		{
			m_triangles[pos] = t;
		}
		void addTriangle(Triangle& t) { m_triangles.push_back(t); }

		void resizeQuadrilaterals(const int size_)
		{
			m_quadrilaterals.resize(size_);
		}

		void addQuadrilateral(const int pos, Quadrilateral& q)
		{
			m_quadrilaterals[pos] = q;
		}

		void addLineVertex(Vertex& v, Normal& n, Attribute& a)
		{
			m_lineVertices.push_back(v);
			m_lineNormals.push_back(n);
			m_lineAttributes.push_back(a);
		}
		void resizeLines(const int size_) { m_lines.resize(size_); }
		void addLine(Line& l) { m_lines.push_back(l); }
		void addLine(const int pos, const int v0, const int v1) { m_lines[pos] = { v0,v1 }; }
		void setBBox(BBox& bbox) {
			m_bbox = BBox{ bbox };
		}


		auto&& getTriangles() { return m_triangles; }
		auto&& getTriangles() const { return m_triangles; }
		auto&& getQuadrilaterals() { return m_quadrilaterals; }
		auto&& getQuadrilaterals() const { return m_quadrilaterals; }
		auto&& getVertices() { return m_vertices; }
		auto&& getVertices() const { return m_vertices; }
		auto&& getNormals() { return m_normals; }
		auto&& getNormals() const { return m_normals; }
		auto&& getAttributes() { return m_attributes; }
		auto&& getAttributes() const { return m_attributes; }

        auto trianglesData() { return m_triangles.data()->data(); }
        auto quadrilateralsData() { return m_quadrilaterals.data()->data(); }
        auto verticesData() { return m_vertices[0].data(); }
        auto normalsData() { return m_normals[0].data(); }
        auto linesData() { return m_lines.data()->data(); }

		auto&& getLineVertices() { return m_lineVertices; }
		auto&& getLineNormals() { return m_lineVertices; }
		auto&& getLineAttributes() { return m_lineVertices; }
		auto&& getLines() { return m_lines; }

		Vertex vertex(const int i) { return m_vertices[i]; }
		const Vertex vertex(const int i) const { return m_vertices[i]; }
		void vertex(const int i, const Vertex& v) { m_vertices[i] = v; }

		Normal normal(const int i) { return m_normals[i]; }
		const Normal normal(const int i) const { return m_normals[i]; }
		
		Triangle triangle(const int i) { return m_triangles[i]; }
		const Triangle triangle(const int i) const { return m_triangles[i]; }
		
		Quadrilateral quadrilateral(const int i) { return m_quadrilaterals[i]; }
		Quadrilateral quadrilateral(const int i) const { return m_quadrilaterals[i]; }
		
		int quadVertex(const int q, const int i) { return m_quadrilaterals[q][i]; }
		int quadVertex(const int q, const int i) const { return m_quadrilaterals[q][i]; }
		
		Attribute attribute(const int i) { return m_attributes[i]; }
		const Attribute attribute(const int i) const { return m_attributes[i]; }
		void attribute(const int i, Attribute& a) {
			m_attributes[i] = a;
		}

		Vertex lineVertex(const int i) { return m_lineVertices[i]; }
		Normal lineNormal(const int i) { return m_lineNormals[i]; }
		Attribute lineAttribute(const int i) { return m_lineAttributes[i]; }
		const Vertex lineVertex(const int i) const { return m_lineVertices[i]; }
		const Normal lineNormal(const int i) const { return m_lineNormals[i]; }
		const Attribute lineAttribute(const int i) const { return m_lineAttributes[i]; }
		size_t nrLineVertices() { return m_lineVertices.size(); }
		Line line(const int i) { return m_lines[i]; }
		const Line line(const int i) const { return m_lines[i]; }

		// modify elements
		void attribute(const int index, const Attribute& a) {
			m_attributes[index] = a;
		}

		size_t nrVertices() { return m_vertices.size(); }
		size_t nrVertices() const { return m_vertices.size(); }
		size_t nrTriangles() { return m_triangles.size(); }
		size_t nrTriangles() const { return m_triangles.size(); }
		size_t nrQuadrilaterals() { return m_quadrilaterals.size(); }
		size_t nrQuadrilaterals() const { return m_quadrilaterals.size(); }
		size_t nrLines() { return m_lines.size(); }
		size_t nrLines() const { return m_lines.size(); }

		void flipNormals()
		{
			for (auto& n : m_normals)
			{
				n.flip();
			}
		}

		double minX() { return m_bbox[0][0]; }
		double minY() { return m_bbox[0][1]; }
		double minZ() { return m_bbox[0][2]; }
		double maxX() { return m_bbox[7][0]; }
		double maxY() { return m_bbox[7][1]; }
		double maxZ() { return m_bbox[7][2]; }
		double minX() const { return m_bbox[0][0]; }
		double minY() const { return m_bbox[0][1]; }
		double minZ() const { return m_bbox[0][2]; }
		double maxX() const { return m_bbox[7][0]; }
		double maxY() const { return m_bbox[7][1]; }
		double maxZ() const { return m_bbox[7][2]; }
		void bbox()
		{
			double xMin = std::numeric_limits<double>::max();
			double xMax = -xMin;
			double yMin = xMin;
			double yMax = -xMin;
			double zMin = xMin;
			double zMax = -xMin;
			for (auto p : m_vertices)
			{
				xMin = (xMin > p[0]) ? p[0] : xMin;
				xMax = (xMax < p[0]) ? p[0] : xMax;
				yMin = (yMin > p[1]) ? p[1] : yMin;
				yMax = (yMax < p[1]) ? p[1] : yMax;
				zMin = (zMin > p[2]) ? p[2] : zMin;
				zMax = (zMax < p[2]) ? p[2] : zMax;
			}
			m_bbox[0] = Vertex{ xMin, yMin, zMin };
			m_bbox[1] = Vertex{ xMax, yMin, zMin };
			m_bbox[2] = Vertex{ xMin, yMax, zMin };
			m_bbox[3] = Vertex{ xMax, yMax, zMin };
			m_bbox[4] = Vertex{ xMin, yMin, zMax };
			m_bbox[5] = Vertex{ xMax, yMin, zMax };
			m_bbox[6] = Vertex{ xMin, yMax, zMax };
			m_bbox[7] = Vertex{ xMax, yMax, zMax };
		}
		void clear()
		{
			m_vertices.clear();
			m_normals.clear();
			m_triangles.clear();
			m_quadrilaterals.clear();
			m_attributes.clear();
			m_lines.clear();
			m_lineVertices.clear();
			m_lineNormals.clear();
			m_lineAttributes.clear();
		}
        void writeObjTriangles(std::string& o_file)
        {
            std::ofstream objF;
            objF.open(o_file.c_str());
            if (!objF.is_open()) {
                std::cout << "ERROR: can't open output file " << std::endl;
            }
            else {
                objF << "#Dual Marching Cubes\n";
                const int nr_v = this->nrVertices();
                for (int i = 0; i < nr_v; i++)
                {
                    Mesh::Vertex v = this->vertex(i);
                    Mesh::Attribute a = this->attribute(i);
                    objF << "v " << v[0] << " " << v[1] << " " << v[2] << " " << a[0] << " " << a[1] << " " << a[2] << std::endl;
                }
                for (auto n : this->getNormals())
                {
                    objF << "vn " << n[0] << " " << n[1] << " " << n[2] << std::endl;
                }
                for (auto t : this->getTriangles())
                {
                    objF << "f " << (t[0] + 1) << "//" << (t[0] + 1) << " " << (t[1] + 1) << "//" << (t[1] + 1) << " " << (t[2] + 1) << "//" << (t[2] + 1) << std::endl;
                }
                objF.close();
            }
        }
        void writeObjQuads(std::string& o_file)
        {
            std::ofstream objF;
            objF.open(o_file.c_str());
            if (!objF.is_open()) {
                std::cout << "ERROR: can't open output file " << std::endl;
            }
            else {
                objF << "#Dual Marching Cubes\n";
                const int nr_v = this->nrVertices();
                for (int i = 0; i < nr_v; i++)
                {
                    Mesh::Vertex v = this->vertex(i);
                    Mesh::Attribute a = this->attribute(i);
                    objF << "v " << v[0] << " " << v[1] << " " << v[2] << " " << a[0] << " " << a[1] << " " << a[2] << std::endl;
                }
                for (auto n : this->getNormals())
                {
                    objF << "vn " << n[0] << " " << n[1] << " " << n[2] << std::endl;
                }
                for (auto q : this->getQuadrilaterals())
                {
                    objF << "f " << (q[0] + 1) << "//" << (q[0] + 1) << " " << (q[1] + 1) << "//" << (q[1] + 1) << " " << (q[2] + 1) << "//" << (q[2] + 1) << " " << (q[3] + 1) << "//" << (q[3] + 1) << std::endl;
                }
                objF.close();
            }
        }
	private:
		std::vector<Vertex> m_vertices;
		std::vector<Normal> m_normals;
		std::vector<Triangle> m_triangles;
		std::vector<Quadrilateral> m_quadrilaterals;
		std::vector<Attribute> m_attributes;

		std::vector<Line> m_lines;
		std::vector<Vertex> m_lineVertices;
		std::vector<Normal> m_lineNormals;
		std::vector<Attribute> m_lineAttributes;
		BBox m_bbox;
	};
}